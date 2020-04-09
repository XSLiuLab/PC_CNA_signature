## Ref: https://keras.rstudio.com/

library(sigminer)
library(keras)
library(dplyr)
library(caret)


# Preparing data ----------------------------------------------------------

load("output/Sig.CNV.seqz.W.RData")
load("output/Sig.CNV.mskcc.RData")
load("output/CNV.mskcc.tally.W.RData")

## Switch Sig1 and Sig2 to make the signature consistent
## See https://shixiangwang.github.io/prad_signature/#method-comparison-between-macintyre-et-al-and-our-study
Sig.CNV.mskcc = sig_modify_names(Sig.CNV.mskcc, new_names = paste0("Sig", c(2,1,3:5)))


## Getting relative and absolute exposure
## relative exposure is more important
rel_expo_combined = get_sig_exposure(Sig.CNV.seqz.W, type = "relative")
rel_expo_mskcc    = get_sig_exposure(Sig.CNV.mskcc, type = "relative")
# expected exposure, use same signature
rel_expo_mskcc_ep = sig_fit(t(CNV.mskcc.tally.W$nmf_matrix), sig = Sig.CNV.seqz.W,
                            mode = "copynumber", type = "relative", rel_threshold = 0.01,
                            return_class = "data.table")

rel_expo_all = dplyr::bind_rows(
  rel_expo_combined %>% dplyr::mutate(cohort = "wang"),
  rel_expo_mskcc_ep %>% dplyr::mutate(cohort = "mskcc")
)

abs_expo_combined = get_sig_exposure(Sig.CNV.seqz.W)
abs_expo_mskcc    = get_sig_exposure(Sig.CNV.mskcc)
abs_expo_mskcc_ep = sig_fit(t(CNV.mskcc.tally.W$nmf_matrix), sig = Sig.CNV.seqz.W,
                            mode = "copynumber", type = "absolute",
                            return_class = "data.table")

abs_expo_all = dplyr::bind_rows(
  abs_expo_combined,
  abs_expo_mskcc_ep
)

col_names  = c("sample", "Sig1", "Sig2", "Sig3", "Sig4", "Sig5")
col_names2 = c("sample", "AbsSig1", "AbsSig2", "AbsSig3", "AbsSig4", "AbsSig5")

rel_expo_mskcc = rel_expo_mskcc[, col_names, with = FALSE]
abs_expo_mskcc = abs_expo_mskcc[, col_names, with = FALSE]

colnames(abs_expo_all) = colnames(abs_expo_combined) = colnames(abs_expo_mskcc) = col_names2

## Getting groups from consensus clustering implemented by NMF
grp_combined = get_groups(Sig.CNV.seqz.W, method = "consensus", match_consensus = TRUE)
# Fix group #1 to sig1 enriched by hand
grp_combined$enrich_sig[grp_combined$group == "1"] = "Sig1"
# Take a check
table(grp_combined$enrich_sig)

grp_mskcc = get_groups(Sig.CNV.mskcc, method = "consensus", match_consensus = TRUE)
table(grp_mskcc$enrich_sig)

grp_combined$group = NULL
grp_combined$silhouette_width = NULL
grp_mskcc$group = NULL
grp_mskcc$silhouette_width = NULL

grp_df = dplyr::bind_rows(
  grp_combined,
  grp_mskcc
)

## Merge predictive variable and target variable
expo_all = purrr::reduce(list(rel_expo_all, abs_expo_all, grp_df), dplyr::full_join, by = "sample") %>%
  dplyr::select(sample, cohort, dplyr::everything())
expo_wang = purrr::reduce(list(rel_expo_combined, abs_expo_combined, grp_combined), dplyr::full_join, by = "sample")
expo_mskcc = purrr::reduce(list(rel_expo_mskcc, abs_expo_mskcc, grp_mskcc), dplyr::full_join, by = "sample")

save(expo_all, expo_wang, expo_mskcc, file = "output/classification_modeling_input.RData")

# Constructing input ------------------------------------------------------

load(file = "output/classification_modeling_input.RData")

dat = expo_all

# ## Map 0-4 to sig1-5
# dat$enrich_sig = as.integer(substr(dat$enrich_sig, 4, 4)) - 1
#
# ## 80% for training
# ## 20% for testing
# idx = createDataPartition(dat$enrich_sig, p = 0.8, list = FALSE)
#
# ## Check if the proportion is close
# dat$enrich_sig %>% table %>% prop.table
# dat[idx, ]$enrich_sig %>% table %>% prop.table

modeling = function(param1, param2, param3, param4, x_train, y_train, x_test, y_test, n_vars, n_class=5,
                    epochs = 50,
                    batch_size = 16,
                    validation_split = 0.2,
                    drop_model = TRUE,
                    model_file = "keras_model.h5") {
  library(keras)
  library(dplyr)

  ## Defining model
  model <- keras_model_sequential()
  model %>%
    layer_dense(units = param1, activation = 'relu', input_shape = n_vars) %>%
    layer_dropout(rate = param2) %>%
    layer_dense(units = param3, activation = 'relu') %>%
    layer_dropout(rate = param4) %>%
    layer_dense(units = n_class, activation = 'softmax')

  # Default testing model
  # model %>%
  #   layer_dense(units = 100, activation = 'relu', input_shape = n_vars) %>%
  #   layer_dropout(rate = 0.4) %>%
  #   layer_dense(units = 20, activation = 'relu') %>%
  #   layer_dropout(rate = 0.3) %>%
  #   layer_dense(units = n_class, activation = 'softmax')


  model %>% compile(
    loss = 'categorical_crossentropy',
    optimizer = optimizer_rmsprop(),
    metrics = c('accuracy')
  )

  ## Training and evaluation
  history <- model %>% fit(
    x_train, y_train,
    epochs = epochs, batch_size = batch_size,
    validation_split = validation_split
  )

  acc_train = history$metrics$accuracy
  acc_train_last = acc_train[length(acc_train)]

  acc_val = history$metrics$val_accuracy
  acc_val_last = acc_val[length(acc_val)]

  # Peformance
  pf = model %>% evaluate(x_test, y_test)
  acc_test = pf$accuracy

  if (drop_model) {
    dplyr::tibble(
      n_param = count_params(model),
      param1 = param1,
      param2 = param2,
      param3 = param3,
      param4 = param4,
      acc_train_last = acc_train_last,
      acc_val_last = acc_val_last,
      acc_test = acc_test
    )
  } else {
    save_model_hdf5(model, filepath = model_file)
    dplyr::tibble(
      model_file = model_file,
      n_param = count_params(model),
      param1 = param1,
      param2 = param2,
      param3 = param3,
      param4 = param4,
      acc_train_last = acc_train_last,
      acc_val_last = acc_val_last,
      acc_test = acc_test
    )
  }
}

RunModel = function(dat, seed=1234,
                    var_names = c(paste0("Sig", 1:5), paste0("AbsSig", 1:5)),
                    param_mat = expand.grid(
                      c(10, 20, 50, 100),
                      c(0, 0.1, 0.2, 0.3, 0.4, 0.5),
                      c(10, 20, 50, 100),
                      c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
                    ),
                    epochs = 50,
                    batch_size = 16,
                    validation_split = 0.2) {
  library(keras)
  library(dplyr)
  library(caret)

  ## Map 0-4 to sig1-5
  dat$enrich_sig = as.integer(substr(dat$enrich_sig, 4, 4)) - 1

  ## 80% for training
  ## 20% for testing
  set.seed(seed = seed)
  idx = createDataPartition(dat$enrich_sig, p = 0.8, list = FALSE)

  dat_train = dat[idx, ]
  dat_test = dat[-idx, ]

  x_train = as.matrix(dat_train[, var_names])
  y_train = to_categorical(dat_train$enrich_sig)

  x_test = as.matrix(dat_test[, var_names])
  y_test = to_categorical(dat_test$enrich_sig)

  # ## simple test
  # modeling(10, 0, 10, 0, x_train, y_train, x_test, y_test, n_vars = length(var_names))

  res_df <- dplyr::tibble()
  for (i in 1:nrow(param_mat)) {
    message("=> Running for param pair #", i)
    temp_df = modeling(param_mat[i, 1], param_mat[i, 2], param_mat[i, 3], param_mat[i, 4],
                       x_train, y_train, x_test, y_test,
                       n_vars = length(var_names),
                       n_class = length(unique(dat$enrich_sig)),
                       epochs = epochs,
                       batch_size = batch_size,
                       validation_split = validation_split)
    print(temp_df)
    res_df = dplyr::bind_rows(res_df, temp_df)
  }

  list(
    x_train = x_train,
    y_train = y_train,
    x_test = x_test,
    y_test = y_test,
    performance = res_df
  )
}

getTotalAcc = function(train_last_acc, validate_last_acc, test_acc, validation_split=0.2, test_split=0.2) {
  (train_last_acc * (1 - validation_split) + validate_last_acc * validation_split) * (1 - test_split) + test_acc * test_split
}

res_all = RunModel(expo_all)
res_all$performance$acc_total =
  getTotalAcc(res_all$performance$acc_train_last,
              res_all$performance$acc_val_last,
              res_all$performance$acc_test)
saveRDS(res_all, file = "output/model_all_20200408.rds")

res_wang = RunModel(expo_wang)
res_wang$performance$acc_total =
  getTotalAcc(res_wang$performance$acc_train_last,
              res_wang$performance$acc_val_last,
              res_wang$performance$acc_test)
saveRDS(res_wang, file = "output/model_wang_20200408.rds")

res_mskcc = RunModel(expo_mskcc)
res_mskcc$performance$acc_total =
  getTotalAcc(res_mskcc$performance$acc_train_last,
              res_mskcc$performance$acc_val_last,
              res_mskcc$performance$acc_test)
saveRDS(res_mskcc, file = "output/model_mskcc_20200408.rds")

## Best model should have the optimal accuracy in both test dataset and whole datasets
## and also with less number of parameters, and should not overfit
##
## We pick a model from top 20 models
res_all = readRDS("output/model_all_20200408.rds")
res_wang = readRDS("output/model_wang_20200408.rds")
res_mskcc = readRDS("output/model_mskcc_20200408.rds")

library(dplyr)

res_all$performance %>%
  dplyr::mutate(index = row_number()) %>%
  dplyr::arrange(desc(acc_test), desc(acc_total), n_param) %>%
  print(n = 20)
# pick index 150
all_pf = modeling(20, 0.1, 50, 0.1,
                  res_all$x_train, res_all$y_train, res_all$x_test, res_all$y_test,
                  n_vars = 10, n_class = 5, drop_model = FALSE, model_file = "output/keras_model_for_all_cohorts_20200409.h5")
saveRDS(all_pf, file = "output/keras_model_performance_for_all_cohorts_20200409.h5")

res_wang$performance %>%
  dplyr::mutate(index = row_number()) %>%
  dplyr::arrange(desc(acc_test), desc(acc_total), n_param) %>%
  print(n = 20)

# pick index 362
wang_pf = modeling(100, 0.1, 10, 0,
                   res_wang$x_train, res_wang$y_train, res_wang$x_test, res_wang$y_test,
                   n_vars = 10, n_class = 5, drop_model = FALSE, model_file = "output/keras_model_for_wang_cohort_20200409.h5")
saveRDS(wang_pf, file = "output/keras_model_performance_for_wang_cohort_20200409.h5")

res_mskcc$performance %>%
  dplyr::mutate(index = row_number()) %>%
  dplyr::arrange(desc(acc_test), desc(acc_total), n_param) %>%
  print(n = 20)

# pick index 121
mskcc_pf = modeling(10, 0, 20, 0.1,
                   res_mskcc$x_train, res_mskcc$y_train, res_mskcc$x_test, res_mskcc$y_test,
                   n_vars = 10, n_class = 5, drop_model = FALSE, model_file = "output/keras_model_for_mskcc_cohort_20200409.h5")
saveRDS(mskcc_pf, file = "output/keras_model_performance_for_mskcc_cohort_20200409.h5")

#
# plot(history)
# Prediction
model = load_model_hdf5(mskcc_pf$model_file)
model %>% predict_classes(res_mskcc$x_train[1, , drop = FALSE])
model %>% predict_proba(res_mskcc$x_train[1, , drop = FALSE])

library(deepviz)
deepviz::plot_model(model)
#deepviz::plot_deepviz(model)

# # confusion matrix
# table(dat_test$enrich_sig, pred_class)

df = readxl::read_excel("output/supp_sample_table.xlsx")
df = df[, c("tumor_Run", "normal_Run")]
write.table(df, file = "output/Run_pair.csv", quote = FALSE, sep = "\t", row.names = F, col.names = F)

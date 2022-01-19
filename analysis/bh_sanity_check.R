#checking BH approach from python with BH R builtin
library(data.table)
df = fread('/oak/stanford/groups/horence/rob/isoform_localizations/SRRS/outputs/gene/Merfish_MOp_peripheral.csv')

#p column is not multiple hypothesis corrected
df$r_bh_p = p.adjust(df$p, "BH")

#bh_p column is from python BH calculation using multitest package
#r_bh_p columns is using the R builtin from p.adjust
print("Correlation between R-calculated and python-calculated BH correct p-values")
print(cor(df$r_bh_p,df$bh_p))


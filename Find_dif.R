library(amplicon)
otu = read.csv("species.csv", row.names = 1)
meta = read.csv("metadata.csv", row.names = 1)

output = compare(data = otu, metadata = meta,
                 group = "Group", compare_pair = "Day0-Control",
                 method = "wilcox", RA = 0.001,
                 pvalue = 0.05, fdr = 0.05)
write.csv(output,"dif_at_baseline.csv")
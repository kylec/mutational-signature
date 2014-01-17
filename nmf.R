#
library(NMF)
#fin="originalGenomes.txt"
fin="alloriginalGenomes.txt"
#fin="prim-originalGenomes.txt"
#fin="recur-originalGenomes.txt"
a= read.table(fin, sep="\t")
aa = as.matrix(a)
#names= system("ls TCGA*.txt | sed 's/.txt//' | cut -d- -f2-4", intern=T)
names = read.table("nmf-pat-names.txt")
colnames(aa)=names[,1]

# estimate ranks
#png("est-10.png")
#png("prim-est-10.png")
png("recur-est-10.png")
est10=nmfEstimateRank(aa, 2:10, nrun=10, "lee")
plot(est10)
dev.off()

est50=nmfEstimateRank(aa, 2:10, nrun=50, "lee")
plot(est50)
dev.off()

# nmf from 2 to 10
nmf10 = nmf(aa,2:10,nrun=10,"lee")
consensusmap(nmf10)

nmf3 = nmf(aa,3,nrun=10,"lee")
consensusmap(nmf3)

nmf4 = nmf(aa,4,nrun=50,"lee")
png("nmf4",width="1200",height="1200")
consensusmap(nmf4)
dev.off()

# extract features
extractFeatures(nmf2)

# extract basis
apply(basis(nmfb5), 1, max) == basis(nmfb5)[,1] 



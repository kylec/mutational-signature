cl=read.table("OV.clin.merged.txt", sep="\t", quote="", header=T, fill=T, stringsAsFactors = FALSE)

cl = cl[,c(2:ncol(cl))]
#rownames(cl) = colnames(cl)
tcl = data.frame(t(cl))
colnames(tcl) = seq(1:ncol(tcl))

sam=read.table("/home/kchang3/data/clinical/clin-patients.txt", stringsAsFactors = FALSE)

s=c("tcga-13-0791", "tcga-13-0913", "tcga-13-1489", "tcga-13-1817", "tcga-13-1819", "tcga-24-1852", "tcga-29-1692", "tcga-29-1704", "tcga-29-1705", "tcga-29-1707", "tcga-29-1710", "tcga-29-1770", "tcga-29-2414", "tcga-61-1916", "tcga-61-2008", "tcga-61-2095")


library(survival)
a = read.table("wgs-we-clinical.txt", header=T)
surv = Surv(a$time, a$status)
mysurvfit <- survfit(Surv(a$time, a$status)~a$sig, data=a)
colors=c("black", "red", "blue", "green" )
plot(mysurvfit, col=colors, main="ova triplets exome")
sigs=c("sig1", "sig2", "sig3", "sig4")
legend("bottomright", sigs, fill=colors)
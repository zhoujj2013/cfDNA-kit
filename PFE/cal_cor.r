args <- commandArgs(TRUE)

if (length(args) < 1) {
	stop("Please input enough args")
}

f = args[1]

dat = read.table(file=f, header=F)

r = cor.test(dat$V1,dat$V2, method="spearman", exact=FALSE)
rr = c(as.numeric(r$estimate), r$p.value)
print(rr)

#names(r$estimate)
#write.table(rr,file="result.txt", sep="\t")
#write.table(rr,file="result.txt", sep="\t", header=F)
#write.table(rr,file="result.txt", sep="\t", header=F)
#write.table(rr,file="result.txt", sep="\t", header="F")
#write.table(rr,file="result.txt", sep="\t", row.names=F)
#write.table(rr,file="result.txt", sep="\t", row.names=F, col.names=F)
#savehistory()

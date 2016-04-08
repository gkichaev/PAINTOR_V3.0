files <- list.files(pattern=".annotations")
freq <- .2

for(f in files){
	input <- read.table(f, header=T)
	rand.annot <- rbinom(nrow(input),1,.2)
	out_frame <- data.frame(input, annot_null=rand.annot)
	write.table(out_frame, paste0(f,".null"), row.names=F, quote=F)
}

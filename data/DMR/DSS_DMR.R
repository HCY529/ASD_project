data.1 = read.table("ASDF_1_trimmed_bismark_bt2.dss.input.txt",header = T)
data.2 = read.table("ASDF_2_trimmed_bismark_bt2.dss.input.txt",header = T)
data.3 = read.table("ASDM_1_trimmed_bismark_bt2.dss.input.txt",header = T)
data.4 = read.table("ASDM_2_trimmed_bismark_bt2.dss.input.txt",header = T)
data.5 = read.table("CTF_1_trimmed_bismark_bt2.dss.input.txt",header = T)
data.6 = read.table("CTF_2_trimmed_bismark_bt2.dss.input.txt",header = T)
data.7 = read.table("CTM_1_trimmed_bismark_bt2.dss.input.txt",header = T)
data.8 = read.table("CTM_2_trimmed_bismark_bt2.dss.input.txt",header = T)
BSobj =  makeBSseqData(list(data.1,data.2,data.3,data.4,data.5,data.6,data.7,date8),c("N1","N2","N3","N4","C1","C2","C3","C4"))[1:1000,]
library(DSS)
BSobj =  makeBSseqData(list(data.1,data.2,data.3,data.4,data.5,data.6,data.7,date8),c("N1","N2","N3","N4","C1","C2","C3","C4"))[1:1000,]
BSobj =  makeBSseqData(list(data.1,data.2,data.3,data.4,data.5,data.6,data.7,date.8),c("N1","N2","N3","N4","C1","C2","C3","C4"))[1:1000,]
BSobj =  makeBSseqData(list(data.1,data.2,data.3,data.4,data.5,data.6,data.7,data.8),c("N1","N2","N3","N4","C1","C2","C3","C4"))[1:1000,]
dmlTest <- DMLtest(BSobj, group1=c("N1","N2","N3","N4"), group2=c("C1","C2","C3","C4"),smoothing=TRUE)
dmls <- callDML(dmlTest.sm, p.threshold=0.001)
head(dmlTest.sm)
dmrs <- callDMR(dmlTest, delta=0.1, p.threshold=0.05)
showOneDMR(dmrs[1,], BSobj)

dmlTest <- DMLtest(BSobj, group1=c("N1","N2","N3","N4"), group2=c("C1","C2","C3","C4"),smoothing=TRUE)
head(dmlTest)
 dmrs <- callDMR(dmlTest, delta=0.1, p.threshold=0.001)
 dmrs <- callDMR(dmlTest, delta=0.1, p.threshold=0.1)
showOneDMR(dmrs[1,], BSobj)
write.table(dmrs,file=paste0("dmrs.txt"),sep = "\t",append = FALSE,row.names = FALSE, col.names = TRUE, quote = FALSE)
write.bed <- function(dmrs,file){
    methy <- data.frame(chrom=unlist(lapply(dmrs$chr, function(x){
        t <- strsplit(x,split = "chr")[[1]][2]
        return(t)
    })),chromStart=dmrs$start,chromEnd=dmrs$end,name=1:nrow(dmrs),score=dmrs$nCG,strand=".",thickStart=dmrs$start,thickEnd=dmrs$end,itemRGB="255,255,255",blockCount=1,blockSizes=dmrs$end-dmrs$start,blockStarts=dmrs$start,signalValue=dmrs$diff.Methy,pValue=-1,qValue=-1)
    write.table(methy,file = file,sep = "\t",row.names = F,col.names = F,quote = F)
}
write.bed(dmrs,paste0(out_dir,"dmrs.bed"))
dmrs_bed<- readPeakFile(paste0(out_dir,"dmrs.bed"))
write.bed(dmrs,paste0("dmrs.bed"))
dmrs_bed<- readPeakFile(paste0("dmrs.bed"))
ibrary(ChIPseeker)
library(ChIPseeker)
write.bed(dmrs,paste0("dmrs.bed"))
dmrs_bed<- readPeakFile(paste0("dmrs.bed"))
head dmrs
write.bed(dmrs,paste0("dmrs.bed"))
dmrs_bed<- readPeakFile(paste0("dmrs.bed"))
write.bed(dmrs,paste0("dmrs.bed"))
dmrs_bed<- readPeakFile(paste0("dmrs.bed"))
less dmrs.txt
showOneDMR(dmrs[2,], BSobj)
write.bed <- function(dmrs,file){
    methy <- data.frame(chrom=unlist(lapply(dmrs$chr, function(x){
        t <- x
        return(t)
    })),chromStart=dmrs$start,chromEnd=dmrs$end,name=1:nrow(dmrs),score=dmrs$nCG,strand=".",thickStart=dmrs$start,thickEnd=dmrs$end,itemRGB="255,255,255",blockCount=1,blockSizes=dmrs$end-dmrs$start,blockStarts=dmrs$start,signalValue=dmrs$diff.Methy,pValue=-1,qValue=-1)
    write.table(methy,file = file,sep = "\t",row.names = F,col.names = F,quote = F)
}
write.bed(dmrs,paste0("dmrs.bed"))
dmrs_bed<- readPeakFile(paste0("dmrs.bed"))
covplot(dmrs_bed,title="dmrs")
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPseeker)
BiocManager::install("GenomicFeatures")
BiocManager::install("GenomicFeatures")
BiocManager::install("GenomicFeatures",force=TRUE)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
BiocManager::install("GenomicFeatures",force=TRUE)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene",force=TRUE)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicFeatures)
upload(GenomicFeatures)
unload(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
q()

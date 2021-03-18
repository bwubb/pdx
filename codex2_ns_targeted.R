library(CODEX2)
##########################################################
# Initialization
##########################################################

if (exists("snakemake")) {
    all_bams<-as.character(snakemake@input)
    head(snakemake@input)
    projectname<-snakemake@params[['project']] #should be dir and not a prefix
    bedFile<-snakemake@params[['bed']]
    #print("YO",snakemake@params[['bed']])
    normal_samples<-read.table(file=snakemake@params[['normals']])
}else{
    args = commandArgs(trailingOnly=TRUE)
    print("No snakmake found")
    projectname<-paste("data","work",args[1],"codex2",sep="/")
    bedFile<-args[2]
    all_bams<-as.character(read.table(file=args[3])$V2)
    normal_samples<-as.character(read.table(file=args[4])$V1)
}
print("Project")
print(projectname)
print("Bed file")
print(bedFile)
print("Bam Files")
head(all_bams)


sampname<-gsub('.ready.bam','',basename(all_bams))
bambedObj <- getbambed(bamdir=all_bams,bedFile=bedFile,sampname=sampname,projectname=projectname)
bamdir <- bambedObj$bamdir; sampname <- bambedObj$sampname
ref <- bambedObj$ref; projectname <- bambedObj$projectname

##########################################################
# Getting GC content and mappability
##########################################################
gc <- getgc(ref)
mapp <- getmapp(ref)

##########################################################
# Getting gene names, needed for targeted sequencing, here generating gene names in silico
##########################################################
#gene=rep(NA,length(ref))
#for(chr in as.matrix(unique(seqnames(ref)))){
#  chr.index=which(seqnames(ref)==chr)
#  gene[chr.index]=paste(chr,'_gene_',ceiling(chr.index/30),sep='')
#}

gene<-read.table(bedFile,header=F)$V4
values(ref) <- cbind(values(ref), DataFrame(gc, mapp, gene))

##########################################################
# Getting depth of coverage
##########################################################
coverageObj <- getcoverage(bambedObj, mapqthres = 20)
Y <- coverageObj$Y
write.csv(Y, file = paste(projectname, 'coverage.csv', sep='/'), quote = FALSE)
head(Y[,1:5])

##########################################################
# Quality control
##########################################################
qcObj <- qc(Y, sampname, ref, cov_thresh = c(20, Inf),length_thresh = c(20, Inf), mapp_thresh = 0.9,gc_thresh = c(20, 80))
Y_qc <- qcObj$Y_qc; sampname_qc <- qcObj$sampname_qc
ref_qc <- qcObj$ref_qc; qcmat <- qcObj$qcmat; gc_qc <- ref_qc$gc
write.table(qcmat, file = paste(projectname, 'qcmat.txt', sep='/'),sep = '\t', quote = FALSE, row.names = FALSE)

##########################################################
# Estimating library size factor for each sample
##########################################################
Y.nonzero <- Y_qc[apply(Y_qc, 1, function(x){!any(x==0)}),]
pseudo.sample <- apply(Y.nonzero,1,function(x){exp(1/length(x)*sum(log(x)))})
N <- apply(apply(Y.nonzero, 2, function(x){x/pseudo.sample}), 2, median)
#how do I plot this?
#plot(N, apply(Y_qc,2,sum), xlab='Estimated library size factor', ylab='Total sum of reads')
#
##########################################################
# Genome-wide normalization using normalize_null
##########################################################
# If there are negative control samples, use normalize_codex2_ns()
# If there are negative control regions, use normalize_codex2_nr()


norm_index<-which(colnames(Y_qc) %in% normal_samples)
normObj.null <- normalize_codex2_ns(Y_qc = Y_qc,gc_qc = gc_qc,K = 1:length(norm_index), norm_index=norm_index, N = N)
Yhat <- normObj.null$Yhat
AIC <- normObj.null$AIC; BIC <- normObj.null$BIC
RSS <- normObj.null$RSS

##########################################################
# Number of latent factors
##########################################################
choiceofK(AIC, BIC, RSS, K = 1:length(norm_index), filename = paste(projectname, 'codex2_null_choiceofK.pdf', sep='/'))

#Causes error, but is shouldnt be a big deal. But it is because snakemake

##########################################################
# CBS segmentation per gene: optinmal for targeted seq
##########################################################
source('segment_targeted.R')
# Available at: https://github.com/yuchaojiang/CODEX2/blob/master/targeted_sequencing/segment_targeted.R
optK=which.max(BIC)
finalcall=matrix(ncol=14,nrow=0)
colnames(finalcall)=c('sample_name','chr','gene','cnv',
                      'st_bp','ed_bp','length_kb',
                      'st_exon','ed_exon','raw_cov',
                      'norm_cov','copy_no','lratio',
                      'mBIC')
for(genei in unique(ref_qc$gene)){
  cat('Segmenting gene',genei,'\n')
  geneindex=which(ref_qc$gene==genei)
  yi=Y_qc[geneindex,, drop=FALSE]
  yhati=Yhat[[optK]][geneindex,, drop=FALSE]
  refi=ref_qc[geneindex]
  finalcalli=segment_targeted(yi, yhati, sampname_qc, refi, genei, lmax=length(geneindex), mode='fraction') 
  finalcall=rbind(finalcall,finalcalli)
}
cn=(as.numeric(as.matrix(finalcall[,'copy_no'])))
cn.filter=(cn<=1.7)|(cn>=2.3) # removing calls with fractional copy numbers close to 2 (for heterogeneous cancer samples)
finalcall=finalcall[cn.filter,]
length_exon=as.numeric(finalcall[,'ed_exon'])-as.numeric(finalcall[,'st_exon'])+1
finalcall=cbind(finalcall[,1:7],length_exon,finalcall[,10:14])

write.table(finalcall, file = paste(projectname,'segments.txt',sep='/'),sep='\t',quote=FALSE,row.names=FALSE)

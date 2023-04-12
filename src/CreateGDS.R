library(data.table)
library(parallel)
library(SNPRelate)
library(EthSEQ)

genoVCF = fread(snakemake@input[['vcf']],data.table = F)
# genoVCF = fread("./out/vcf/head",data.table = F)
sif = fread(snakemake@input[['sif']],data.table = F)
# sif=fread('./data/sif.tsv',data.table = F)
subset = fread(snakemake@input[['header']],header=F,data.table = F)
# subset = fread('./out/all_hg38_ext_header.vcf',data.table = F,header = F)
ncores = snakemake@threads

gdsFile = snakemake@output[['gds']]

sif = sif[match(subset[-(1:9)],sif$Sample),]

# library(GenomicRanges)
# variants = makeGRangesFromDataFrame(vcf,seqnames.field = "V1",start.field = "V2",end.field = "V2",keep.extra.columns = T)
# kits = fread("/shares/CIBIO-Storage/BCGLAB/EthSEQ_salmon/all_hg38_ext.bed")
# kits$V1 = gsub("chr","",kits$V1)
# kitsRanges = makeGRangesFromDataFrame(kits,seqnames.field = "V1",start.field = "V2",end.field = "V3")
# olaps = findOverlaps(variants,kitsRanges)
# variants_sub = variants[queryHits(olaps),]
# idx = match(variants_sub$V3,map[,2])
# map = map[idx,]
# vcf = vcf[idx,]
# ped = ped[,c(1:6,idx+6)]
# 
# sif = read.csv("/shares/CIBIO-Storage/BCGLAB/Resources/1000GP/igsr_samples.tsv",sep="\t",as.is=T)
# sif=sif[sif$Sample%in%rownames(ped),]
# 
# ped[1:10,1:10]
# cat(all(vcf[,1]==map[,1]),"\n")
# cat(all(vcf[,2]==map[,4]),"\n")
# cat(all(vcf[,3]==map[,2]),"\n")

res = mclapply(10:ncol(genoVCF),function(i)
{
  tmp = genoVCF[,i]
  tmp[which(tmp=="1|1")] = 2
  tmp[which(tmp=="0|1")] = 1
  tmp[which(tmp=="1|0")] = 1
  tmp[which(tmp=="0|0")] = 0
  tmp[which(!tmp%in%c(0,1,2))] = 3
  return(as.integer(tmp))
},mc.cores=ncores)
geno = do.call(cbind,res)
# geno = matrix(unlist(res),nrow = length(res[[1]]),byrow = F)
# ped = cbind(ped[,1:6],geno)
# samples = intersect(sif$Sample,ped[,2])

# ped = ped[which(ped[,2]%in%samples),]
# sif = sif[which(sif$Sample%in%samples),]
# all(sif$Sample==ped[,2])

# idx = order(vcf$V1,vcf$V2)
# vcf = vcf[idx,]
# map = map[idx,]
# ped = ped[,c(1:6,idx+6)]
# all(vcf$V2==map[,4])

#write(sign[idx],paste("/elaborazioni/sharedCO/Precision_Medicine/Ethnicity/TAGC_sharing/",labout,".positions",sep=""))
#fwrite(ped,paste("/elaborazioni/sharedCO/CO_Shares/Code/ethseq/Models/",labout,".Model.ped",sep=""),sep="\t",quote = F,col.names = F,row.names = F,nThread=5)
#fwrite(map,paste("/elaborazioni/sharedCO/CO_Shares/Code/ethseq/Models/",labout,".Model.map",sep=""),sep="\t",quote = F,col.names = F,row.names = F,nThread=5)
#fwrite(vcf,paste("/elaborazioni/sharedCO/CO_Shares/Code/ethseq/VCF/",labout,".Model.vcf",sep=""),sep="\t",quote = F,col.names = F,row.names = F,nThread=5)
#fwrite(sif,paste("/elaborazioni/sharedCO/CO_Shares/Code/ethseq/Models/",labout,".Model.txt",sep=""),sep="\t",quote = F,col.names = T,row.names = F,nThread=5)

# geno = ped[,7:ncol(ped)]
# res = mclapply(1:ncol(geno),function(i)
# {
#   tmp = geno[,i]
#   tmp[which(!tmp%in%c("A A","A B","B B"))] = "3"
#   tmp[which(tmp=="A A")] = "0"
#   tmp[which(tmp=="A B")] = "1"
#   tmp[which(tmp=="B B")] = "2"
#   return(as.numeric(tmp))
# },mc.cores=20)
# 
# geno = matrix(as.numeric(unlist(res)),nrow = length(res[[1]]),byrow = F)

mafs = 1-apply(geno,1,function(x) (length(which(x==0))*2+length(which(x==1)))/(length(which(x!=3))*2))
idx = which(mafs>0.5)
for(i in idx)
{
  tmp = geno[i,]
  tmp[which(tmp==0)] = 4
  tmp[which(tmp==2)] = 0
  tmp[which(tmp==4)] = 2
  geno[i,] = tmp
}


snp.allele = rep("A/B",nrow(geno))
snp.allele[idx] = "B/A"

geno = as.matrix(t(geno))
##########
####################moved at the bottom
##########

# genofile <- snpgdsOpen(paste("/shares/CIBIO-Storage/BCGLAB/EthSEQ_salmon/VCF_model/Exonic_hg38.Model4.gds",sep=""),readonly = F)

# Add your sample annotation
#x = read.gdsn(index.gdsn(genofile, "sample.annot"))
sex = sif$Gender
# unique(sex)
sex[which(sex%in%c(1,'male'))] = "M"
sex[which(sex%in%c(2,'female'))] = "F"
samp.annot <- data.frame(pop.group = sif$Pop,sex = sex)
# add.gdsn(genofile, "sample.annot", samp.annot,replace=T)

# sign = paste(read.gdsn(index.gdsn(genofile, "snp.rs.id")),read.gdsn(index.gdsn(genofile, "snp.chromosome")),read.gdsn(index.gdsn(genofile, "snp.position")),sep="-")
# sign.vcf = paste(vcf$V3,vcf$V1,vcf$V2,sep="-")
# isort = match(sign,sign.vcf)
# vcf = vcf[isort,]
# sign.vcf = paste(vcf$V3,vcf$V1,vcf$V2,sep="-")
# all(sign.vcf==sign)
# 
# add.gdsn(genofile, "snp.ref", vcf$V4)
# add.gdsn(genofile, "snp.alt", vcf$V5)


snpgdsCreateGeno(gdsFile,
                 genmat = geno,
                 sample.id = sif$Sample,
                 snp.id = 1:ncol(geno),
                 snp.rs.id = genoVCF[,3],
                 snp.chromosome = genoVCF[,1],
                 snp.position = genoVCF[,2],
                 snp.allele = snp.allele,
                 other.vars = list("sample.annot"=samp.annot,"snp.ref"=genoVCF[,4],"snp.alt"=genoVCF[,5]),
                 snpfirstdim=FALSE)

# snpgdsClose(genofile)

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

# Add your sample annotation
sex = sif$Gender
sex[which(sex%in%c(1,'male'))] = "M"
sex[which(sex%in%c(2,'female'))] = "F"
if(length(unique(sif$Pop))==1){
  samp.annot <- data.frame(pop.group = sif$Sub.pop,sex = sex)
} else {
  samp.annot <- data.frame(pop.group = sif$Pop,sex = sex)
}

snpgdsCreateGeno(gdsFile,
                 genmat = t(geno),
                 sample.id = sif$Sample,
                 snp.id = paste(genoVCF[,1],genoVCF[,2],genoVCF[,3],sep=":"),
                 snp.rs.id = genoVCF[,3],
                 snp.chromosome = genoVCF[,1],
                 snp.position = genoVCF[,2],
                 snp.allele = snp.allele,
                 other.vars = list("sample.annot"=samp.annot,"snp.ref"=genoVCF[,4],"snp.alt"=genoVCF[,5]),
                 snpfirstdim=TRUE)

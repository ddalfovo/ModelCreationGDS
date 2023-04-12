import os

input = snakemake.input
output = snakemake.output
params = snakemake.params


vcf = input[1]
subSamples = input[2] 
bed = input[3]
dbsnp = input[4]

out_vcf = output[0]

minAF = params[0]
# minAF = 0.01
with open(subSamples) as file:
    subs = file.readlines()[0].strip().split("\t")

head = os.popen("zgrep -P 'CHROM\tPOS' -m 1 " + vcf)
header = head.readlines()[0].strip().split("\t")

idx = [i for i in range(len(header)) if header[i] in subs]


vcfTabix = os.popen("tabix " + vcf + " -R " + bed + " | awk '{if(length($4)==1 && length($5)==1) print}'")

vcfGeno = []
vcfPosition = []
for ln in vcfTabix:
    vcf_line = ln.strip().split('\t')
    geno = [vcf_line[i] for i in idx]
    if not geno[8]=='GT':
        vcfFormat = geno[8].split(':')
        GTidx = vcfFormat.index("GT")
        geno = geno[:8] + ["GT"] + [i.split(":")[GTidx] for i in geno[9:]]
    vcfGeno.append(geno)
    vcfPosition.append("\t".join(geno[:2] + [geno[1]])+"\n" )

out_pos = out_vcf + ".positions"
outF = open(out_pos,'w')
outF.writelines(vcfPosition)
outF.close()

dbsnpTabix = os.popen("tabix " + dbsnp + " -R " + out_pos)

dbsnpRes = {}
for ln in dbsnpTabix:
    dbsnp_line = ln.strip().split('\t')
    keyPos = dbsnp_line[0] + ":" + dbsnp_line[1]
    keyRef = len(dbsnp_line[3])
    keyAlt = len(dbsnp_line[4])
    if keyRef==1 & keyAlt==1:
        dbsnpRes[keyPos] = dbsnp_line[2]
os.remove(out_pos)


rsIDs = [dbsnpRes[i[0]+":"+i[1]] if i[0]+":"+i[1] in dbsnpRes else "." for i in vcfGeno]
def retAF(info):
    d = dict(x.split("=") for x in info.split(';') if len(x.split("="))==2)
    num = float(d['AF'])
    return minAF<num<(1-minAF)

vcfFinal = []
for i in range(len(vcfGeno)):
    vcfGeno[i][2] = rsIDs[i]
    if retAF(vcfGeno[i][7]):
        vcfFinal.append("\t".join(vcfGeno[i])+'\n')

outFile = open(out_vcf,'w')
outFile.writelines(vcfFinal)
outFile.close()


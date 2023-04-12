from pathlib import Path
import pandas as pd
import os
# --- Variable Declarations ---- #
# --- Importing Configuration Files --- #
configfile: "paths.yaml"
ncores = 20

vcfs = pd.read_csv(config['vcf'],header=None,names=['vcf'])
beds = pd.read_csv(config['bed'],header=None,names=['bed'])
dbsnp = pd.read_csv(config['dbsnp'],header=None,names=['dbsnp'])
subSamples = pd.read_csv(config['subSamples'],header=None,names=['code'])

VCF = set(vcfs['vcf'])
BED = set(beds['bed'])
DBSNP = set(dbsnp['dbsnp'])
POP = set(subSamples['code'])

runR = "Rscript --vanilla"
logAll = "2>&1"

def removeExtension(input):
    return([os.path.splitext(os.path.basename(Path(i)))[0] for i in input])

def extractName(files):
    return([os.path.splitext(os.path.basename(Path(i)).strip('.gz'))[0] for i in files])

dictVCF = {os.path.splitext(os.path.basename(Path(i)).strip('.gz'))[0]:i for i in vcfs['vcf']}

# def getOut():
#     with open(config['subSamples']) as f:
#         return

# --- Rules --- #
include: config["rules"] + "dag.smk"


# rule all:
#     input:
#         expand(config['out'] + config['assembly'] + '/models/{modelName}.gds',modelName=extractName(BED))
rule all:
    input:
        expand(config['out'] + config['assembly'] + '.{popList}/models/{modelName}.gds',modelName=extractName(BED), popList = POP)


rule createModelGDS:
    input:
        vcf = config['out'] + config['assembly'] + ".{popList}/vcf/{modelName}.vcf.gz",
        sif = config['sif'],
        header = config['out'] + config['assembly'] + '.{popList}/{modelName}_header.vcf',
    output:
        gds = config['out'] + config['assembly'] + '.{popList}/models/{modelName}.gds'
    singularity:
        config['singularity']
    threads:
        20
    script:
        config['src'] + 'CreateGDS_T.R'


rule parseBED:
    input:
        'data/bedFiles' + '/{modelName}.bed',
    output:
        temp(config['tmp'] + '{modelName}.bed')
    params:
        extP = config['extendProbes']
    shell:
        "sed 's/^chr//' {input} | awk '{{FS=OFS=\"\\t\"}} {{print($1,$2-{params.extP},$3+{params.extP})}}' OFS='\\t' | grep -P '^\d+\t\d+\t\d+' | bedtools sort -i stdin | bedtools merge -i stdin > {output}"


rule subsetVCFs:
    input:
        script = 'src/subsetVCF.py',
        v = lambda wc: dictVCF[wc.vcf],
        subset = config['out'] + config['assembly'] + '.{popList}/{modelName}_header.vcf',
        bed = config['tmp'] + config['assembly'] + '/{modelName}.bed',
        dbsnp = expand('{dbsnp}', dbsnp=DBSNP)
    output:
        temp(config['tmp'] + "tempVCFfiles/{modelName}.{popList}/{vcf}.filtered.vcf")
    params:
        minAF = config['minAF']
    threads:
        11
    singularity:
        config['singularity']
    script:
        config['src'] + "subsetVCF.py"


rule mergeSubsetVCFs:
    input:
        expand(config['tmp'] + 'tempVCFfiles/{{modelName}}.{{popList}}/{vcf}.filtered.vcf', vcf=extractName(VCF)),
        header = config['out'] + config['assembly'] + '.{popList}/{modelName}_header.vcf',
    output:
        config['out'] + config['assembly'] + ".{popList}/vcf/{modelName}.vcf"
    params:
        config['tmp'] + '/tempVCFfiles/{modelName}.{popList}/*.filtered.vcf'
    shell:
        "ls -v {params} | xargs cat > {output}.temp && cat {input.header} {output}.temp > {output} && rm {output}.temp"

rule compressVCF:
    input:
        config['out'] + config['assembly'] + ".{popList}/vcf/{modelName}.vcf"
    output:
        config['out'] + config['assembly'] + ".{popList}/vcf/{modelName}.vcf.gz"
    shell:
        "gzip {input}"


rule subHeaderVCF:
    input:
        v = vcfs['vcf'][0],
        popFile = 'data/popLists/samples4model{popList}.txt',
    output:
        path = config['out'] + config['assembly'] + '.{popList}/{modelName}_header.vcf'
    run:
        import os
        with open(input['popFile']) as file:
            subs = [line.strip() for line in file]
        head = os.popen("zgrep -P 'CHROM\tPOS' -m 1 " + input['v'])
        header = head.readlines()[0].strip().split("\t")
        idx_name = [header[i] for i in range(9)]
        idx_name = idx_name + [header[i] for i in range(len(header)) if header[i] in subs]
        f = open(output['path'],'w')
        f.write("\t".join(idx_name) + '\n')



# clean_output   : delete all tmp files
rule clean_output:
    shell:
        "rm -rdf tempVCFfiles/*"

# Rule to copy the gds in a specific folder, adding the assembly and pop informations to the file name.
rule rename_gds:
    params:
        outFolder = config['out']
    script:
        config['src'] + "renameGDS.R"
setwd('/shares/CIBIO-Storage/BCGLAB/Tools/GDSmodelCreation')
out = 'out/'
# out = snakemake@params[['outFolder']]

dirs = list.dirs(out,recursive = F,full.names = F)
dirs = dirs[grep("^hg*",dirs)]
for(i in 1:length(dirs)) {
  gdsFiles = list.files(file.path(out,dirs[i],"models"),full.names = T)
  splitted = strsplit(dirs[i],"\\.")[[1]]
  assembly = splitted[1]
  pop = splitted[2]
  for (gds in gdsFiles){
    system(paste0("cp ",gds," ",out,gsub(".gds","",basename(gds)),".",assembly,".",pop,".Model.gds"))
  }
}

setwd('/proj/hsantosj/users/abhattac/ELGAN_GWAS/eQTL')
require(bigsnpr)
require(data.table)
require(MOSTWAS)

tempFolder = paste0('temp',i,'/')
if (dir.exists(tempFolder)){
  system(paste0('rm -r ',tempFolder))
}
if (!dir.exists(tempFolder)){
  dir.create(tempFolder)
}
require(bigsnpr)
require(data.table)
require(MOSTWAS)
geneLocs = fread('ELGANgenelocs.tsv')
geneLocs = subset(geneLocs,chr %in% c(1:22))
int = floor(seq(0,nrow(geneLocs),length.out=25))
geneList = geneLocs$geneid[int[i]+1:int[i+1]]

snpObj = snp_attach(snp_readBed('ELGAN_tot.bed',
                                backingfile = paste0(tempFolder,
                                                     'temp',i,'_SnpObj.bk')))

mediator = fread('ELGAN_mediators.tsv')
exp = fread('ELGAN_exp.tsv')
colnames(exp) = colnames(mediator)
mediator = rbind(mediator,exp)
mediator = mediator[!duplicated(mediator$Mediator),]
medLocs = fread('ELGANmediatorlocs_3col.tsv')
medLocs$right = medLocs$pos + 1
geneLocs = fread('ELGANgenelocs.tsv')
colnames(medLocs) = colnames(geneLocs)
medLocs = rbind(geneLocs,medLocs)
medLocs = medLocs[!duplicated(medLocs$geneid),]
rm(geneLocs)
covariates = fread('ELGAN_covariates.tsv')
qtlTra = fread('ELGAN_dis_eqtl_top.tsv')
qtMed = fread('ELGAN_loc_snp_to_med.tsv')
qtlTra_parts = paste0('Fold',c(1,2,3),'ELGAN_dis_snp_to_med.tsv')
qtMed_parts = paste0('Fold',c(1,2,3),'ELGAN_loc_snp_to_med.tsv')

for (g in geneList){
  
  print(g)
  if (!paste0(g,'.wgt.med.RData') %in% list.files('DePMAModels/')){
    DePMA(geneInt = g,
      snpObj,
      mediator,
      medLocs,
      covariates,
      cisDist = 1e6,
      qtlTra,
      qtMed,
      h2Pcutoff = .1,
      dimNumeric = 20,
      verbose = F,
      seed = 1218,
      sobel = F,
      nperms = 1000,
      k = 5,
      parallel = F,
      parType = 'no',
      prune = F,
      ldThresh = .5,
      cores = 5,
      qtlTra_parts,
      qtMed_parts,
      modelDir = 'DePMAModels/',
      tempFolder = tempFolder,
      R2Cutoff = 0.01)
    fff = paste0(tempFolder,list.files(tempFolder))
    fff = fff[!(fff %in%  paste0(tempFolder,'temp',
                                 i,'_SnpObj.bk',
                                 c('.bk','.rds')))]
    file.remove(fff)
    }
  
}
system(paste0('rm -r ',tempFolder))
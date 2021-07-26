setwd('/proj/hsantosj/users/abhattac/ELGAN_GWAS/eQTL')
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
int = floor(seq(0,nrow(geneLocs),length.out=15))
geneList = geneLocs$geneid[int[i]+1:int[i+1]]

snpObj = snp_attach(snp_readBed('ELGAN_tot.bed',
		backingfile = paste0(tempFolder,'temp',i,'_SnpObj.bk')))
mediator = fread('ELGAN_mediators.tsv')
exp = fread('ELGAN_exp.tsv')
colnames(exp) = colnames(mediator)
mediator = rbind(mediator,exp)
mediator = mediator[!duplicated(mediator$Mediator),]
medLocs = fread('ELGANmediatorlocs_3col.tsv')
medLocs$right = medLocs$pos + 1
colnames(medLocs) = colnames(geneLocs)
medLocs = rbind(geneLocs,medLocs)
medLocs = medLocs[!duplicated(medLocs$geneid),]
rm(geneLocs)
covariates = fread('ELGAN_covariates.tsv')
qtlFull = fread('ELGAN_dis_medqtl_top10.tsv')

for (g in geneList){

        print(g)
	if (!paste0(g,'.wgt.med.RData') %in% list.files('MeTWASModels/')){
        MeTWAS(geneInt = g,
               snpObj = snpObj,
               mediator = mediator,
               medLocs = medLocs,
               covariates = covariates,
               dimNumeric = 20,
               qtlFull = qtlFull,
               h2Pcutoff = .1,
               numMed = 10,
               seed = 1218,
               k = 5,
               cisDist = 1e6,
               parallel = F,
               prune = F,
               ldThresh = .5,
               cores = 5,
               verbose = F,
               R2Cutoff = .01,
               modelDir = 'MeTWASModels/',
               tempFolder = tempFolder)
        fff = paste0(tempFolder,list.files(tempFolder))
        fff = fff[!(fff %in%  paste0(tempFolder,'temp',
                                       i,'_SnpObj.bk',
                                c('.bk','.rds')))]
        file.remove(fff)}


        }
system(paste0('rm -r ',tempFolder))

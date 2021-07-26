require(data.table)
require(MOSTWAS)
setwd('/proj/hsantosj/users/abhattac/GWASSummaryStats/ELGAN/coordChange/secondPass')
fff.ss = paste('/proj/hsantosj/users/abhattac/GWASSummaryStats/ELGAN/coordChange/secondPass',
               list.files(),sep='/')

setwd('/proj/hsantosj/users/abhattac/TWAS_res/secondPass')
fff = list.files()

genelocs = fread('/proj/hsantosj/users/abhattac/ELGAN_GWAS/eQTL/ELGANgenelocs.tsv')
snps = fread('/proj/hsantosj/users/abhattac/ELGAN_GWAS/eQTL/ELGAN_snps.tsv')

model.folder = '/proj/hsantosj/users/abhattac/ELGAN_GWAS/eQTL/finModels/'
outFile = '/proj/hsantosj/users/abhattac/secondPassDistal.tsv'
for (i in 1:length(fff)){
    
    print(fff[i])
    twas = fread(fff[i])
    twas = subset(twas, q < .05 & Permute.P < .05)
    sumStats = fread(fff.ss[i])
    sumStats$GenPos = paste(sumStats$CHR_HG38,
                            sumStats$POS_HG38,sep=':')
    if (nrow(twas) > 0){
    for (j in 1:nrow(twas)){
        print(twas$Gene[j])
        qqq = addedLastTest(wgt = paste0(model.folder,
                                         twas$Gene[j],'.wgt.med.RData'),
                            snps = snps,
                            sumStats = sumStats,
                            snpAnnot = NULL,
                            beta = 'BETA',
                            se = 'SE',
                            chr = 'CHR_HG38',
                            pos = 'POS_HG38',
                            ref = 'A2',
                            pval = 'P',
                            R2cutoff = .01,
                            locChrom = genelocs$chr[genelocs$geneid ==
                                                        twas$Gene[j]])
        df = data.frame(Trait = fff[i],
                        Gene = twas$Gene[j],
                        Distal.Z = ifelse(length(qqq)==1,
                                          'NA',qqq$Z.Dist),
                        Distal.P = ifelse(length(qqq)==1,
                                          'NA',qqq$P.Dist))
        fwrite(df,
               outFile,
               sep='\t',
               quote=F,
               append=T,
               row.names=F)
    }}
    
}
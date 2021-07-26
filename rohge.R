require(data.table)
setwd('C:/Users/Arjun/OneDrive - University of North Carolina at Chapel Hill/DOHaD TWAS/Results and Data/twasresults')
fff= paste('C:/Users/Arjun/OneDrive - University of North Carolina at Chapel Hill/DOHaD TWAS/Results and Data/twasresults',
           list.files()[-23],sep='/')
sample.size = c(10555,
                21309,
                298142,
                84689,
                10678,
                9228,
                13960,
                10799,
                9916,
                16469,
                806834,
                694649,
                269867,
                55374,
                52848,
                46350,
                51710,
                72517,
                480359,
                9725,
                150064,
                14307,
                13511,
                19000,
                24013,
                420473,
                420473,
                412960,
                420473,
                135088,
                420473,
                419799)

sample.size = append(sample.size,
                     rep(400000,
                         8))

setwd('C:/Users/Arjun/OneDrive - University of North Carolina at Chapel Hill/DOHaD TWAS/FinalModels/')
qqq = list.files()
FILE = ID = HSQ = c()

for (i in 1:length(qqq)){
    print(qqq[i])
    FILE = c(FILE,qqq[i])
    ID = c(ID,strsplit(qqq[i],'.wgt.')[[1]][1])
    load(qqq[i])
    HSQ = c(HSQ,h2)
}
df.shell = data.frame(FILE = FILE,
                      ID = ID,
                      HSQ = HSQ)

setwd('C:/Users/Arjun/OneDrive - University of North Carolina at Chapel Hill/DOHaD TWAS')
gl = fread('ELGANgenelocs.tsv')
colnames(gl) = c('ID','CHR','P0','P1')

library(liftOver)
library(rtracklayer)
path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
ch = import.chain(path)
df <- data.frame(chr=paste0('chr',gl$CHR), 
                 start=gl$P0, 
                 end=gl$P1, 
                 score=1:nrow(gl),
                 id = gl$ID)

gr  = makeGRangesFromDataFrame(df, 
                               ignore.strand=TRUE,
                               keep.extra.columns=TRUE)
cur19 = as.data.frame(unlist(liftOver(gr, ch)))
cur19 = cur19[order(cur19$id,cur19$width,decreasing = T),]
cur19 = cur19[!duplicated(cur19$id),]
gl = data.frame(ID = cur19$id,
                CHR = sapply(strsplit(as.character(cur19$seqnames),'chr'),
                             function(x) x[2]),
                P0 = cur19$start,
                P1 = cur19$end)
gl = merge(gl,df.shell,by='ID')


require(RHOGE)
setwd('C:/Users/Arjun/OneDrive - University of North Carolina at Chapel Hill/DOHaD TWAS/Results and Data/twasresults/')
ge.cor.df = as.data.frame(matrix(nrow=0,
                                 ncol = 7))
colnames(ge.cor.df) = c('RHOGE','SE','TSTAT','DF','P','Trait1','Trait2')
for (i in 1:length(fff)){
    for (j in (i+1):length(fff)){
        print(c(i,j))
        t1 = fread(fff[i])
        t2 = fread(fff[j])
        n1 = sample.size[i]
        n2 = sample.size[j]
        colnames(t1)[1] = colnames(t2)[1] = 'ID'
        colnames(t1)[2] = colnames(t2)[2] = 'TWAS.Z'
        colnames(t1)[3] = colnames(t2)[3] = 'TWAS.P'
        t1 = merge(t1,gl,by='ID')
        t2 = merge(t2,gl,by='ID')
        ge_cor_res <- tryCatch(rhoge.gw(t1, t2, n1, n2),
                               error=function(e){
                                   data.frame(RHOGE = NA,
                                              SE = NA,
                                              TSTAT = NA,
                                              DF = NA,
                                              P = NA)
                               })
        
        
        tr.sp1 = strsplit(fff[i],'/')[[1]]
        tr1 = strsplit(tr.sp1[length(tr.sp1)],'_TWAS')[[1]][1]
        
        tr.sp2 = strsplit(fff[j],'/')[[1]]
        tr2 = strsplit(tr.sp2[length(tr.sp2)],'_TWAS')[[1]][1]
        
        ge_cor_res$Trait1 = tr1
        ge_cor_res$Trait2 = tr2
        ge.cor.df = rbind(ge.cor.df,ge_cor_res)
    }
}
fwrite(ge.cor.df,
       'C:/Users/Arjun/OneDrive - University of North Carolina at Chapel Hill/DOHaD TWAS/gecor_fromexp.tsv',
       col.names=T,
       row.names=F,
       sep='\t',
       quote=F)

geh2 = geh2.se = geh2.sig = geh2.sig.se = trait = c()
geth2 = function(z,n,m,l){
    return((mean(z^2)-1) * (m/(l*n)))
}
require(bootstrap)
for (i in 1:length(fff)){
    print(i)
    a = fread(fff[i])
    l = 1.4
    n = sample.size[i]
    m = nrow(a)
    hhh = geth2(a$Z,n,m,l)
    geh2 = c(geh2,hhh)
    hhh.se = sd(jackknife(a$Z,geth2,
                          m=m,l=l,n=n)$jack.values)
    geh2.se = c(geh2.se,hhh.se)
    
    a = subset(a,q<.05)
    m = nrow(a)
    hhh.sig = geth2(a$Z,n,m,l)
    geh2.sig = c(geh2.sig,hhh.sig)
    hhh.sig.se = sd(jackknife(a$Z,geth2,
                              m=m,l=l,n=n)$jack.values)
    geh2.sig.se = c(geh2.sig.se,hhh.sig.se)
    
    tr.sp = strsplit(fff[i],'/')[[1]]
    tr = strsplit(tr.sp[length(tr.sp)],'_TWAS')[[1]][1]
    trait = c(trait,tr)
}
fwrite(data.frame(Trait = trait,
                  h2ge = geh2,
                  se = geh2.se,
                  h2ge.sig = geh2.sig,
                  se.sig = geh2.sig.se),
       'C:/Users/Arjun/OneDrive - University of North Carolina at Chapel Hill/DOHaD TWAS/h2_fromexp.tsv',
       col.names=T,
       row.names=F,
       sep='\t',
       quote=F)

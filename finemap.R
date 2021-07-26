#### FOCUS IMPLEMENTATION
require(vroom)
require(data.table)
require(GenomicRanges)
require(Matrix)
snps = fread('/proj/hsantosj/users/abhattac/ELGAN_GWAS/eQTL/ELGAN_snps.tsv')
wgt.folder = '/proj/hsantosj/users/abhattac/ELGAN_GWAS/eQTL/FinalModels_ELGAN'
require(xlsx)
tot.res = read.xlsx('/proj/hsantosj/users/abhattac/SupplementalTables.xlsx',
                    sheetIndex=5)[1:253,1:18]

traits = c("BMI",
           'Body fat percentage',
           "Cholesterol",
           "HDL",
           "Triglycerides",
           'Waist-hip ratio, BMI-adjusted',
           'Diastolic blood pressure',
           "Hypertension")

estimate_cor = function(wmat,ldmat,intercept=F){
    wcov = t(wmat) %*% ldmat %*% wmat
    scale = diag(1/sqrt(diag(wcov)))
    wcor = scale %*% wcov %*% scale
    if (intercept){
        inter = scale %*% t(wmat) %*% ldmat
        return(list(wcor,
                    inter))
    } else {
        return(list(wcor,NA))
    }
}

bayes_factor = function(zscores,
                        idx_set,
                        wcor,
                        prior_chisq = 40,
                        prb = 1e-3,
                        use_log = T){
    m = length(zscores)
    nc = length(idx_set)
    cur_chi2 = prior_chisq/nc
    cur_wcor = wcor[idx_set,idx_set]
    cur_zscors = zscores[idx_set]
    if (nc > 1){
        sss = svd(cur_wcor)
        cur_U = sss$u
        cur_EIG = sss$d
        rm(sss)
    } else {
        cur_U = 1
        cur_EIG = 1
    }
    scaled_chisq = cur_zscors^2

    cur_bf = .5 * -1*sum(log(1 + cur_chi2 %*% cur_EIG)) +
        .5 * sum((cur_chi2 / (1 + cur_chi2 %*% cur_EIG)) * scaled_chisq) +
        nc * log(prb) + (m-nc) * log(1-prb)
    if (use_log){
        return(cur_bf)
    } else {
        return(exp(cur_bf))
    }

}


get_resid = function(zscores, swld, wcor){
    m = nrow(wcor)
    p = ncol(swld)
    intercept = swld %*% rep(1,p)
    wcor_inv = MASS::ginv(as.matrix(wcor))
    rank = Matrix::rankMatrix(wcor_inv)[1]
    numer = t(intercept) %*% wcor_inv %*% zscores
    denom = t(intercept) %*% wcor_inv %*% intercept
    alpha = numer/denom
    resid = zscores -
        (t(intercept) * alpha)
    s2 = (resid %*% wcor_inv %*%
              t(resid))/(rank-1)
    inter_se = sqrt(s2/denom)
    inter_z = alpha/inter_se
    return(list(resid,
                inter_z))
}

file.remove('credible_set.tsv')
for (t in 1:length(traits)){
    tr = traits[t]
    res = subset(tot.res,Trait == tr)
    res = res[order(res$Chromosome,res$Start),]
    chr.table = table(res$Chromosome)
    res = subset(res,Chromosome %in%
                     as.numeric(names(which(table(res$Chromosome)>1))))
    chr.un = unique(res$Chromosome)
    keep.genes = c()
    for (c in chr.un){

        res.cur = subset(res,Chromosome == c)
        res.cur = res.cur[order(res.cur$Start),]
        for (i in 1:(nrow(res.cur)-1)){
            if (res.cur$End[i] > res.cur$Start[i+1] - 1e6){
                keep.genes = unique(c(keep.genes,
                                      c(res.cur$Gene[c(i,i+1)])))
            }
        }
    }

    res = subset(res,Gene %in% keep.genes)
    all.snps = c()
    omega = c()
    gene = c()
    snp.chr = c()
    chr.un = unique(res$Chromosome)
    if (length(chr.un) >= 1){
    for (c in chr.un){
        res = subset(res,Chromosome == c)
        for (i in 1:nrow(res)){
            load(paste(wgt.folder,
                       paste0(res$Gene[i],'.wgt.med.RData'),
                       sep = '/'))
            all.snps = c(all.snps,
                         as.character(Model$SNP))
            omega = c(omega,
                      as.numeric(Model$Effect))
            gene = c(gene,
                     rep(res$Gene[i],nrow(Model)))
            snp.chr = c(snp.chr,
                        as.numeric(Model$Chromosome))

        }
        tot.df = data.frame(SNP = all.snps,
                            Gene = gene,
                            Effect = omega,
                            Chromosome = snp.chr)
        model.df = as.data.frame(matrix(nrow = length(unique(all.snps)),
                                        ncol = nrow(res)+1))
        colnames(model.df) = c('SNP',res$Gene)
        model.df$SNP = as.character(unique(all.snps))
        for (g in 1:nrow(res)){
            print(res$Gene[g])
            cur.tot.df = subset(tot.df,Gene == res$Gene[g])
            cur.tot.df$SNP = as.character(cur.tot.df$SNP)
            for (i in 1:nrow(model.df)){
                print(i)
                w = which(cur.tot.df$SNP == model.df$SNP[i])
                model.df[i,g+1] = ifelse(length(w) != 0,
                                         cur.tot.df$Effect[w],
                                         0)
            }
        }

        model.df$Chromosome = c
        for (i in 1:nrow(model.df)){
            rrr = subset(tot.df,SNP == model.df$SNP[i])
            model.df$Chromosome[i] = rrr$Chromosome[1]
        }
        snp.set = as.data.frame(subset(snps,
                                       SNP %in% model.df$SNP))
        snp.set = snp.set[match(model.df$SNP,snp.set$SNP),]
        snpMat = as.matrix(snp.set[,-1])
        rm(snp.set)
        snpMat = Matrix(t(scale(t(snpMat))), sparse = T)
        V = tcrossprod(snpMat) / ncol(snpMat)
        if (any(model.df$Chromsome == c)){
        V[model.df$Chromosome == c,model.df$Chromosome != c] = 0
        V[model.df$Chromosome != c,model.df$Chromosome == c] = 0
        } else {
            V = diag(diag(V),nrow=nrow(V))
        }
        Omega = Matrix(as.matrix(model.df[,-c(1,ncol(model.df))]))
        zscores = res$Z
        m = length(zscores)
        wcor = estimate_cor(Omega,V,intercept=T)[[1]]
        swld = estimate_cor(Omega,V,intercept=T)[[2]]
        null_res = m * log(1 - 1e-3)
        marginal = m * log(1 - 1e-3)
        comb_list = list()
        for (n in 1:min(3,length(zscores))){
            comb_list = c(comb_list,
                          combn(1:length(zscores),n,simplify=F))
        }

        pips = rep(0,length(zscores))
        zscores = get_resid(zscores,swld,wcor)[[1]]
        for (j in 1:length(comb_list)){
            subset = comb_list[[j]]
            local = bayes_factor(zscores,
                                 idx_set = subset,
                                 wcor = wcor)
            marginal = log(exp(local) + exp(marginal))
            for (idx in subset){
                if (pips[idx] == 0){
                    pips[idx] = local
                } else {
                    pips[idx] = log(exp(pips[idx]) + exp(local))
                }
            }
            print(pips)
            print(marginal)
        }
        
        pips = exp(pips - marginal)
        null_res = exp(null_res - marginal)
        res$pip = pips
        res = res[order(res$pip,decreasing = T),]
        npost = res$pip/sum(res$pip)
        csum = cumsum(npost)
        res$in_cred_set = csum < .9
        fwrite(res,'credible_set.tsv',append=T,
               row.names= F,quote=F,sep='\t')

        res = subset(tot.res,Trait == tr)
        res = subset(res,Gene %in% keep.genes)
    }
    }
}

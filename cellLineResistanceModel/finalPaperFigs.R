##process cytokine data
library(amlresistancenetworks)
library(dplyr)


##first run gilteritinib data
protData<-querySynapseTable('syn22156807')%>%mutate(Gene=unlist(Gene))%>%
  dplyr::rename(sample='Sample')%>%
  dplyr::rename(LogRatio='value')
#readRDS(system.file('gilteritinibData.Rds',package='amlresistancenetworks'))
#prot.univ<-unique(gilt.data$Gene)
phosData<-querySynapseTable('syn22156809')%>%subset(!is.nan(value))%>%
  mutate(Gene=unlist(Gene))%>%
  mutate(site=unlist(site))%>%
  dplyr::rename(LogRatio='value')%>%
  dplyr::rename(sample='Sample')

syn<-synapseLogin()
targ.data<-read.csv2(syn$get('syn22214245')$path,sep='\t',header=T)%>%
  dplyr::rename(Gene='T..Category')%>%
  tidyr::pivot_longer(-Gene,names_to='sample',values_to='value')%>%
  mutate(value=as.numeric(value))%>%
  rowwise()%>%
  mutate(samp2=stringr::str_replace(sample,'_3June20.+',''))%>%
  mutate(Patient=stringr::str_replace(samp2,'.*Ex16_',''))%>%
  dplyr::select(-c(samp2,sample))

norm.data<-targ.data%>%tidyr::pivot_wider(values_from=value,names_from=Patient)%>%
  tidyr::pivot_longer(-c(Gene,pool),values_to='value',names_to='Patient')%>%
  rowwise()%>%
  mutate(logRatio=log10(value)-log10(pool))%>%
  dplyr::select(-c(pool,value))%>%
  dplyr::rename(value='logRatio')


protMat<-protData%>%dplyr::select(sample,Gene,LogRatio)%>%
  tidyr::pivot_wider(values_from=LogRatio,names_from=sample,
                     values_fn=list(LogRatio=mean),values_fill=list(LogRatio=0.0))%>%
  tibble::column_to_rownames('Gene')

phosMat<-phosData%>%dplyr::select(sample,site,LogRatio)%>%
  tidyr::pivot_wider(values_from=LogRatio,names_from=sample,
                     values_fn=list(LogRatio=mean),values_fill=list(LogRatio=0.0))%>%
  tibble::column_to_rownames('site')

kindat<-mapPhosphoToKinase(dplyr::rename(phosData,LogFoldChange='LogRatio')%>%rename(Sample='sample'))

##
#' @param dat.table
plotAllData<-function(dat.table,vars=c('sample','cellLine','ligand','treatment')){
  library(ggfortify)
  met<-dat.table%>%dplyr::select(vars)%>%
    distinct()
  #%>%
  #  tibble::column_to_rownames('sample')
    
  mat<-dat.table%>%dplyr::select(Gene,LogRatio,sample)%>%
    distinct()%>%
    mutate(LogRatio=as.numeric(LogRatio))%>%
    tidyr::pivot_wider(names_from='sample',values_from='LogRatio',
                       values_fn=list(LogRatio=function(x) mean(x,na.rm=T)),values_fill=list(LogRatio=0))%>%
    tibble::column_to_rownames('Gene')
  
  autoplot(prcomp(t(mat)),data=met,colour='treatment',shape='cellLine')
 
}

#' plotKSEAinHeatmap
#' Plots results of KSEA in heatmap across all kinase data
#' @param kindat
#' @param ksres
#' @param vars to plot
plotKSEAinHeatmap<-function(kindat,ksres,vars){
  
  ksres%>%subset(p.value<0.05)%>%select(Kinase.Gene,Condition)%>%
    distinct()%>%
    group_by(Condition)%>%
    group_map(~ plotKinDat(kindat,prefix=.y$Condition,vars=vars,genelist=.x$Kinase.Gene))
  
}

#' plots all kinase activity in a heatmap
#' @param kindat list of kinases and their summarized activity in each sample
#' @param prefix
#' @param vars clinical variables
#' @param genelist filter for these genes
plotKinDat<-function(kindat,prefix='all',vars=c('sample','time','ligand','treatment'),genelist=c()){
  library(pheatmap)
  ##create matrix of kinase scores
  if(length(genelist)>0)
    kindat<-kindat%>%subset(Kinase%in%genelist)
  print(genelist)
  
  mat <-kindat%>%
    ungroup()%>%
    tidyr::pivot_wider(-c(meanNKINscore,numSubstr),
                        values_from=meanLFC,
                        names_from=Sample,
                      values_fn=list(meanLFC=mean))%>%
    tibble::column_to_rownames('Kinase')
  
  kinAts<-kindat%>%ungroup()%>%dplyr::select(Kinase,numSubstr)%>%distinct()%>%
    group_by(Kinase)%>%summarize(substrates=mean(numSubstr))%>%
        tibble::remove_rownames()%>%
  tibble::column_to_rownames('Kinase')
  
  sampAts<-phosData%>%dplyr::select(vars)%>%
    distinct()%>%
        tibble::remove_rownames()%>%
    tibble::column_to_rownames('sample')
  
  #sampAts$TimePoint=as.factor(sampAts$TimePoint)
  all.vars<-sort(apply(mat,1,var),decreasing=T)
  all.vars<-all.vars[which(all.vars!=0)]
  if(length(genelist)==0){
    vars=names(all.vars)[1:150]
  }else{
    vars=intersect(genelist,names(all.vars))
  }
  
 pheatmap(mat[vars,],cellwidth = 8,cellheight=8,clustering_distance_cols = 'correlation',
          clustering_distance_rows = 'correlation',
          annotation_row = kinAts,annotation_col=sampAts,
          file=paste0(prefix,'KinaseHeatmap.pdf'),height=20,width=8) 
}



#' plot all the KSEA 
#' @param condList
#' @return data frame
doAllKSEAplots<-function(condList,pdat=phosData){
  
  gene.to.site<-dplyr::select(pdat,Gene,site,Peptide)%>%distinct()%>%
    dplyr::mutate(residue=stringr::str_replace(site,paste0(Gene,'-'),''))%>%
    dplyr::mutate(residue=stringr::str_replace_all(residue,"([STY])", ";\\1"))%>%
    dplyr::mutate(residue=stringr::str_replace(residue,"^;", ""))%>%
    dplyr::mutate(residue=stringr::str_replace_all(residue,"([sty])", ""))
  
  full.df<-purrr::map_df(names(condList),.f=function(clName){ 
    condList[[clName]]%>%
      tibble::rownames_to_column('site')%>%
      left_join(gene.to.site)%>%
      dplyr::select(Gene,Peptide,residue,value='logFC',p_adj='adj.P.Val')%>%
      amlresistancenetworks::computeKSEA(.,prefix=clName,0.05)%>%
      mutate(Condition=clName)%>%
      as.data.frame()
  })
  return(full.df)
  
}

  
  plots=list(plotAllData(protData),plotAllData(phosData))
  cowplot::plot_grid(plotlist=plots,labels=c("Bulk Proteomics",'Phosphoprotomics'),nrow=2)
  ggsave('pcaOfSamples.png')
  

      clinvars = c("sample","cellLine","treatment","ligand")
    plotKinDat(kindat,prefix='giltResistance',
                                    vars=clinvars)


  ##what are we doing again?
  summary<-protData%>%dplyr::select(clinvars)%>%distinct()%>%rowwise()%>%
    mutate(Condition=stringr::str_c(cellLine,treatment,ligand,sep='_'))
  print(summary)
  
  
   earlyLateProt<-list(gilt_late_early_flt3 =limmaTwoFactorDEAnalysis(protMat,
                                                        filter(summary,Condition=='MOLM14_Early Gilteritinib_FLT3')$sample,
                                                        filter(summary,Condition=='MOLM14_Late Gilteritinib_FLT3')$sample),
                      gilt_late_early_fgf=limmaTwoFactorDEAnalysis(protMat,
                                                      filter(summary,Condition=='MOLM14_Early Gilteritinib_FGF2')$sample,
                                                      filter(summary,Condition=='MOLM14_Late Gilteritinib_FGF2')$sample),
                      gilt_late_early_combined=limmaTwoFactorDEAnalysis(protMat,
                                                           filter(summary,Condition%in%c('MOLM14_Early Gilteritinib_FGF2','MOLM14_Early Gilteritinib_FLT3',
                                                                                         'MV411_Early Gilteritinib_FLT3','MV411_Early Gilteritinib_FGF2'))$sample,
                                                           filter(summary,Condition%in%c('MOLM14_Late Gilteritinib_FGF2','MOLM14_Late Gilteritinib_FLT3',
                                                                                         'MV411_Late Gilteritinib_FGF2','MV411_Late Gilteritinib_FLT3'))$sample))
  
  
earlyLateMolmPhos<-list(gilt_early_flt3_molm14 =limmaTwoFactorDEAnalysis(phosMat,
                                                    filter(summary,Condition=='MOLM14_None_None')$sample,
                                                    filter(summary,Condition=='MOLM14_Early Gilteritinib_FLT3')$sample),
                 gilt_early_fgf_molm14=limmaTwoFactorDEAnalysis(phosMat,
                                                  filter(summary,Condition=='MOLM14_None_None')$sample,
                                                  filter(summary,Condition=='MOLM14_Early Gilteritinib_FGF2')$sample),
                 gilt_late_flt3_molm14=limmaTwoFactorDEAnalysis(phosMat,
                                                  filter(summary,Condition%in%c('MOLM14_None_None'))$sample,
                                                  filter(summary,Condition%in%c('MOLM14_Late Gilteritinib_FLT3'))$sample),
                  gilt_late_fgf_molm14=limmaTwoFactorDEAnalysis(phosMat,
                                                    filter(summary,Condition%in%c('MOLM14_None_None'))$sample,
                                                    filter(summary,Condition%in%c('MOLM14_Late Gilteritinib_FGF2'))$sample),
                 gilt_late_vs_early_fgf_moml14=limmaTwoFactorDEAnalysis(phosMat,
                                                                        filter(summary,Condition%in%c('MOLM14_Early Gilteritinib_FGF2'))$sample,
                                                                        filter(summary,Condition%in%c('MOLM14_Late Gilteritinib_FGF2'))$sample),
                  gilt_late_vs_early_flt3_molm14=limmaTwoFactorDEAnalysis(phosMat,
                                                                    filter(summary,Condition%in%c('MOLM14_Early Gilteritinib_FLT3'))$sample,
                                                                          filter(summary,Condition%in%c('MOLM14_Late Gilteritinib_FLT3'))$sample))


earlyLateMv411Phos<-list( gilt_early_flt3_mv411 =limmaTwoFactorDEAnalysis(phosMat,
                                                             filter(summary,Condition=='MV411_None_None')$sample,
                                                             filter(summary,Condition=='MV411_Early Gilteritinib_FLT3')$sample),
                         gilt_early_fgf_mv411=limmaTwoFactorDEAnalysis(phosMat,
                                                           filter(summary,Condition=='MV411_None_None')$sample,
                                                           filter(summary,Condition=='MV411_Early Gilteritinib_FGF2')$sample),
                         gilt_late_flt3_mv411=limmaTwoFactorDEAnalysis(phosMat,
                                                           filter(summary,Condition%in%c('MV411_None_None'))$sample,
                                                           filter(summary,Condition%in%c('MV411_Late Gilteritinib_FLT3'))$sample),
                         gilt_late_fgf_mv411=limmaTwoFactorDEAnalysis(phosMat,
                                                          filter(summary,Condition%in%c('MV411_None_None'))$sample,
                                                          filter(summary,Condition%in%c('MV411_Late Gilteritinib_FGF2'))$sample),
                         gilt_late_vs_early_fgf_mv411=limmaTwoFactorDEAnalysis(phosMat,
                                                                        filter(summary,Condition%in%c('MV411_Early Gilteritinib_FGF2'))$sample,
                                                                        filter(summary,Condition%in%c('MV411_Late Gilteritinib_FGF2'))$sample),
                         gilt_late_vs_early_flt3_mv411=limmaTwoFactorDEAnalysis(phosMat,
                                                                          filter(summary,Condition%in%c('MV411_Early Gilteritinib_FLT3'))$sample,
                                                                          filter(summary,Condition%in%c('MV411_Late Gilteritinib_FLT3'))$sample))

earlyLateCombPhos<-list( gilt_early_vs_parental_comb =limmaTwoFactorDEAnalysis(phosMat,
                                                              filter(summary,Condition%in%c("MV411_None_None","MOLM14_None_None"))$sample,
                                                              filter(summary,Condition%in%c('MOLM14_Early Gilteritinib_FLT3','MV411_Early Gilteritinib_FLT3',
                                                                                            'MOLM14_Early Gilteritinib_FGF2','MV411_Early Gilteritinib_FGF2'))$sample),
                     gilt_late_vs_parental_comb=limmaTwoFactorDEAnalysis(phosMat,
                                                            filter(summary,Condition%in%c("MV411_None_None","MOLM14_None_None"))$sample,
                                                            filter(summary,Condition%in%c('MOLM14_Late Gilteritinib_FLT3','MV411_Late Gilteritinib_FLT3',
                                                                                          'MOLM14_Late Gilteritinib_FGF2','MV411_Late Gilteritinib_FGF2'))$sample),
                     gilt_late_vs_early_comb=limmaTwoFactorDEAnalysis(phosMat,
                                               filter(summary,Condition%in%c('MOLM14_Early Gilteritinib_FLT3','MV411_Early Gilteritinib_FLT3',
                                                                             'MOLM14_Early Gilteritinib_FGF2','MV411_Early Gilteritinib_FGF2'))$sample,
                                               filter(summary,Condition%in%c('MOLM14_Late Gilteritinib_FLT3','MV411_Late Gilteritinib_FLT3',
                                                                             'MOLM14_Late Gilteritinib_FGF2','MV411_Late Gilteritinib_FGF2'))$sample))


  
ph2<-doAllKSEAplots(earlyLateCombPhos)
plotKSEAinHeatmap(kindat,ph2,clinvars)

ph3<-doAllKSEAplots(earlyLateMv411Phos)
plotKSEAinHeatmap(kindat,ph3,clinvars)

ph4<-doAllKSEAplots(earlyLateMolmPhos)
plotKSEAinHeatmap(kindat,ph4,clinvars)
##plot single kinase/substrate expression of mapk3, mapk1, and mapk8
p5<-kindat%>%
  left_join(dplyr::rename(summary,Sample='sample'))%>%
  subset(Kinase%in%c('CDC7','AURKB'))%>%
  ggplot(aes(x=as.factor(ligand),y=meanLFC,fill=treatment))+
  geom_boxplot()+
  facet_grid(~Kinase)+scale_fill_viridis_d()+facet_grid(~cellLine)+
  ggtitle("Estimated Kinase Activity")
ggsave('estimatedGiltAurbActivity.png',p5,width=10)

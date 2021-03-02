##process cytokine data
library(amlresistancenetworks)
library(dplyr)
library(nationalparkcolors)
pal<-park_palette('Saguaro',5)

##first run gilteritinib data
protData<-querySynapseTable('syn24189419')%>%
  mutate(Gene=unlist(Gene))%>%
    dplyr::rename(time='Time (minutes)')%>%
  dplyr::rename(sample='Sample')
  #dplyr::rename(LogRatio='value')
#readRDS(system.file('gilteritinibData.Rds',package='amlresistancenetworks'))
#prot.univ<-unique(gilt.data$Gene)
phosData<-querySynapseTable('syn24189487')%>%#subset(!is.nan(LogRatio))%>%
  mutate(Gene=unlist(Gene))%>%
  mutate(site=unlist(site))%>%
  dplyr::rename(time='Time (minutes)')%>%
 # dplyr::rename(LogRatio='value')%>%
  dplyr::rename(sample='Sample')


kindat<-mapPhosphoToKinase(dplyr::rename(phosData,LogFoldChange='LogRatio',
                                         Sample='sample'))

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
                       values_fn=list(LogRatio=function(x) mean(x,na.rm=T)),
                       values_fill=list(LogRatio=0))%>%
    tibble::remove_rownames()%>%
    tibble::column_to_rownames('Gene')
  
  autoplot(prcomp(t(mat)),data=met,shape='Treatment',colour='time')
 
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



protMat<-protData%>%dplyr::select(sample,Gene,LogRatio)%>%
  tidyr::pivot_wider(values_from=LogRatio,names_from=sample,
                     values_fn=list(LogRatio=mean),values_fill=list(LogRatio=0.0))%>%
  tibble::column_to_rownames('Gene')

phosMat<-phosData%>%dplyr::select(sample,site,LogRatio)%>%
  tidyr::pivot_wider(values_from=LogRatio,names_from=sample,
                     values_fn=list(LogRatio=mean),values_fill=list(LogRatio=0.0))%>%
  tibble::column_to_rownames('site')


library(ggplot2)
library(ggalluvial)
library(ggridges)
#' plotConditionsInFlow
#' Tries to compare multiple logfc values
plotConditionsInFlow<-function(condList,title='',pvalThresh=0.05,upDown='logFC'){
  
  full.df<-purrr::map_df(names(condList),.f=function(clName){ 
    condList[[clName]]%>%
      tibble::rownames_to_column('Gene')%>%
      dplyr::select(Gene,logFC,adj.P.Val)%>%
      mutate(Condition=clName)%>%
      mutate(direction=ifelse(logFC>0,'UpReg','DownReg'))%>%
      mutate(signif=ifelse(adj.P.Val<pvalThresh,'Significant','Non-significant'))
  })
  
  p1<-ggplot(full.df,aes(x=logFC,y=Condition))+geom_density_ridges_gradient(aes(fill=signif))
  
  res.df<-full.df%>%
    subset(adj.P.Val<pvalThresh)
  if(nrow(res.df)==0)
    return(p1)
  
  p2<-res.df%>%#subset(Gene%in%changing$Gene)%>%
    ggplot(aes(x=Condition,stratum=direction,alluvium=Gene,fill=direction,label=Condition,alpha=0.5))+
    geom_flow(stat='alluvium',lode.guidance='frontback')+
    geom_stratum()+
    theme_minimal()+
    viridis::scale_fill_viridis(3,discrete=T)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  p<-cowplot::plot_grid(p1,p2,nrow=1)+ggtitle(title)
  return(p)
}



#' plot all the GO 
#' @param condList
#' @return data frame
doAllGOplots<-function(condList){
  
  full.df<-purrr::map_df(names(condList),.f=function(clName){ 
    condList[[clName]]%>%
      tibble::rownames_to_column('Gene')%>%
      dplyr::select(Gene,value='logFC')%>%
      amlresistancenetworks::plotOldGSEA(.,prefix=clName,0.05)%>%
      as.data.frame()
  })
  return(full.df)
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

#' build networks from data frame
#' @param data.res
#' @param gene.col
#' @param weight.col
#' @param condition.col
#' @return network list?
runNetworksFromDF<-function(data,gene.col='Kinase.Gene',
                              weight.col='aveSubstrateLog2FC',
                              condition.col='Condition',
                            extra.col=c('Substrate.Gene','Source','log2FC'),
                              signif=0.05){
  res = data%>%
   # dplyr::select(cond=condition.col,value=weight.col,Gene=gene.col,p.value)%>%
    mutate(signif=p.value<signif)%>%
      dplyr::select(c(condition.col,weight.col,gene.col,'signif',extra.col))%>%distinct()%>%
    dplyr::rename(cond=condition.col,value=weight.col,Gene=gene.col)%>%
    dplyr::select(c('cond','Gene','value',extra.col,'signif'))%>%
    group_by(cond)%>%
    group_map(~ amlresistancenetworks::computeProteinNetwork(.x),keep=TRUE)
  return(res)
}


  clinvars = c("sample","time","Treatment","Ligand")
  ##what are we doing again?
  summary<-protData%>%dplyr::select(clinvars)%>%distinct()%>%rowwise()%>%
    mutate(Condition=stringr::str_c(Treatment,Ligand,time,sep='_'))
  print(summary)
  doPlots=FALSE
if(doPlots){

  plotKinDat(kindat,prefix='giltTimeCourse',vars=c('sample','Treatment','Ligand','time'))
  
  plots=list(plotAllData(protData,vars=c('sample','Treatment','Ligand','time')),plotAllData(phosData,vars=c('sample','Treatment','Ligand','time')))
  cowplot::plot_grid(plotlist=plots,labels=c("Bulk Proteomics",'Phosphoprotomics'),
                     nrow=2)
  ggsave('pcaOfGiltTimeCourse.png')
  

  
   earlyTimeProt<-list(thirty_vs_zero_none=limmaTwoFactorDEAnalysis(protMat,
                                                        filter(summary,Condition=='none_None_0')$sample,
                                                        filter(summary,Condition=='Gilteritinib_None_30')$sample),
                      thirs_vs_zero_fgf2=limmaTwoFactorDEAnalysis(protMat,
                                                      filter(summary,Condition=='none_FGF2_0')$sample,
                                                      filter(summary,Condition=='Gilteritinib_FGF2_30')$sample),
                      thirty_vs_zero_flt3=limmaTwoFactorDEAnalysis(protMat,
                                                           filter(summary,Condition%in%c('none_FLT3_0'))$sample,
                                                           filter(summary,Condition%in%c('Gilteritinib_FLT3_30'))$sample),
                      thirty_vs_zero_combined=limmaTwoFactorDEAnalysis(protMat,
                                                                       filter(summary,Condition%in%c('none_FLT3_0','none_FGF2_0'))$sample,
                                                                       filter(summary,Condition%in%c('Gilteritinib_FLT3_30','Gilteritinib_FGF2_30'))$sample))
   earlyTimePhos<-list(thirty_vs_zero_none=limmaTwoFactorDEAnalysis(phosMat,
                                                                    filter(summary,Condition=='none_None_0')$sample,
                                                                    filter(summary,Condition=='Gilteritinib_None_30')$sample),
                       thirs_vs_zero_fgf2=limmaTwoFactorDEAnalysis(phosMat,
                                                                   filter(summary,Condition=='none_FGF2_0')$sample,
                                                                   filter(summary,Condition=='Gilteritinib_FGF2_30')$sample),
                       thirty_vs_zero_flt3=limmaTwoFactorDEAnalysis(phosMat,
                                                                    filter(summary,Condition%in%c('none_FLT3_0'))$sample,
                                                                    filter(summary,Condition%in%c('Gilteritinib_FLT3_30'))$sample),
                       thirty_vs_zero_combined=limmaTwoFactorDEAnalysis(phosMat,
                                                                       filter(summary,Condition%in%c('none_FLT3_0','none_FGF2_0'))$sample,
                                                                       filter(summary,Condition%in%c('Gilteritinib_FLT3_30','Gilteritinib_FGF2_30'))$sample))

   lateTimePhos<-list(threeHr_vs_zero_none=limmaTwoFactorDEAnalysis(phosMat,
                                                                   filter(summary,Condition=='none_None_0')$sample,
                                                                   filter(summary,Condition=='Gilteritinib_None_180')$sample),
                      threeHr_vs_zero_fgf2=limmaTwoFactorDEAnalysis(phosMat,
                                                                  filter(summary,Condition=='none_FGF2_0')$sample,
                                                                  filter(summary,Condition=='Gilteritinib_FGF2_180')$sample),
                      threeHr_vs_zero_flt3=limmaTwoFactorDEAnalysis(phosMat,
                                                                   filter(summary,Condition%in%c('none_FLT3_0'))$sample,
                                                                   filter(summary,Condition%in%c('Gilteritinib_FLT3_180'))$sample),
                      threeHr_vs_zero_combined=limmaTwoFactorDEAnalysis(phosMat,
                                                                       filter(summary,Condition%in%c('none_FLT3_0','none_FGF2_0'))$sample,
                                                                       filter(summary,Condition%in%c('Gilteritinib_FLT3_180','Gilteritinib_FGF2_180'))$sample))
   
   
      earlyLatePhos<-list(threeHr_vs_thirty_none=limmaTwoFactorDEAnalysis(phosMat,
                                                                    filter(summary,Condition=='Gilteritinib_None_30')$sample,
                                                                    filter(summary,Condition=='Gilteritinib_None_180')$sample),
                       threeHr_vs_thirty_fgf2=limmaTwoFactorDEAnalysis(phosMat,
                                                                   filter(summary,Condition=='Gilteritinib_FGF2_30')$sample,
                                                                   filter(summary,Condition=='Gilteritinib_FGF2_180')$sample),
                       threeHr_vs_thirty_flt3=limmaTwoFactorDEAnalysis(phosMat,
                                                                    filter(summary,Condition%in%c('Gilteritinib_FLT3_30'))$sample,
                                                                    filter(summary,Condition%in%c('Gilteritinib_FLT3_180'))$sample),
                                                 threeHr_vs_thirty_combined=limmaTwoFactorDEAnalysis(phosMat,
                                                                       filter(summary,Condition%in%c('Gilteritinib_FLT3_30','Gilteritinib_FGF2_30'))$sample,
                                                                       filter(summary,Condition%in%c('Gilteritinib_FLT3_180','Gilteritinib_FGF2_180'))$sample))

      earlyLateProt<-list(threeHr_vs_thirty_none=limmaTwoFactorDEAnalysis(protMat,
                                                                          filter(summary,Condition=='Gilteritinib_None_30')$sample,
                                                                          filter(summary,Condition=='Gilteritinib_None_180')$sample),
                          threeHr_vs_thirty_fgf2=limmaTwoFactorDEAnalysis(protMat,
                                                                          filter(summary,Condition=='Gilteritinib_FGF2_30')$sample,
                                                                          filter(summary,Condition=='Gilteritinib_FGF2_180')$sample),
                          threeHr_vs_thirty_flt3=limmaTwoFactorDEAnalysis(protMat,
                                                                          filter(summary,Condition%in%c('Gilteritinib_FLT3_30'))$sample,
                                                                          filter(summary,Condition%in%c('Gilteritinib_FLT3_180'))$sample),
                          threeHr_vs_thirty_combined=limmaTwoFactorDEAnalysis(protMat,
                                                                       filter(summary,Condition%in%c('Gilteritinib_FLT3_30','Gilteritinib_FGF2_30'))$sample,
                                                                       filter(summary,Condition%in%c('Gilteritinib_FLT3_180','Gilteritinib_FGF2_180'))$sample))

      

p1<-plotConditionsInFlow(earlyTimeProt,title='30 min results',0.05)
ggsave("thirtyMinProt.png",width=11,height=6)
#gp<-doAllGOplots(earlyTimeProt)


p3<-plotConditionsInFlow(earlyTimePhos,title='30 min results',0.05)
ggsave('thirtyMinPhos.png',p3,width=11,height=6)
ph3<-doAllKSEAplots(earlyTimePhos)
plotKSEAinHeatmap(kindat,ph3,clinvars)

p2<-plotConditionsInFlow(earlyLateProt,title='180 min results',0.05)
ggsave("lateVsthirtyMinProt.png",width=11,height=6)
#gp<-doAllGOplots(earlyTimeProt)


p4<-plotConditionsInFlow(earlyLatePhos,title='180 min results',0.05)
ggsave('lateVsthirtyMinPhos.png',p4,width=11,height=6)
ph4<-doAllKSEAplots(earlyLatePhos)
plotKSEAinHeatmap(kindat,ph4,clinvars)

p5<-plotConditionsInFlow(lateTimePhos,title='180 min vs 0',0.05)
ggsave('lateVsZeroPhos.png',p5)
ph5<-doAllKSEAplots(lateTimePhos)
plotKSEAinHeatmap(kindat,ph5,clinvars)

#nets<-ph3%>%mutate(Condition=stringr::str_c(Condition,'_phos'))%>%
#xs  runNetworksFromDF()

#resdf<-do.call(rbind,lapply(names(earlyLateProt),function(x) data.frame(earlyLateProt[[x]],Condition=x)))

#pnets<-resdf%>%mutate(Condition=stringr::str_c(Condition,'_prot'))%>%
#  dplyr::rename(p.value='adj.P.Val')%>%
#  runNetworksFromDF(.,gene.col='featureID',weight.col='logFC',condition.col='Condition',extra.col=c('AveExpr','t','B','P.Value'),signif=0.01)

#mcp1Resistnetworks<-runNetworksFromDF(ph3)

##plot single kinase/substrate expression of mapk3, mapk1, and mapk8
p5<-kindat%>%
  left_join(dplyr::rename(summary,Sample='sample'))%>%
  subset(Kinase%in%c('CDC7','AURKB','MTOR'))%>%
  ggplot(aes(x=time,y=meanLFC,color=Treatment,shape=Kinase))+
  geom_jitter()+
  facet_grid(~Ligand)+scale_color_viridis_d()+
  ggtitle("Estimated Kinase Activity")
ggsave('estimatedAurbActivityGiltTime.png',p5,width=10)

t0<-summary%>%
    subset(time==0)%>%mutate(t2='DMSO')%>%
  select(-Treatment)%>%rename(Treatment='t2')
t1<-summary%>%
  subset(time==0)%>%mutate(t2='Gilteritinib')%>%
  select(-Treatment)%>%rename(Treatment='t2')

newsumm<-summary%>%
  subset(time!=0)%>%rbind(t0)%>%rbind(t1)

p6<-kindat%>%
    left_join(dplyr::rename(newsumm,Sample='sample'))%>%
    subset(Kinase%in%c('CDC7','AURKB','MTOR'))%>% 
    group_by(Ligand,time,Kinase,Treatment)%>%
    summarize(mlfc=mean(meanLFC),sd=sd(meanLFC))%>%
    ggplot(aes(x=time,y=mlfc,color=Kinase,shape=Treatment,linetype=Treatment))+
  geom_line(aes(size=.25))+geom_point(aes(size=.4))+scale_color_manual(values=pal)+
  facet_grid(~Ligand)+
   geom_errorbar(aes(ymin=mlfc-sd, ymax=mlfc+sd))
ggsave('timeCourseKinaseActivity.pdf',p6,width=12)

}
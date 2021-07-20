library(amlresistancenetworks)
library(dplyr)
library(kableExtra)
library(tibble)
library(leapr)
library(ggplot2)
library(scales)
library(gridExtra)
library(readr)
library(rstudioapi)
source("../../../Util/synapseUtil.R")

############## This script generates the KSEA heatmaps seen in file 14.

########################### Setting up ##########################

## Saving plots to the folder in which this script is located!!
current_path <- getActiveDocumentContext()$path 
setwd(dirname(current_path))

# corrected for loading mass and plexID, long time and subtle
p.site.corrected <- querySynapseTable("syn25706631")

phos.dat <- p.site.corrected

phos.dat <- p.site.corrected
phos.dat$Peptide <- sub("^.*@(.*)$", "\\1", phos.dat$Accession)
phos.dat$site <- sub("^(.*)@.*$", "\\1", phos.dat$Accession)
phos.dat$Gene <- sub("(.*)-.*$", "\\1", phos.dat$site)

# from here we can select just the metadata of interest
metadata.columns = c('Sample','BeatAML Patient ID','Plex','Loading Mass', 'Description')
summary <- phos.dat%>%
  select(metadata.columns)%>%
  distinct()
summary$Description <- as.character(summary$Description)
summary$Type <- case_when(summary$`BeatAML Patient ID` == "Healthy Donor Stroma" ~ "Healthy Donor Stroma",
                          summary$`BeatAML Patient ID` == "Cell line" ~ "Cell line")
summary[is.na(summary$Type), "Type"] <- "Treated"

## Elie has corrected the time period annotation. Adding a pre-treatment, early, and late period.
correct.time.annotation <- as.data.frame(read_csv("../../Updated time period annotation.csv"))
summary$Period <- correct.time.annotation[summary$Sample, "Classification"]
summary$Period <- case_when(grepl("early", summary$Period) ~ "Early",
                            grepl("Late", summary$Period) ~ "Late",
                            grepl("pre-treatment", summary$Period) ~ "Pre-Treatment")

##select the samples of interest
healthy <- subset(summary, Type == 'Healthy Donor Stroma')%>%
  select(Sample)%>%
  unique()

patients <- subset(summary, Type == "Treated") %>%
  select(Sample, `BeatAML Patient ID`, Period, Plex) %>%
  dplyr::rename(ID = `BeatAML Patient ID`) %>%
  mutate(Plex = as.character(Plex)) %>%
  unique()

##then we spread the proteomics data into a matrix
phos.mat <- phos.dat%>%
  select(LogRatio, Sample, site)%>%
  tidyr::pivot_wider(values_from='LogRatio',names_from='Sample',
                     values_fn=list(LogRatio=mean),values_fill=0.0)%>%
  tibble::column_to_rownames('site')

gene.to.site<-dplyr::select(phos.dat,Gene,site,Peptide)%>%distinct()%>%
  dplyr::mutate(residue=stringr::str_replace(site,paste0(Gene,'-'),''))%>%
  dplyr::mutate(residue=stringr::str_replace_all(residue,"([STY])", ";\\1"))%>%
  dplyr::mutate(residue=stringr::str_replace(residue,"^;", ""))%>%
  dplyr::mutate(residue=stringr::str_replace_all(residue,"([sty])", ""))

####################### Getting KSEA Tables for every patient #######################

pvalue.cutoff <- 0.05

## There are two comparisons we will make, using the patients for which we have the appropriate data.
## We begin with the Early vs Pre-Treatment comparison. For this we have 6 patients we can use.
all.patient.ids <- unique(patients$ID)
pre.treatment.ids <- patients %>%
  dplyr::filter(Period == "Pre-Treatment") %>%
  dplyr::select(ID)

pre.treatment.ids <- pre.treatment.ids$ID

res <- list()

for (patient.id in pre.treatment.ids){
  print(patient.id)
  
  pretreatment <- patients %>%
    filter(ID == patient.id, Period == "Pre-Treatment")
  early <- patients %>%
    filter(ID == patient.id, Period == "Early")
  
  means.pretreatment <- phos.mat[, pretreatment$Sample]
  if (!is.null(dim(means.pretreatment))) {
    means.pretreatment <- apply(means.pretreatment, 1, FUN = mean, na.rm = T)
  }
  
  means.early <- phos.mat[, early$Sample]
  if (!is.null(dim(means.early))) {
    means.early <- apply(means.early, 1, FUN = mean, na.rm = T)
  }
  
  df.ksea <- data.frame(value = means.early - means.pretreatment, 
                        site = rownames(phos.mat), 
                        row.names = rownames(phos.mat)) %>%
    filter(!is.na(value))
  
  df.ksea <- arrange(df.ksea, desc(value))
  genelist = df.ksea$value
  names(genelist) = df.ksea$site
  print(head(genelist))
  genelist <- sort(genelist, decreasing = TRUE)
  
  
  ## Use cutoff pvalue of 1 to obtain full table
  prefix <- paste(patient.id, "KSEA Phosphosite Plus + NetworKIN")
  ksea <- df.ksea %>%      
    left_join(gene.to.site) %>%
    dplyr::select(Gene,Peptide, residue, value) %>%
    mutate(p_adj = 1) %>%
    computeKSEA(prefix = prefix, ksea_FDR = 1) 
  res[[prefix]] <- ksea
}

res.ksea.pre <- list()

## selecting only relevant columns

for (patient.id in pre.treatment.ids){
  prefix <- paste(patient.id, "KSEA Phosphosite Plus + NetworKIN")
  ksea <- res[[prefix]]
  
  x <- ksea %>%
    select(Kinase.Gene, z.score, p.value) %>%
    dplyr::rename(Kinase = Kinase.Gene, 
                  `z score` = z.score) %>%
    unique() %>%
    mutate(Database = "Phosphosite Plus + NetworKIN", ID = patient.id)
  rownames(x) <- x$Kinase
  
  res.ksea.pre[[patient.id]] <- x
}

save.image("KSEA enrichment Early vs Pre-Treatment Phosphosite Plus + NetworKIN.RData")

################### Combining the KSEA's ####################


IDs <- c("5087", "5104", "5145", "5174", "5210", "5180")


## The number of kinases taken from among the most enriched
number.chosen <- 50

kinases.present <- rownames(res.ksea.pre[["5087"]])
ksea.all <- cbind("5087" = res.ksea.pre[["5087"]][kinases.present, ],
                  "5104" = res.ksea.pre[["5104"]][kinases.present, ], 
                  "5145" = res.ksea.pre[["5145"]][kinases.present, ],
                  "5174" = res.ksea.pre[["5174"]][kinases.present, ], 
                  "5210" = res.ksea.pre[["5210"]][kinases.present, ],
                  "5180" = res.ksea.pre[["5180"]][kinases.present, ])
rownames(ksea.all) <- ksea.all$`5087.Kinase`
enrichment.cols <- which(grepl("z score", colnames(ksea.all)))
significance.cols <- which(grepl("p.value", colnames(ksea.all)))

ksea.enrichment <- ksea.all[, c(enrichment.cols)] %>%
  setNames(IDs) %>%
  as.matrix()
ksea.significance <- ksea.all[, c(significance.cols)] %>%
  setNames(IDs) %>%
  as.matrix()

## BH correction of p-values for better assessment of significance.
p.values <- as.list(ksea.significance) %>%
  p.adjust(method = "BH") %>%
  matrix(dimnames = list(rownames(ksea.significance), colnames(ksea.significance)),
         ncol = 6)
ksea.significance <- p.values

## Only the top 50 most enriched pathways are selected.
top.kinases <- apply(ksea.enrichment, 1, function(x) max(abs(x)))  %>%
  sort(decreasing = TRUE) %>%
  names() %>%
  .[1:number.chosen]
## Then select only significant pathways
sig.kinases <- which(apply(ksea.significance, 1, min) < pvalue.cutoff) %>%
  names()

ksea.significance <- ksea.significance[intersect(sig.kinases, top.kinases), ]
ksea.heatmap <- ksea.enrichment[intersect(sig.kinases, top.kinases), ]


####################### Heatmaps ###############################


library(gplots)

## Pathway names will wrap
rownames(ksea.heatmap) <- rownames(ksea.heatmap) %>%
  sapply(function(y) paste(strwrap(y, 38), collapse = "\n"), 
         USE.NAMES = FALSE)


## Adding indicator of significance for heatmap. * if significant, blank otherwise.
helper <- function(x) {
  out <- ""
  if (x < pvalue.cutoff){
    out <- "*"
  }
  return(out)
}

ksea.significance.cutoff <- lapply(ksea.significance, Vectorize(helper)) %>%
  matrix(ncol = 6)

png(filename = "Kinase heatmap Early vs Pre-Treatment - Phosphosite Plus + NetworKIN.png", 
    width = 1500, height = 2000)

p <- heatmap.2(ksea.heatmap, 
               cellnote = ksea.significance.cutoff,
               scale = 'none', 
               trace = 'none', 
               col = bluered, 
               notecol = "black",
               key = TRUE, 
               keysize = 1, 
               #( "bottom.margin", "left.margin", "top.margin", "right.margin" )
               key.par = list(mar = c(3,0,3,0), cex = 1.75),
               margins = c(8,47), 
               lhei = c(1,9), 
               lwid = c(1,7), 
               cexRow = 2,
               cexCol = 2.6,
               notecex = 2,
               hclustfun = function(x, ...) hclust(x, method = "ward.D", ...),
               distfun = function(x, ...) dist(x, method = "minkowski", p = 2, ...))
dev.off()


####################################### RESET ###################################

### Clears all objects in environment, as we compare Early to Late below.
rm(list = ls())

### Load code from the Setup section above.

####################### Getting KSEA Tables for every patient #######################

pvalue.cutoff <- 0.05

## There are two comparisons we will make, using the patients for which we have the appropriate data.
## Here we have the Late vs Late comparison. For this we have 2 patients we can use.
all.patient.ids <- unique(patients$ID)
late.treatment.ids <- patients %>%
  dplyr::filter(Period == "Late") %>%
  dplyr::select(ID)

late.treatment.ids <- late.treatment.ids$ID

res <- list()

for (patient.id in late.treatment.ids){
  print(patient.id)
  
  early <- patients %>%
    filter(ID == patient.id, Period == "Early")
  late <- patients %>%
    filter(ID == patient.id, Period == "Late")
  
  means.early <- phos.mat[, early$Sample]
  if (!is.null(dim(means.early))) {
    means.early <- apply(means.early, 1, FUN = mean, na.rm = T)
  }  
  
  means.late <- phos.mat[, late$Sample]
  if (!is.null(dim(means.late))) {
    means.late <- apply(means.late, 1, FUN = mean, na.rm = T)
  }
  
  df.ksea <- data.frame(value = means.late - means.early, 
                        site = rownames(phos.mat), 
                        row.names = rownames(phos.mat)) %>%
    filter(!is.na(value))
  
  df.ksea <- arrange(df.ksea, desc(value))
  genelist = df.ksea$value
  names(genelist) = df.ksea$site
  print(head(genelist))
  genelist <- sort(genelist, decreasing = TRUE)
  
  
  ## Use cutoff pvalue of 1 to obtain full table
  prefix <- paste(patient.id, "KSEA Phosphosite Plus + NetworKIN")
  ksea <- df.ksea %>%      
    left_join(gene.to.site) %>%
    dplyr::select(Gene,Peptide, residue, value) %>%
    mutate(p_adj = 1) %>%
    computeKSEA(prefix = prefix, ksea_FDR = 1) 
  res[[prefix]] <- ksea
}

res.ksea.late <- list()

## selecting only relevant columns

for (patient.id in late.treatment.ids){
  prefix <- paste(patient.id, "KSEA Phosphosite Plus + NetworKIN")
  ksea <- res[[prefix]]
  
  x <- ksea %>%
    select(Kinase.Gene, z.score, p.value) %>%
    dplyr::rename(Kinase = Kinase.Gene, 
                  `z score` = z.score) %>%
    unique() %>%
    mutate(Database = "Phosphosite Plus + NetowrKIN", ID = patient.id)
  rownames(x) <- x$Kinase
  
  res.ksea.late[[patient.id]] <- x
}

save.image("KSEA enrichment Late vs Early Phosphosite Plus + NetworKIN.RData")

################### Combining the KSEA's ####################


IDs <- c("5029", "5104")


## The number of kinases taken from among the most enriched
number.chosen <- 50

kinases.present <- rownames(res.ksea.late[["5029"]])
ksea.all <- cbind("5029" = res.ksea.late[["5029"]][kinases.present, ],
                  "5104" = res.ksea.late[["5104"]][kinases.present, ])
rownames(ksea.all) <- ksea.all$`5029.Kinase`
enrichment.cols <- which(grepl("z score", colnames(ksea.all)))
significance.cols <- which(grepl("p.value", colnames(ksea.all)))

ksea.enrichment <- ksea.all[, c(enrichment.cols)] %>%
  setNames(IDs) %>%
  as.matrix()
ksea.significance <- ksea.all[, c(significance.cols)] %>%
  setNames(IDs) %>%
  as.matrix()

## BH correction of p-values for better assessment of significance.
p.values <- as.list(ksea.significance) %>%
  p.adjust(method = "BH") %>%
  matrix(dimnames = list(rownames(ksea.significance), colnames(ksea.significance)),
         ncol = 2)
ksea.significance <- p.values

## Only the top 50 most enriched pathways are selected.
top.kinases <- apply(ksea.enrichment, 1, function(x) max(abs(x)))  %>%
  sort(decreasing = TRUE) %>%
  names() %>%
  .[1:number.chosen]
## Then select only significant pathways
sig.kinases <- which(apply(ksea.significance, 1, min) < pvalue.cutoff) %>%
  names()

ksea.significance <- ksea.significance[intersect(sig.kinases, top.kinases), ]
ksea.heatmap <- ksea.enrichment[intersect(sig.kinases, top.kinases), ]


####################### Heatmaps ###############################


library(gplots)

## Pathway names will wrap
rownames(ksea.heatmap) <- rownames(ksea.heatmap) %>%
  sapply(function(y) paste(strwrap(y, 38), collapse = "\n"), 
         USE.NAMES = FALSE)


## Adding indicator of significance for heatmap. * if significant, blank otherwise.
helper <- function(x) {
  out <- ""
  if (x < pvalue.cutoff){
    out <- "*"
  }
  return(out)
}

ksea.significance.cutoff <- lapply(ksea.significance, Vectorize(helper)) %>%
  matrix(ncol = 2)

png(filename = "Kinase heatmap Late vs Early - Phosphosite Plus + NetworKIN.png", 
    width = 1500, height = 2000)

p <- heatmap.2(ksea.heatmap, 
               cellnote = ksea.significance.cutoff,
               scale = 'none', 
               trace = 'none', 
               col = bluered, 
               notecol = "black",
               key = TRUE, 
               keysize = 1, 
               #( "bottom.margin", "left.margin", "top.margin", "right.margin" )
               key.par = list(mar = c(3,0,3,0), cex = 1.75),
               margins = c(8,47), 
               lhei = c(1,9), 
               lwid = c(1,7), 
               cexRow = 2,
               cexCol = 2.6,
               notecex = 2,
               hclustfun = function(x, ...) hclust(x, method = "ward.D", ...),
               distfun = function(x, ...) dist(x, method = "minkowski", p = 2, ...))
dev.off()



































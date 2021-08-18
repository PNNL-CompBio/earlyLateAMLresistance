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

############## This script generates the GSEA heatmaps seen in file 14.

## Saving plots to the folder in which this script is located!!
current_path <- getActiveDocumentContext()$path 
setwd(dirname(current_path))

source("../../../Util/synapseUtil.R")

########################### Setting up ##########################

# corrected for loading mass and plexID, long time and subtle
p.site.corrected <- querySynapseTable("syn26047149")

phos.dat <- p.site.corrected
phos.dat$Peptide <- sub("^.*@(.*)$", "\\1", phos.dat$Accession)
phos.dat$site <- sub("^(.*)@.*$", "\\1", phos.dat$Accession)
phos.dat$Gene <- sub("(.*)-.*$", "\\1", phos.dat$site)

## We have the metadata selected from the long form table
summary <- p.site.corrected[p.site.corrected$Accession == "AAK1-T606t@K.VQTT*PPPAVQGQK.V", 
                            c("Sample", "BeatAML Patient ID", 
                              "Plex", "Loading Mass", "Description", 
                              "Period", "Treatment")]

rownames(summary) <- summary$Sample

## Organizing time points into the right order (not alphanumeric)
time.points <- c("P2 Pre-Study", "P2 Cycle 1 Day 1", "P2 Cycle 1 Day 28", "P2 Cycle 2 Day 28", "P2 Cycle 3 Day 28", 
                 "P1B Cycle 1 Day 1", "P1B Cycle 1 Day 28", "P1B Cycle 2 Day 28", "P1B Cycle 3 Day 28", "P1B Cycle 6 Day 28",
                 "P1B Cycle 9 Day 28", "P1B Conc Cycle 9 Day 28") %>%
  factor(.)


summary$Description <- factor(summary$Description, levels = time.points)
summary <- summary[sort(summary$Sample), ]
summary$`BeatAML Patient ID` <- factor(summary$`BeatAML Patient ID`, levels = unique(summary$`BeatAML Patient ID`))
summary$Treatment <- make.names(summary$Treatment)

##select the samples of interest
healthy <- subset(summary, `BeatAML Patient ID` == 'Healthy Donor Stroma')

pre.treatment <- subset(summary, Period == "Pre-Treatment")

##then we spread the proteomics data into a matrix
phos.mat <- phos.dat%>%
  select(LogRatio, Sample, site)%>%
  tidyr::pivot_wider(values_from='LogRatio',names_from='Sample',
                     values_fn=list(LogRatio=mean))%>%
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
pre.treatment.ids <- summary %>%
  dplyr::filter(Period == "Pre-Treatment") %>%
  select(`BeatAML Patient ID`)

pre.treatment.ids <- as.character(pre.treatment.ids$`BeatAML Patient ID`)

res <- list()

for (patient.id in pre.treatment.ids){
  print(patient.id)
  
  pretreatment <- summary %>%
    filter(`BeatAML Patient ID` == patient.id, Treatment == "PRE")
  gilt <- summary %>%
    filter(`BeatAML Patient ID` == patient.id, Treatment == "GILT")
  
  means.pretreatment <- phos.mat[, pretreatment$Sample]
  if (!is.null(dim(means.pretreatment))) {
    means.pretreatment <- apply(means.pretreatment, 1, FUN = mean, na.rm = T)
  }
  
  means.gilt <- phos.mat[, gilt$Sample]
  if (!is.null(dim(means.gilt))) {
    means.gilt <- apply(means.gilt, 1, FUN = mean, na.rm = T)
  }
  
  df.ksea <- data.frame(value = means.gilt - means.pretreatment, 
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

save.image("KSEA enrichment Gilteritnib vs Pre-Treatment Phosphosite Plus + NetworKIN.RData")

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

png(filename = "Kinase heatmap Gilteritinib vs Pre-Treatment - Phosphosite Plus + NetworKIN.png", 
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

### Clears all objects in environment, as we compare Gilteritinib + Decitabine to Gilteritinib below.
rm(list = ls())

### Load code from the Setup section above.

####################### Getting KSEA Tables for every patient #######################

pvalue.cutoff <- 0.05

## There are two comparisons we will make, using the patients for which we have the appropriate data.
## Here we have the Late vs Late comparison. For this we have 2 patients we can use.
gilt.dec.ids <- summary %>%
  dplyr::filter(Treatment == "GILT...DEC") %>%
  dplyr::select(`BeatAML Patient ID`)

gilt.dec.ids <- as.character(gilt.dec.ids$`BeatAML Patient ID`) %>%
  unique()

IDs <- gilt.dec.ids

res <- list()

for (patient.id in gilt.dec.ids){
  print(patient.id)
  
  gilt <- summary %>%
    filter(`BeatAML Patient ID` == patient.id, Treatment == "GILT")
  gilt.dec <- summary %>%
    filter(`BeatAML Patient ID` == patient.id, Treatment == "GILT...DEC")
  
  means.gilt <- phos.mat[, gilt$Sample]
  if (!is.null(dim(means.gilt))) {
    means.gilt <- apply(means.gilt, 1, FUN = mean, na.rm = T)
  }
  
  means.gilt.dec <- phos.mat[, gilt.dec$Sample]
  if (!is.null(dim(means.gilt.dec))) {
    means.gilt.dec <- apply(means.gilt.dec, 1, FUN = mean, na.rm = T)
  }
  
  df.ksea <- data.frame(value = means.gilt.dec - means.gilt, 
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

for (patient.id in gilt.dec.ids){
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

save.image("KSEA enrichment Gilteritinib + Decitabine vs Gilteritinib Phosphosite Plus + NetworKIN.RData")

################### Combining the KSEA's ####################


IDs <- gilt.dec.ids


## The number of kinases taken from among the most enriched
number.chosen <- 50

kinases.present <- rownames(res.ksea.late[["5029"]])
ksea.all <- cbind("5029" = res.ksea.late[["5029"]][kinases.present, ],
                  "5104" = res.ksea.late[["5104"]][kinases.present, ],
                  "5174" = res.ksea.late[["5174"]][kinases.present, ])
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
         ncol = 3)
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
  matrix(ncol = 3)

png(filename = "Kinase heatmap Gilteritinib + Decitabine vs Gilteritinib - Phosphosite Plus + NetworKIN.png", 
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



































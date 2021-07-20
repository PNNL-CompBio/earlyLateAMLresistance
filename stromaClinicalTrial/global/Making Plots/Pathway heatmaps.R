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

############## This script generates the GSEA heatmaps seen in file 14.

########################### Setting up ##########################

## Saving plots to the folder in which this script is located!!
current_path <- getActiveDocumentContext()$path 
setwd(dirname(current_path))

# corrected for loading mass and plexID, long time and subtle
g.gene.corrected <- querySynapseTable("syn25706561")

prot.dat <- g.gene.corrected

# from here we can select just the metadata of interest
metadata.columns = c('Sample','BeatAML Patient ID','Plex','Loading Mass', 'Description')
summary <- prot.dat%>%
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
colnames(prot.dat)[[1]] <- "Gene"
prot.mat <- prot.dat%>%
  select(LogRatio,Sample,Gene)%>%
  tidyr::pivot_wider(values_from='LogRatio',names_from='Sample',
                     values_fn=list(LogRatio=mean),values_fill=0.0)%>%
  tibble::column_to_rownames('Gene')

### Pathway libraries

library(tidyr)
library(leapr)
data("krbpaths")

idx.kegg <- grepl("^KEGG_", krbpaths$names)
names.kegg <- krbpaths$names[idx.kegg]
names.kegg <- sub("KEGG_", "", names.kegg)
names.kegg <- gsub("_", " ", names.kegg)
desc.kegg <- krbpaths$desc[idx.kegg]
sizes.kegg <- krbpaths$sizes[idx.kegg]
Max <- max(sizes.kegg)
matrix.kegg <- krbpaths$matrix[idx.kegg, 1:Max]
keggpaths <- list(names = names.kegg,
                  desc = desc.kegg,
                  sizes = sizes.kegg,
                  matrix = matrix.kegg)

kegg.term.2.gene <- as.data.frame(keggpaths$matrix) %>%
  mutate(term = keggpaths$names) %>%
  pivot_longer(!term, names_to = "Column", values_to = "gene") %>%
  filter(!(gene == "null")) %>%
  select(term, gene)

kegg.term.2.name <- data.frame(term = keggpaths$names, name = keggpaths$names)

idx.reactome <- grepl("^REACTOME_", krbpaths$names)
names.reactome <- krbpaths$names[idx.reactome]
names.reactome <- sub("REACTOME_", "", names.reactome)
names.reactome <- gsub("_", " ", names.reactome)
desc.reactome <- krbpaths$desc[idx.reactome]
sizes.reactome <- krbpaths$sizes[idx.reactome]
Max <- max(sizes.reactome)
matrix.reactome <- krbpaths$matrix[idx.reactome, 1:Max]
reactomepaths <- list(names = names.reactome,
                      desc = desc.reactome,
                      sizes = sizes.reactome,
                      matrix = matrix.reactome)

reactome.term.2.gene <- as.data.frame(reactomepaths$matrix) %>%
  mutate(term = reactomepaths$names) %>%
  pivot_longer(!term, names_to = "Column", values_to = "gene") %>%
  filter(!(gene == "null")) %>%
  filter(!(gene == "")) %>%
  select(term, gene)

reactome.term.2.name <- data.frame(term = reactomepaths$names, name = reactomepaths$names)


####################### Getting GSEA Tables for every patient Early vs Pre-Treatment #######################

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
  
  means.pretreatment <- prot.mat[, pretreatment$Sample]
  if (!is.null(dim(means.pretreatment))) {
    means.pretreatment <- apply(means.pretreatment, 1, FUN = mean, na.rm = T)
  }
  
  means.early <- prot.mat[, early$Sample]
  if (!is.null(dim(means.early))) {
    means.early <- apply(means.early, 1, FUN = mean, na.rm = T)
  }
  
  df.gsea <- data.frame(value = means.early - means.pretreatment, 
                        Gene = rownames(prot.mat), 
                        row.names = rownames(prot.mat)) %>%
    filter(!is.na(value))
  
  df.gsea <- arrange(df.gsea, desc(value))
  genelist = df.gsea$value
  names(genelist) = df.gsea$Gene
  print(head(genelist))
  genelist <- sort(genelist, decreasing = TRUE)
  
  
  ## Use cutoff pvalue of 1 to obtain full table
  prefix <- paste(patient.id, "KEGG")
  gr <- clusterProfiler::GSEA(unlist(genelist), TERM2GENE = kegg.term.2.gene, eps = 1e-20,
                              TERM2NAME = kegg.term.2.name, pAdjustMethod = "BH", pvalueCutoff = 1,
                              nPermSimple = 2000)
  res[[prefix]] <- filter(as.data.frame(gr))
  
  prefix <- paste(patient.id, "REACTOME")
  gr <- clusterProfiler::GSEA(unlist(genelist), TERM2GENE = reactome.term.2.gene, eps = 1e-20,
                              TERM2NAME = reactome.term.2.name, pAdjustMethod = "BH", pvalueCutoff = 1,
                              nPermSimple = 2000)
  res[[prefix]] <- filter(as.data.frame(gr))
  
  prefix <- paste(patient.id, "Biological Process Ontology")
  res[[prefix]] <- plotOldGSEA(df.gsea, prefix, width = 13, gsea_FDR = 1, eps = 1e-20,
                               nPermSimple = 2000)
}

res.kegg.pre <- list()
res.reactome.pre <- list()
res.gobp.pre <- list()

## selecting only relevant columns

for (patient.id in pre.treatment.ids){
  prefix <- paste(patient.id, "KEGG")
  gsea <- res[[prefix]]
  
  x <- gsea %>%
    select(Description, NES, pvalue) %>%
    dplyr::rename(Pathway = Description, 
                  `Net Enrichment Score` = NES) %>%
    mutate(Database = "KEGG", ID = patient.id)
  
  res.kegg.pre[[patient.id]] <- x
  
  prefix <- paste(patient.id, "REACTOME")
  gsea <- res[[prefix]]
  
  x <- gsea %>%
    select(Description, NES, pvalue) %>%
    dplyr::rename(Pathway = Description, 
                  `Net Enrichment Score` = NES) %>%
    mutate(Database = "REACTOME", ID = patient.id)
  
  res.reactome.pre[[patient.id]] <- x
  
  prefix <- paste(patient.id, "Biological Process Ontology")
  gsea <- res[[prefix]]
  
  x <- gsea %>%
    select(Description, NES, pvalue) %>%
    dplyr::rename(Pathway = Description, 
                  `Net Enrichment Score` = NES) %>%
    mutate(Database = "Gene Ontology BP", ID = patient.id)
  
  res.gobp.pre[[patient.id]] <- x
}

save.image("GSEA enrichment - Early vs Pre-treatment for all patients.RData")

################### Combining the GSEA's Early vs Pre-Treatment####################


IDs <- c("5087", "5104", "5145", "5174", "5210", "5180")


##### KEGG

## The number of pathways taken from among the most enriched
number.chosen <- 50

pathways.present <- rownames(res.kegg.pre[["5087"]])
kegg.all <- cbind("5087" = res.kegg.pre[["5087"]][pathways.present, ],
                  "5104" = res.kegg.pre[["5104"]][pathways.present, ], 
                  "5145" = res.kegg.pre[["5145"]][pathways.present, ],
                  "5174" = res.kegg.pre[["5174"]][pathways.present, ], 
                  "5210" = res.kegg.pre[["5210"]][pathways.present, ],
                  "5180" = res.kegg.pre[["5180"]][pathways.present, ])
enrichment.cols <- which(grepl("Enrichment", colnames(kegg.all)))
significance.cols <- which(grepl("pvalue", colnames(kegg.all)))

kegg.enrichment <- kegg.all[, c(enrichment.cols)] %>%
  setNames(IDs) %>%
  as.matrix()
kegg.significance <- kegg.all[, c(significance.cols)] %>%
  setNames(IDs) %>%
  as.matrix()

## BH correction of p-values for better assessment of significance.
p.values <- as.list(kegg.significance) %>%
  p.adjust(method = "BH") %>%
  matrix(dimnames = list(rownames(kegg.significance), colnames(kegg.significance)),
         ncol = 6)
kegg.significance <- p.values

## Only the top 50 most enriched pathways are selected.
top.pathways <- apply(kegg.enrichment, 1, function(x) max(abs(x)))  %>%
  sort(decreasing = TRUE) %>%
  names() %>%
  .[1:number.chosen]
## Then select only significant pathways
sig.pathways <- which(apply(kegg.significance, 1, min) < pvalue.cutoff) %>%
  names()

kegg.significance <- kegg.significance[intersect(sig.pathways, top.pathways), ]
kegg.heatmap <- kegg.enrichment[intersect(sig.pathways, top.pathways), ]
  

##### REACTOME

## The number of pathways taken from among the most enriched
number.chosen <- 50


pathways.present <- rownames(res.reactome.pre[["5087"]])
reactome.all <- cbind("5087" = res.reactome.pre[["5087"]][pathways.present, ],
                  "5104" = res.reactome.pre[["5104"]][pathways.present, ], 
                  "5145" = res.reactome.pre[["5145"]][pathways.present, ],
                  "5174" = res.reactome.pre[["5174"]][pathways.present, ], 
                  "5210" = res.reactome.pre[["5210"]][pathways.present, ],
                  "5180" = res.reactome.pre[["5180"]][pathways.present, ])
enrichment.cols <- which(grepl("Enrichment", colnames(reactome.all)))
significance.cols <- which(grepl("pvalue", colnames(reactome.all)))

reactome.enrichment <- reactome.all[, c(enrichment.cols)] %>%
  setNames(IDs) %>%
  as.matrix()
reactome.significance <- reactome.all[, c(significance.cols)] %>%
  setNames(IDs) %>%
  as.matrix()

## BH correction of p-values for better assessment of significance.
p.values <- as.list(reactome.significance) %>%
  p.adjust(method = "BH") %>%
  matrix(dimnames = list(rownames(reactome.significance), colnames(reactome.significance)),
         ncol = 6)
reactome.significance <- p.values

## Only the top 50 most enriched pathways are selected.
top.pathways <- apply(reactome.enrichment, 1, function(x) max(abs(x)))  %>%
  sort(decreasing = TRUE) %>%
  names() %>%
  .[1:number.chosen]
## Then select only significant pathways
sig.pathways <- which(apply(reactome.significance, 1, min) < pvalue.cutoff) %>%
  names()

reactome.significance <- reactome.significance[intersect(sig.pathways, top.pathways), ]
reactome.heatmap <- reactome.enrichment[intersect(sig.pathways, top.pathways), ]


##### GENE ONTOLOGY BP

## The number of pathways taken from among the most enriched
number.chosen <- 50

pathways.present <- rownames(res.gobp.pre[["5087"]])
gobp.all <- cbind("5087" = res.gobp.pre[["5087"]][pathways.present, ],
                  "5104" = res.gobp.pre[["5104"]][pathways.present, ], 
                  "5145" = res.gobp.pre[["5145"]][pathways.present, ],
                  "5174" = res.gobp.pre[["5174"]][pathways.present, ], 
                  "5210" = res.gobp.pre[["5210"]][pathways.present, ],
                  "5180" = res.gobp.pre[["5180"]][pathways.present, ])
enrichment.cols <- which(grepl("Enrichment", colnames(gobp.all)))
significance.cols <- which(grepl("pvalue", colnames(gobp.all)))

gobp.enrichment <- gobp.all[, c(1, enrichment.cols)] %>%
  remove_rownames() %>%
  column_to_rownames(var = "5087.Pathway") %>%
  setNames(IDs) %>%
  as.matrix()
gobp.significance <- gobp.all[, c(1, significance.cols)] %>%
  remove_rownames() %>%
  column_to_rownames(var = "5087.Pathway") %>%
  setNames(IDs) %>%
  as.matrix()

## BH correction of p-values for better assessment of significance.
p.values <- as.list(gobp.significance) %>%
  p.adjust(method = "BH") %>%
  matrix(dimnames = list(rownames(gobp.significance), colnames(gobp.significance)),
         ncol = 6)
gobp.significance <- p.values

## Only the top 50 most enriched pathways are selected.
top.pathways <- apply(gobp.enrichment, 1, function(x) max(abs(x)))  %>%
  sort(decreasing = TRUE) %>%
  names() %>%
  .[1:number.chosen]
## Then select only significant pathways
sig.pathways <- which(apply(gobp.significance, 1, min) < pvalue.cutoff) %>%
  names()

gobp.significance <- gobp.significance[intersect(sig.pathways, top.pathways), ]
gobp.heatmap <- gobp.enrichment[intersect(sig.pathways, top.pathways), ]

####################### Heatmaps ###############################


library(gplots)

## Pathway names will wrap
rownames(kegg.heatmap) <- rownames(kegg.heatmap) %>%
  sapply(function(y) paste(strwrap(y, 38), collapse = "\n"), 
         USE.NAMES = FALSE)
rownames(reactome.heatmap) <- rownames(reactome.heatmap) %>%
  sapply(function(y) paste(strwrap(y, 50), collapse = "\n"), 
         USE.NAMES = FALSE)
rownames(gobp.heatmap) <- rownames(gobp.heatmap) %>%
  sapply(function(y) paste(strwrap(y, 50), collapse = "\n"), 
         USE.NAMES = FALSE)


## Adding indicator of significance for heatmap. * if significant, blank otherwise.
helper <- function(x) {
  out <- ""
  if (x < pvalue.cutoff){
    out <- "*"
  }
  return(out)
}

kegg.significance.cutoff <- lapply(kegg.significance, Vectorize(helper)) %>%
  matrix(ncol = 6)
reactome.significance.cutoff <- lapply(reactome.significance, Vectorize(helper)) %>%
  matrix(ncol = 6)
gobp.significance.cutoff <- lapply(gobp.significance, Vectorize(helper)) %>%
  matrix(ncol = 6)


## KEGG
png(filename = "Pathway heatmap Early vs Pre-Treatment - KEGG.png", width = 1380, height = 2500)
p <- heatmap.2(kegg.heatmap, 
               cellnote = kegg.significance.cutoff,
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


## REACTOME
png(filename = "Pathway heatmap Early vs Pre-Treatment - REACTOME.png", width = 1380, height = 2600)
p <- heatmap.2(reactome.heatmap, 
               cellnote = reactome.significance.cutoff,
               scale = 'none', 
               trace = 'none', 
               col = bluered, 
               notecol = "black",
               key = TRUE, 
               keysize = 1, 
               #( "bottom.margin", "left.margin", "top.margin", "right.margin" )
               key.par = list(mar = c(3,0,3,0), cex = 1.75),
               margins = c(8,55), 
               lhei = c(1,9), 
               lwid = c(1,7), 
               cexRow = 1.85,
               cexCol = 2.6,
               notecex = 2,
               hclustfun = function(x, ...) hclust(x, method = "ward.D", ...),
               distfun = function(x, ...) dist(x, method = "minkowski", p = 2, ...))
dev.off()


## GENE ONTOLOGY BIOLOGICAL PROCESS
png(filename = "Pathway heatmap Early vs Pre-Treatment - GENE ONTOLOGY.png", width = 1380, height = 2600)
p <- heatmap.2(gobp.heatmap, 
               cellnote = gobp.significance.cutoff,
               scale = 'none', 
               trace = 'none', 
               col = bluered, 
               notecol = "black",
               key = TRUE, 
               keysize = 1, 
               #( "bottom.margin", "left.margin", "top.margin", "right.margin" )
               key.par = list(mar = c(3,0,3,0), cex = 1.75),
               margins = c(8,55), 
               lhei = c(1,9), 
               lwid = c(1,7), 
               cexRow = 2,
               cexCol = 2.6,
               notecex = 2,
               hclustfun = function(x, ...) hclust(x, method = "ward.D", ...),
               distfun = function(x, ...) dist(x, method = "minkowski", p = 2, ...))
dev.off()


############################### RESET #############################################

### Clears all objects in environment, as we compare Early to Late below.
rm(list = ls())

### Load code from the Setup section above.


####################### Getting GSEA Tables for every patient Early vs Late #######################

pvalue.cutoff <- 0.05

## There are two comparisons we will make, using the patients for which we have the appropriate data.
## Here we have the Early vs Late comparison, for which we can use 2 patients.

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
  
  means.early <- prot.mat[, early$Sample]
  if (!is.null(dim(means.early))) {
    means.early <- apply(means.early, 1, FUN = mean, na.rm = T)
  }
  
  means.late <- prot.mat[, late$Sample]
  if (!is.null(dim(means.late))) {
    means.late <- apply(means.late, 1, FUN = mean, na.rm = T)
  }
  
  df.gsea <- data.frame(value = means.late - means.early, 
                        Gene = rownames(prot.mat), 
                        row.names = rownames(prot.mat)) %>%
    filter(!is.na(value))
  
  df.gsea <- arrange(df.gsea, desc(value))
  genelist = df.gsea$value
  names(genelist) = df.gsea$Gene
  print(head(genelist))
  genelist <- sort(genelist, decreasing = TRUE)
  
  ## Use cutoff pvalue of 1 to obtain full table
  prefix <- paste(patient.id, "KEGG")
  gr <- clusterProfiler::GSEA(unlist(genelist), TERM2GENE = kegg.term.2.gene, eps = 1e-20,
                              TERM2NAME = kegg.term.2.name, pAdjustMethod = "BH", pvalueCutoff = 1,
                              nPermSimple = 2000)
  res[[prefix]] <- filter(as.data.frame(gr))
  
  prefix <- paste(patient.id, "REACTOME")
  gr <- clusterProfiler::GSEA(unlist(genelist), TERM2GENE = reactome.term.2.gene, eps = 1e-20,
                              TERM2NAME = reactome.term.2.name, pAdjustMethod = "BH", pvalueCutoff = 1,
                              nPermSimple = 2000)
  res[[prefix]] <- filter(as.data.frame(gr))
  
  prefix <- paste(patient.id, "Biological Process Ontology")
  res[[prefix]] <- plotOldGSEA(df.gsea, prefix, width = 13, gsea_FDR = 1, eps = 1e-20,
                               nPermSimple = 2000)
}

res.kegg.late <- list()
res.reactome.late <- list()
res.gobp.late <- list()

## selecting only relevant columns

for (patient.id in late.treatment.ids){
  prefix <- paste(patient.id, "KEGG")
  gsea <- res[[prefix]]
  
  x <- gsea %>%
    select(Description, NES, pvalue) %>%
    dplyr::rename(Pathway = Description, 
                  `Net Enrichment Score` = NES) %>%
    mutate(Database = "KEGG", ID = patient.id)
  
  res.kegg.late[[patient.id]] <- x
  
  prefix <- paste(patient.id, "REACTOME")
  gsea <- res[[prefix]]
  
  x <- gsea %>%
    select(Description, NES, pvalue) %>%
    dplyr::rename(Pathway = Description, 
                  `Net Enrichment Score` = NES) %>%
    mutate(Database = "REACTOME", ID = patient.id)
  
  res.reactome.late[[patient.id]] <- x
  
  prefix <- paste(patient.id, "Biological Process Ontology")
  gsea <- res[[prefix]]
  
  x <- gsea %>%
    select(Description, NES, pvalue) %>%
    dplyr::rename(Pathway = Description, 
                  `Net Enrichment Score` = NES) %>%
    mutate(Database = "Gene Ontology BP", ID = patient.id)
  
  res.gobp.late[[patient.id]] <- x
}

save.image("GSEA enrichment - Late vs Early for all patients.RData")

################### Combining the GSEA's Early vs Late ####################


IDs <- c("5029", "5104")


##### KEGG

## The number of pathways taken from among the most enriched
number.chosen <- 50

pathways.present <- rownames(res.kegg.late[["5029"]])
kegg.all <- cbind("5029" = res.kegg.late[["5029"]][pathways.present, ],
                  "5104" = res.kegg.late[["5104"]][pathways.present, ])
enrichment.cols <- which(grepl("Enrichment", colnames(kegg.all)))
significance.cols <- which(grepl("pvalue", colnames(kegg.all)))

kegg.enrichment <- kegg.all[, c(enrichment.cols)] %>%
  setNames(IDs) %>%
  as.matrix()
kegg.significance <- kegg.all[, c(significance.cols)] %>%
  setNames(IDs) %>%
  as.matrix()

## BH correction of p-values for better assessment of significance.
p.values <- as.list(kegg.significance) %>%
  p.adjust(method = "BH") %>%
  matrix(dimnames = list(rownames(kegg.significance), colnames(kegg.significance)),
         ncol = 2)
kegg.significance <- p.values

## Only the top 50 most enriched pathways are selected.
top.pathways <- apply(kegg.enrichment, 1, function(x) max(abs(x)))  %>%
  sort(decreasing = TRUE) %>%
  names() %>%
  .[1:number.chosen]
## Then select only significant pathways
sig.pathways <- which(apply(kegg.significance, 1, min) < pvalue.cutoff) %>%
  names()

kegg.significance <- kegg.significance[intersect(sig.pathways, top.pathways), ]
kegg.heatmap <- kegg.enrichment[intersect(sig.pathways, top.pathways), ]


##### REACTOME

## The number of pathways taken from among the most enriched
number.chosen <- 50


pathways.present <- rownames(res.reactome.late[["5029"]])
reactome.all <- cbind("5029" = res.reactome.late[["5029"]][pathways.present, ],
                      "5104" = res.reactome.late[["5104"]][pathways.present, ])
enrichment.cols <- which(grepl("Enrichment", colnames(reactome.all)))
significance.cols <- which(grepl("pvalue", colnames(reactome.all)))

reactome.enrichment <- reactome.all[, c(enrichment.cols)] %>%
  setNames(IDs) %>%
  as.matrix()
reactome.significance <- reactome.all[, c(significance.cols)] %>%
  setNames(IDs) %>%
  as.matrix()

## BH correction of p-values for better assessment of significance.
p.values <- as.list(reactome.significance) %>%
  p.adjust(method = "BH") %>%
  matrix(dimnames = list(rownames(reactome.significance), colnames(reactome.significance)),
         ncol = 2)
reactome.significance <- p.values

## Only the top 50 most enriched pathways are selected.
top.pathways <- apply(reactome.enrichment, 1, function(x) max(abs(x)))  %>%
  sort(decreasing = TRUE) %>%
  names() %>%
  .[1:number.chosen]
## Then select only significant pathways
sig.pathways <- which(apply(reactome.significance, 1, min) < pvalue.cutoff) %>%
  names()

reactome.significance <- reactome.significance[intersect(sig.pathways, top.pathways), ]
reactome.heatmap <- reactome.enrichment[intersect(sig.pathways, top.pathways), ]


##### GENE ONTOLOGY BP

## The number of pathways taken from among the most enriched
number.chosen <- 50

pathways.present <- rownames(res.gobp.late[["5029"]])
gobp.all <- cbind("5029" = res.gobp.late[["5029"]][pathways.present, ],
                  "5104" = res.gobp.late[["5104"]][pathways.present, ])
enrichment.cols <- which(grepl("Enrichment", colnames(gobp.all)))
significance.cols <- which(grepl("pvalue", colnames(gobp.all)))

gobp.enrichment <- gobp.all[, c(1, enrichment.cols)] %>%
  remove_rownames() %>%
  column_to_rownames(var = "5029.Pathway") %>%
  setNames(IDs) %>%
  as.matrix()
gobp.significance <- gobp.all[, c(1, significance.cols)] %>%
  remove_rownames() %>%
  column_to_rownames(var = "5029.Pathway") %>%
  setNames(IDs) %>%
  as.matrix()

## BH correction of p-values for better assessment of significance.
p.values <- as.list(gobp.significance) %>%
  p.adjust(method = "BH") %>%
  matrix(dimnames = list(rownames(gobp.significance), colnames(gobp.significance)),
         ncol = 2)
gobp.significance <- p.values

## Only the top 50 most enriched pathways are selected.
top.pathways <- apply(gobp.enrichment, 1, function(x) max(abs(x)))  %>%
  sort(decreasing = TRUE) %>%
  names() %>%
  .[1:number.chosen]
## Then select only significant pathways
sig.pathways <- which(apply(gobp.significance, 1, min) < pvalue.cutoff) %>%
  names()

gobp.significance <- gobp.significance[intersect(sig.pathways, top.pathways), ]
gobp.heatmap <- gobp.enrichment[intersect(sig.pathways, top.pathways), ]


####################### Heatmaps ###############################


library(gplots)

## Pathway names will wrap
rownames(kegg.heatmap) <- rownames(kegg.heatmap) %>%
  sapply(function(y) paste(strwrap(y, 38), collapse = "\n"), 
         USE.NAMES = FALSE)
rownames(reactome.heatmap) <- rownames(reactome.heatmap) %>%
  sapply(function(y) paste(strwrap(y, 50), collapse = "\n"), 
         USE.NAMES = FALSE)
rownames(gobp.heatmap) <- rownames(gobp.heatmap) %>%
  sapply(function(y) paste(strwrap(y, 50), collapse = "\n"), 
         USE.NAMES = FALSE)


## Adding indicator of significance for heatmap. * if significant, blank otherwise.
helper <- function(x) {
  out <- ""
  if (x < pvalue.cutoff){
    out <- "*"
  }
  return(out)
}

kegg.significance.cutoff <- lapply(kegg.significance, Vectorize(helper)) %>%
  matrix(ncol = 2)
reactome.significance.cutoff <- lapply(reactome.significance, Vectorize(helper)) %>%
  matrix(ncol = 2)
gobp.significance.cutoff <- lapply(gobp.significance, Vectorize(helper)) %>%
  matrix(ncol = 2)


## KEGG
png(filename = "Pathway heatmap Early vs Late - KEGG.png", width = 1380, height = 2500)
p <- heatmap.2(kegg.heatmap, 
               cellnote = kegg.significance.cutoff,
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


## REACTOME
png(filename = "Pathway heatmap Early vs Late - REACTOME.png", width = 1380, height = 2600)
p <- heatmap.2(reactome.heatmap, 
               cellnote = reactome.significance.cutoff,
               scale = 'none', 
               trace = 'none', 
               col = bluered, 
               notecol = "black",
               key = TRUE, 
               keysize = 1, 
               #( "bottom.margin", "left.margin", "top.margin", "right.margin" )
               key.par = list(mar = c(3,0,3,0), cex = 1.75),
               margins = c(8,55), 
               lhei = c(1,9), 
               lwid = c(1,7), 
               cexRow = 1.85,
               cexCol = 2.6,
               notecex = 2,
               hclustfun = function(x, ...) hclust(x, method = "ward.D", ...),
               distfun = function(x, ...) dist(x, method = "minkowski", p = 2, ...))
dev.off()


## GENE ONTOLOGY BIOLOGICAL PROCESS
png(filename = "Pathway heatmap Early vs Late - GENE ONTOLOGY.png", width = 1380, height = 2600)
p <- heatmap.2(gobp.heatmap, 
               cellnote = gobp.significance.cutoff,
               scale = 'none', 
               trace = 'none', 
               col = bluered, 
               notecol = "black",
               key = TRUE, 
               keysize = 1, 
               #( "bottom.margin", "left.margin", "top.margin", "right.margin" )
               key.par = list(mar = c(3,0,3,0), cex = 1.75),
               margins = c(8,55), 
               lhei = c(1,9), 
               lwid = c(1,7), 
               cexRow = 2,
               cexCol = 2.6,
               notecex = 2,
               hclustfun = function(x, ...) hclust(x, method = "ward.D", ...),
               distfun = function(x, ...) dist(x, method = "minkowski", p = 2, ...))
dev.off()




































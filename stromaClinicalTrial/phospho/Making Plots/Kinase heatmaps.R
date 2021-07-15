library(amlresistancenetworks)
library(dplyr)
library(kableExtra)
library(tibble)
library(leapr)
library(ggplot2)
library(scales)
library(gridExtra)
library(rstudioapi)
source("../Util/synapseUtil.R")

############## This script generates the KSEA heatmaps seen in file 14.

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

# for exploration purposes, we create a new variable containing the type of sample,
# and another indicating the time group of the sample.
summary$Type <- case_when(summary$`BeatAML Patient ID` == "Healthy Donor Stroma" ~ "Healthy Donor Stroma",
                          summary$`BeatAML Patient ID` == "Cell line" ~ "Cell line")
summary[is.na(summary$Type), "Type"] <- "Treated"
summary$Period <- case_when(grepl("Pre",summary$Description) | grepl("Day 1", summary$Description) ~ "Early",
                            grepl("Day 28",summary$Description) ~ "Late")
rownames(summary) <- summary$Sample


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

all.patient.ids <- unique(patients$ID)
res <- list()

for (patient.id in all.patient.ids){
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
  prefix <- paste(patient.id, "KSEA Phosphosite Plus")
  ksea <- df.ksea %>%      
    left_join(gene.to.site) %>%
    dplyr::select(Gene,Peptide, residue, value) %>%
    mutate(p_adj = 1) %>%
    computeKSEA(prefix = prefix, ksea_FDR = 1) 
  res[[prefix]] <- ksea
}

res.ksea <- list()

## selecting only relevant columns

for (i in 1:7){
  patient.id <- all.patient.ids[i]
  
  prefix <- paste(patient.id, "KSEA Phosphosite Plus")
  ksea <- res[[prefix]]
  
  x <- ksea %>%
    select(Kinase.Gene, z.score, p.value) %>%
    dplyr::rename(Kinase = Kinase.Gene, 
                  `z score` = z.score) %>%
    unique() %>%
    mutate(Database = "Phosphosite Plus", ID = patient.id)
  rownames(x) <- x$Kinase
  
  res.ksea[[patient.id]] <- x
}

save.image("KSEA enrichment for all patients Phosphosite Plus.RData")

################### Combining the KSEA's ####################


IDs <- c("5029", "5087", "5104", "5145", "5174", "5210", "5180")


## The number of kinases taken from among the most enriched
number.chosen <- 50

kinases.present <- rownames(res.ksea[["5029"]])
ksea.all <- cbind("5029" = res.ksea[["5029"]][kinases.present, ], 
                  "5087" = res.ksea[["5087"]][kinases.present, ],
                  "5104" = res.ksea[["5104"]][kinases.present, ], 
                  "5145" = res.ksea[["5145"]][kinases.present, ],
                  "5174" = res.ksea[["5174"]][kinases.present, ], 
                  "5210" = res.ksea[["5210"]][kinases.present, ],
                  "5180" = res.ksea[["5180"]][kinases.present, ])
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
         ncol = 7)
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


## Adding indicator of significance for heatmap. T if significant, blank otherwise.
helper <- function(x) {
  out <- ""
  if (x < pvalue.cutoff){
    out <- "*"
  }
  return(out)
}

ksea.significance.cutoff <- lapply(ksea.significance, Vectorize(helper)) %>%
  matrix(ncol = 7)

png(filename = "Kinase heatmap all patients cutoff at 0.05 - Phosphosite Plus.png", 
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
               distfun = function(x, ...) dist(x, method = "minkowski", p = 3, ...))
dev.off()



































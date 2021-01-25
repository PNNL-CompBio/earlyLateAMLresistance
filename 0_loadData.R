#' loadTreatedData
#' 
#' 
#' 
#
library(amlresistancenetworks)
require(readxl)
require(dplyr)
require(tidyr)
global.syn<-'syn24181246'
phos.syn<- 'syn24181242'
phos.stoich <-'syn24181286'

metadata.syn <- 'syn24185272'

syn <- synapseLogin()


#' parse the meta data, which i updated manually before upload
#' @param metadata.syn: synapseID of file
#' @return data frame
parseMetadata<-function(metadata.syn){
  tab <- readxl::read_xlsx(syn$get(metadata.syn)$path)%>%
    dplyr::select(cellLine='Sample description',Ligand,Treatment, 
                  'Time (minutes)', Sample='SampleName')%>%
    distinct()
  return(tab)
}

#' parse the proteomics matrix
#' @param global.syn synapse id of global crosstab file
#' @param metadata.syn synapse id of metadata
#' @return data.frame
parseProt<-function(global.syn, metadata.syn){
  mat <- read.csv(syn$get(global.syn)$path,sep='\t',header=T)
  
  mat<-mat%>%tidyr::pivot_longer(cols=c(5:ncol(mat)),names_to='Sample',
                        values_to='LogRatio')%>%
    left_join(parseMetadata(metadata.syn))%>%
        subset(!is.na(LogRatio))

  #  tidyr::replace_na(list(LogRatio=-1000000))
  return(mat)
  
}

#' parse the phosphoproteomics matrix
#' @param phos.syn synapse id of phospho crosstab file
#' @param metadata.syn synapse id of metadata
parsePhos<-function(phos.syn, metadata.syn){
   mat <- read.csv(syn$get(phos.syn)$path,sep='\t',header=T)
   
  mat<-mat%>%tidyr::pivot_longer(cols=c(6:ncol(mat)),names_to='Sample',
                        values_to='LogRatio')%>%
     left_join(parseMetadata(metadata.syn))%>%
    subset(!is.na(LogRatio))
    #tidyr::replace_na(list(LogRatio=-1000000))

}

##first grab the files and merge them with their metadata
prot<-parseProt(global.syn, metadata.syn)
phos<-parsePhos(phos.syn,metadata.syn)
stoich.phos<-parsePhos(phos.stoich,metadata.syn)


##now store the tables to synapse!!
p1<-synTableStore(prot,'Gilteritinib treated MOLM14 Proteomics')

p2<-synTableStore(phos,'Gilteritinib treated MOLM14 Phosphoproteomics')
p3<-synTableStore(stoich.phos,'Gilteritinib treated MOLM14 Site-corrected Phosphoproteomics')

##todo: set permissions YOU NEED TO DO THIS MANUALLY RIGHT NOW!!!

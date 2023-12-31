% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MMRFVariant_functions.R
\name{VariantFinder_PlotVariantsbyGene}
\alias{VariantFinder_PlotVariantsbyGene}
\title{VariantFinder_PlotVariantsbyGene}
\usage{
VariantFinder_PlotVariantsbyGene(
  variant.ann,
  ListGene = NULL,
  topN = 20,
  filenm = "PlotVariantsbyGene_heatmap",
  height = 10,
  width = 10
)
}
\arguments{
\item{variant.ann}{is the dataframe of annotated variants downloaded from MMRF-Commpass Researcher Gateway 
(i.e. MMRF_CoMMpass_IA14a_All_Canonical_Variants file) and imported into environment}

\item{ListGene}{is the list of the genes to analyze (optional).}

\item{topN}{is the top number of genes whose to visualize the variants occurrence number}

\item{filenm}{is the name of the png file. If filenm is Null, the plot is draw but it is not saved.}

\item{height}{Image height}

\item{width}{Image width}
}
\value{
heatmap of the annotated variants occurrences
}
\description{
draws heatmap of the annotated variants occurrences
}
\examples{
variant.ann<- data.frame(public_id=c("MMRF_0000","MMRF_0001",
                                     "MMRF_0002","MMRF_0003",
                                     "MMRF_0004","MMRF_0005",
                                     "MMRF_0006","MMRF_0007",
                                     "MMRF_0008",""),                  
                 dbSNP=c(rep("rs755588843",2),rep("rs569344016",5),rep("rs2066497",2),rep(".",1)),                                                    
                 Effect=c(rep("intragenic_variant",3),
                           rep("missense_variant",2),
                           rep("intron_variant",1),
                           rep("5_prime_UTR_variant",4)),
                  Gene=c(rep("PRDM16",3),
                           rep("AGO1",2),
                           rep("FPGT-TNNI3K",1),
                           rep("TNNI3K",4)), 
                 REF=c(rep("C",3),
                           rep("G",2),
                           rep("A",1),
                           rep("T",4)),                            
                 ALT=c(rep("GGCCT",3),
                           rep("G",2),
                           rep("T",1),
                           rep("A",4)),    
                  Biotype=c(rep("protein_coding",6),
                           rep("antisense",2),
                           rep("processed_pseudogene",2)),   
                           
                 Impact= c(rep("MODERATE",2),rep("MODIFIER",2),
                            rep("LOW",3),rep("HIGH",2),rep("MODIFIER",1)),
                 feature_type= c(rep("ENST00000388718",2),rep("ENST00000344616",2),
                            rep("ENST00000431492",3),rep("ENST00000390268",2),rep("ENST00000316407",1)),
                            
                 SIFT= c(rep("0.035,0.035,0.057,0.057,0.035,0.042,0.04,0.058",2),rep("0.002,0.002,0.001,0.002",2),
                            rep("0.614,0.614,0.781",6)),           
                            
                 Polyphen2=c(rep("0.021,0.986,0.884,0.977",2),rep("0.99",2),
                            rep("0.614,0.781",6))           
                                                              
                                 
 )





variants.plot<-VariantFinder_PlotVariantsbyGene(variant.ann,ListGene,topN=10,)
variants.plot<-VariantFinder_PlotVariantsbyGene(variant.ann,topN=10)
}

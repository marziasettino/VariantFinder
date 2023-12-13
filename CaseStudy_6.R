library(VariantFinder)
library(dplyr)
library(DT)
library(ggplot2)
library(stringr)
library(ggpubr)
library(survminer)
library(survival)
library(formattable)



#-----------dataset from MMRF-RG---------

#patient<-MMRF_CoMMpass_IA15_PER_PATIENT
#trt<-MMRF_CoMMpass_IA15_STAND_ALONE_TRTRESP
#variant.ann<-MMRF_CoMMpass_IA15a_All_Canonical_Variants

#--------Example-----------
#patient<-patient.example
#trt<-trt.example
#variant.ann<-variant.ann.example


#(GRCh37)

#hg19





# Analyze the topN recurrent variants in the complete samples cohort 


## Heatmap of the N# of variants occurrence in the complete samples cohort



variants.plot.all<-VariantFinder_PlotVariantsbyGene(variant.ann,height=10, width=20,topN=10,
                                                  
                                                  filenm="PlotVariantsbyGene_heatmapAll")




  
  
  ## Get the list of variants found in the in the complete samples cohort ranked by the occurrence number
  
 
ListSNPs<-VariantFinder_GetVariantsbyGene(variant.ann)

ListSNPs<-ListSNPs[order(ListSNPs$count, decreasing = TRUE),]




datatable(head(as.data.frame(ListSNPs),10), 
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)



impact.table<-VariantFinder_GetImpact(variant.ann,ListSNPs$dbSNP)

impact.table<-dplyr::select(impact.table,dbSNP,Gene,REF,ALT,feature,Effect,
                                    SIFT_Impact,Polyphen_Impact,Impact)

  
  ## Perform the Impact-Effect plot
  

plot.impact.effect.all<-VariantFinder_PlotbyEffectImpact(variant.ann,ListSNPs$dbSNP,topN=3,height=20, 
                                                       width=30, filenm="PlotbyEffectImpactAll")

  
  # Prioritizing variants in a known Multiple Myeloma (MM) disease-causing gene set
  
# The six gene-signature including the genes KRAS, NRAS, TP53, FAM46C, DIS3, BRAF have a high recurrence rate and may play important roles in the pathogenesis, progression and prognosis of MM. 
#This case study shows how VariantFinder performs the integrative analysis of variants in that six gene-signature to prioritize pathogenic variants involved in the MM.


## Set the list of gene "ListGene" to explore 



ListGene<-c("KRAS", "NRAS","TP53","FAM46C","DIS3","BRAF")


## Draw the Heatmap of the N# of variants occurring in "ListGene" 




variants.plot<-VariantFinder_PlotVariantsbyGene(variant.ann,ListGene,height=15,
                                              width=20,topN=50,
                                              filenm="PlotVariantsbyGene_heatmap")








## Get the list of Variants found in the gene set <ListGene>


ListSNPs<-VariantFinder_GetVariantsbyGene(variant.ann, ListGene)

#head(ListSNPs,20)

## Plot the Impact-Effect of Variants 

plot.impact.effect<-VariantFinder_PlotbyEffectImpact(variant.ann,ListSNPs$dbSNP,topN=50,height=30, 
                                                   width=15, filenm="PlotbyEffectImpact")




## Perform the impact table (ordered by ascending SIFT and descending Poliphen-2)


impact.table<-VariantFinder_GetImpact(variant.ann,ListSNPs$dbSNP)

#For semplification purposes, we visualize a subset of columns and rows




impact.table.sub<-dplyr::select(impact.table,dbSNP,Gene,REF,ALT,feature,Effect,
                                SIFT_Impact,Polyphen_Impact,Impact)

head(unique(impact.table.sub),10)





ListSNPs_NRAS<-VariantFinder_GetVariantsbyGene(variant.ann,"NRAS")



datatable(as.data.frame(ListSNPs_NRAS), 
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)




#We take into account only the top five SNPs by occurrence (see <KRAS_SNPs.tab>)



impact.table_NRAS<-VariantFinder_GetImpact(variant.ann,ListSNPs_NRAS$dbSNP)

#For semplification purposes, we visualize a subset of columns and rows


impact.table_NRAS<-dplyr::select(impact.table_NRAS,dbSNP,Gene,REF,ALT,
                                 feature,Effect,SIFT_Impact,Polyphen_Impact,Impact)

head(unique(impact.table_NRAS),10)  



## Plot the Impact-Effect of SNPs in NRAS gene



plot.impact.effect_NRAS<-VariantFinder_PlotbyEffectImpact(variant.ann,ListSNPs_NRAS$dbSNP,topN=50,height=30, 
                                                        width=15, filenm="PlotbyEffectImpact_NRAS")






## KM Survival curves in NRAS gene
#(!) SNPs in <ListSNPs_NRAS> are discarded if:
#  a) only a group with respect to FilterBy parameter is found
#b) pvalue is>=0.05



NRAS_SNPs.treat<-c("rs11554290","rs121913254","rs121913237", "rs121913255")

NRAS_surv.treatment<-VariantFinder_SurvivalKM(patient,  
                                            trt,
                                            variant.ann,
                                            NRAS_SNPs.treat,
                                            FilterBy="Treatment", 
                                            filename="KM_Plot_NRAS_treatment",
                                            xlim = c(100,3000),
                                            height=22,
                                            width=12,
                                            conf.range = FALSE,
                                            color = c("Dark2"))








NRAS_surv.Effect<-VariantFinder_SurvivalKM(patient,  #no significant results are found (all pvalue>0.05)
                                         trt,
                                         variant.ann,
                                         ListSNPs_NRAS,
                                         FilterBy="Effect", 
                                         filename="KM_Plot_NRAS_effect",
                                         xlim = c(100,100),
                                         height=22,
                                         width=12,
                                         conf.range = FALSE,
                                         color = c("Dark2"))


# see (*)

NRAS_SNPs.stage<-c("rs121913254")
NRAS_surv.Stage<-VariantFinder_SurvivalKM(patient,  
                                        trt,
                                        variant.ann,
                                        NRAS_SNPs.stage,
                                        FilterBy="Stage", 
                                        filename="KM_Plot_NRAS_stage",
                                        xlim = c(100,3000),
                                        height=22,
                                        width=12,
                                        conf.range = FALSE,
                                        color = c("Dark2"))






# see (*)
NRAS_SNPs.bestresp<-c("rs11554290","rs121913254","rs121434595", "rs121913237")
NRAS_surv.Bestresp<-VariantFinder_SurvivalKM(patient,  
                                           trt,
                                           variant.ann,
                                           NRAS_SNPs.bestresp,
                                           FilterBy="Bestresp", 
                                           filename="KM_Plot_NRAS_bestresp",
                                           xlim = c(100,3000),
                                           height=22,
                                           width=12,
                                           conf.range = FALSE,
                                           color = c("Dark2"))







# see (*)

NRAS_surv.Gender<-VariantFinder_SurvivalKM(patient,  #All SNPs have pvalue<=0.05
                                         trt,
                                         variant.ann,
                                         ListSNPs_NRAS,
                                         FilterBy="Gender", 
                                         filename="KM_Plot_NRAS_gender",
                                         xlim = c(100,3000),
                                         height=22,
                                         width=12,
                                         conf.range = FALSE,
                                         color = c("Dark2"))




# see (*)

NRAS_surv.Biotype<-VariantFinder_SurvivalKM(patient,  #All SNPs have have only a group with respect to FilterBy parameter 
                                          trt,
                                          variant.ann,
                                          ListSNPs_NRAS,
                                          FilterBy="Biotype", 
                                          filename="KM_Plot_NRAS_biotype",
                                          xlim = c(100,3000),
                                          height=22,
                                          width=12,
                                          conf.range = FALSE,
                                          color = c("Dark2"))



# see (*)
NRAS_SNPs.ethnicity<-c("rs11554290","rs121913254","rs121913237")
NRAS_surv.Ethnicity<-VariantFinder_SurvivalKM(patient,  
                                            trt,
                                            variant.ann,
                                            NRAS_SNPs.ethnicity,
                                            FilterBy="Ethnicity", 
                                            filename="KM_Plot_NRAS_ethnicity",
                                            xlim = c(100,3000),
                                            height=22,
                                            width=12,
                                            conf.range = FALSE,
                                            color = c("Dark2"))








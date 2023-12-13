## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  warning=FALSE, 
  message=FALSE,
  comment = "#>"
)

## ---- echo = FALSE,hide=TRUE, message=FALSE,warning=FALSE---------------------
devtools::load_all(".")

## ----message=FALSE, warning=FALSE, include=FALSE------------------------------
library(dplyr)
library(DT)
library(ggplot2)
library(stringr)
library(ggpubr)

## ----results = 'hide', message=FALSE, warning=FALSE, eval = F-----------------
#  patient<-MMRF_CoMMpass_IA15_PER_PATIENT
#  
#  trt<-MMRF_CoMMpass_IA15_STAND_ALONE_TRTRESP
#  
#  variant.ann<-MMRF_CoMMpass_IA15a_All_Canonical_Variants
#  

## ----results = 'hide', message=FALSE, warning=FALSE, eval = F-----------------
#  variants.plot.all<-MMRFVariant_PlotVariantsbyGene(variant.ann,height=10, width=20,topN=10,
#  
#                                                    filenm="PlotVariantsbyGene_heatmapAll")
#  

## ----figurename="PlotVariantsbyGene_heatmap", echo=FALSE, fig.cap="Heatmap of the N# of variants occurrence in the complete cohort", out.width = '90%'----
knitr::include_graphics("imgs/PlotVariantsbyGene_heatmap_all.png")

## ----results = 'hide', message=FALSE, warning=FALSE, eval = F-----------------
#  ListSNPs.all<-MMRFVariant_GetVariantsbyGene(variant.ann)
#  
#  
#  ListSNPs.all<-ListSNPs.all[order(ListSNPs.all$count, decreasing = TRUE),]
#  
#  ListSNPs.all.10<-head(unique(ListSNPs.all),10)
#  

## ----echo=TRUE, message=FALSE, warning=FALSE----------------------------------
datatable(as.data.frame(ListSNPs.all.10), 
              options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
              rownames = FALSE)




## ----results = 'hide', message=FALSE, warning=FALSE, eval = F-----------------
#  
#  impact.table.all<-MMRFVariant_GetImpact(variant.ann,ListSNPs.all$dbSNP)
#  
#  impact.table.all.sub<-dplyr::select(impact.table.all,dbSNP,Gene,REF,ALT,feature,Effect,
#                                  SIFT_Impact,Polyphen_Impact,Impact)
#  
#  head(unique(impact.table.all.sub),10)
#  
#  

## ----figurename="ImpactTableAll", echo=FALSE, fig.cap="Impact table of each SNP in ListGene", out.width = '90%'----
knitr::include_graphics("imgs/ImpactTableAll.png")

## ----results = 'hide', message=FALSE, warning=FALSE, eval = F-----------------
#  
#  ListSNPs.all.100<-head(ListSNPs.bycount.all$dbSNP,100)
#  
#  plot.impact.effect.all<-MMRFVariant_PlotbyEffectImpact(variant.ann,ListSNPs.all.100,topN=3,height=20,
#                                                     width=30, filenm="PlotbyEffectImpactAll")
#  

## ----figurename="ImpactTableAll", echo=FALSE, fig.cap="Impact table of each SNP in ListGene", out.width = '90%'----
knitr::include_graphics("imgs/PlotbyEffectImpactAll.png")

## ----figurename="workflow", echo=FALSE, fig.cap="The workflow describes graphically step by step the procedure carried out to perform this case of study ", out.width = '99%'----
knitr::include_graphics("imgs/workflow.png")

## ----results = 'hide', message=FALSE, warning=FALSE, eval = F-----------------
#  ListGene<-c("KRAS", "NRAS","TP53","FAM46C","DIS3","BRAF")
#  
#  

## ----results = 'hide', message=FALSE, warning=FALSE, eval = F-----------------
#  
#  variants.plot<-MMRFVariant_PlotVariantsbyGene(variant.ann,ListGene,height=15,
#                                                width=20,topN=50,
#                                                filenm="PlotVariantsbyGene_heatmap")
#  
#  
#  
#  
#  
#  
#  

## ----figurename="workflow", echo=FALSE, fig.cap="The workflow describes graphically step by step the procedure carried out to perform this case of study ", out.width = '99%'----
knitr::include_graphics("imgs/PlotVariantsbyGene_heatmap.png")

## ----results = 'hide', message=FALSE, warning=FALSE, eval = F-----------------
#  
#  ListSNPs<-MMRFVariant_GetVariantsbyGene(variant.ann, ListGene)
#  
#  head(ListSNPs,20)

## ----echo=TRUE, message=FALSE, warning=FALSE----------------------------------
datatable(as.data.frame(ListSNPs), 
              options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
              rownames = FALSE)




## ----results = 'hide', message=FALSE, warning=FALSE, eval = F-----------------
#  plot.impact.effect<-MMRFVariant_PlotbyEffectImpact(variant.ann,ListSNPs,topN=50,height=30,
#                                                     width=15, filenm="PlotbyEffectImpact")
#  

## ----figurename="workflow", echo=FALSE, fig.cap="MMRFVariant_PlotbyEffectImpact plot ", out.width = '99%'----
knitr::include_graphics("imgs/PlotbyEffectImpact.png")

## ----results = 'hide', message=FALSE, warning=FALSE, eval = F-----------------
#  impact.table<-MMRFVariant_GetImpact(variant.ann,ListSNPs)
#  

## ----results = 'hide', message=FALSE, warning=FALSE, eval = F-----------------
#  
#  
#  impact.table.sub<-dplyr::select(impact.table,dbSNP,Gene,REF,ALT,feature,Effect,
#                                  SIFT_Impact,Polyphen_Impact,Impact)
#  
#  head(unique(impact.table.sub),10)
#  
#  

## ----figurename="ImpactTable", echo=FALSE, fig.cap="Impact table of each SNP in ListGene", out.width = '90%'----
knitr::include_graphics("imgs/ImpactTable.png")

## ----results = 'hide', message=FALSE, warning=FALSE, eval = F-----------------
#  
#  
#  ListSNPs_NRAS<-MMRFVariant_GetVariantsbyGene(variant.ann,"NRAS")
#  
#  

## ----echo=TRUE, message=FALSE, warning=FALSE----------------------------------
datatable(as.data.frame(ListSNPs_NRAS), 
              options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
              rownames = FALSE)




## ----results = 'hide', message=FALSE, warning=FALSE, eval = F-----------------
#  
#  
#  impact.table_NRAS<-MMRFVariant_GetImpact(variant.ann,ListSNPs_NRAS$dbSNP)

## ----results = 'hide', message=FALSE, warning=FALSE, eval = F-----------------
#  impact.table_NRAS<-dplyr::select(impact.table_NRAS,dbSNP,Gene,REF,ALT,
#                                  feature,Effect,SIFT_Impact,Polyphen_Impact,Impact)
#  
#  head(unique(impact.table_NRAS),10)
#  
#  

## ----figurename="ImpactTable_NRAS", echo=FALSE, fig.cap="Impact table of each SNP in NRAS", out.width = '90%'----
knitr::include_graphics("imgs/ImpactTable_NRAS.png")

## ----results = 'hide', message=FALSE, warning=FALSE, eval = F-----------------
#  
#  
#  
#  plot.impact.effect_NRAS<-MMRFVariant_PlotbyEffectImpact(variant.ann,ListSNPs_NRAS,topN=50,height=30,
#                                                           width=15, filenm="PlotbyEffectImpact_NRAS")
#  
#  

## ----figurename="ImpactTable", echo=FALSE, fig.cap="Impact table of each SNP in ListGene", out.width = '90%'----
knitr::include_graphics("imgs/PlotbyEffectImpact_NRAS.png")

## ----results = 'hide', message=FALSE, warning=FALSE, eval = F-----------------
#  
#  NRAS_SNPs.treat<-c("rs11554290","rs121913254","rs121913237", "rs121913255")
#  
#  NRAS_surv.treatment<-MMRFVariant_SurvivalKM(patient,
#                                              trt,
#                                              variant.ann,
#                                              NRAS_SNPs.treat,
#                                              FilterBy="Treatment",
#                                              filename="KM_Plot_NRAS_treatment",
#                                              xlim = c(100,3000),
#                                              height=22,
#                                              width=12,
#                                              conf.range = FALSE,
#                                              color = c("Dark2"))
#  
#  
#  
#  
#  
#  
#  
#  
#  NRAS_surv.Effect<-MMRFVariant_SurvivalKM(patient,  #no significant results are found (all pvalue>0.05)
#                                           trt,
#                                           variant.ann,
#                                           ListSNPs_NRAS,
#                                           FilterBy="Effect",
#                                           filename="KM_Plot_NRAS_effect",
#                                           xlim = c(100,100),
#                                           height=22,
#                                           width=12,
#                                           conf.range = FALSE,
#                                           color = c("Dark2"))
#  
#  
#  # see (*)
#  
#  NRAS_SNPs.stage<-c("rs121913254")
#  NRAS_surv.Stage<-MMRFVariant_SurvivalKM(patient,
#                                          trt,
#                                          variant.ann,
#                                          NRAS_SNPs.stage,
#                                          FilterBy="Stage",
#                                          filename="KM_Plot_NRAS_stage",
#                                          xlim = c(100,3000),
#                                          height=22,
#                                          width=12,
#                                          conf.range = FALSE,
#                                          color = c("Dark2"))
#  
#  
#  
#  
#  
#  
#  # see (*)
#  NRAS_SNPs.bestresp<-c("rs11554290","rs121913254","rs121434595", "rs121913237")
#  NRAS_surv.Bestresp<-MMRFVariant_SurvivalKM(patient,
#                                             trt,
#                                             variant.ann,
#                                             NRAS_SNPs.bestresp,
#                                             FilterBy="Bestresp",
#                                             filename="KM_Plot_NRAS_bestresp",
#                                             xlim = c(100,3000),
#                                             height=22,
#                                             width=12,
#                                             conf.range = FALSE,
#                                             color = c("Dark2"))
#  
#  
#  
#  
#  
#  
#  
#  # see (*)
#  
#  NRAS_surv.Gender<-MMRFVariant_SurvivalKM(patient,  #All SNPs have pvalue<=0.05
#                                           trt,
#                                           variant.ann,
#                                           ListSNPs_NRAS,
#                                           FilterBy="Gender",
#                                           filename="KM_Plot_NRAS_gender",
#                                           xlim = c(100,3000),
#                                           height=22,
#                                           width=12,
#                                           conf.range = FALSE,
#                                           color = c("Dark2"))
#  
#  
#  
#  
#  # see (*)
#  
#  NRAS_surv.Biotype<-MMRFVariant_SurvivalKM(patient,  #All SNPs have have only a group with respect to FilterBy parameter
#                                            trt,
#                                            variant.ann,
#                                            ListSNPs_NRAS,
#                                            FilterBy="Biotype",
#                                            filename="KM_Plot_NRAS_biotype",
#                                            xlim = c(100,3000),
#                                            height=22,
#                                            width=12,
#                                            conf.range = FALSE,
#                                            color = c("Dark2"))
#  
#  
#  
#  # see (*)
#  NRAS_SNPs.ethnicity<-c("rs11554290","rs121913254","rs121913237")
#  NRAS_surv.Ethnicity<-MMRFVariant_SurvivalKM(patient,
#                                              trt,
#                                              variant.ann,
#                                              NRAS_SNPs.ethnicity,
#                                              FilterBy="Ethnicity",
#                                              filename="KM_Plot_NRAS_ethnicity",
#                                              xlim = c(100,3000),
#                                              height=22,
#                                              width=12,
#                                              conf.range = FALSE,
#                                              color = c("Dark2"))
#  
#  
#  
#  
#  
#  
#  
#  
#  
#  
#  

## ----figurename="ImpactTable", echo=FALSE, fig.cap="KM Survival curves in NRAS gene by ethnicity", out.width = '90%'----
knitr::include_graphics("imgs/KM_Surv_ethnicity.png")

## ----figurename="ImpactTable", echo=FALSE, fig.cap="KM Survival curves in NRAS gene by bestresp", out.width = '90%'----
knitr::include_graphics("imgs/KM_Surv_bestresp.png")

## ----figurename="ImpactTable", echo=FALSE, fig.cap="KM Survival curves in NRAS gene by stage", out.width = '90%'----
knitr::include_graphics("imgs/KM_Surv_stage.png")

## ----figurename="ImpactTable", echo=FALSE, fig.cap="KM Survival curves in NRAS gene by treatment", out.width = '90%'----
knitr::include_graphics("imgs/KM_Surv_treatment.png")


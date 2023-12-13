## ----results = 'hide', message=FALSE, warning=FALSE, eval = F-----------------
#  NRAS_SNPs<-head(ListSNPs_NRAS,5) #["rs2157615","rs61731685","rs370560636","rs116293337"]
#  
#  surv<-MMRFVariant_SurvivalKM(patient,
#                                trt,
#                                variant.ann,
#                                NRAS_SNPs,
#                                FilterBy="Treatment",
#                                filename=NULL,
#                                xlim = c(100,3000),
#                                conf.range = FALSE,
#                                color = c("Dark2"))
#  
#  
#  
#  

## ----figurename2, echo=FALSE, fig.cap="KM survival curves are drawn for each SNP found in the NRAS gene in the case of the FilterBy parameter is set to Treatment", out.width = '99%'----
knitr::include_graphics("imgs/KM_Surv_treatment.png")

## ----results = 'hide', message=FALSE, warning=FALSE, eval = F-----------------
#  NRAS_SNPs<-head(ListSNPs_NRAS,5) #["rs2157615","rs61731685","rs370560636","rs116293337"]
#  
#  surv<-MMRFVariant_SurvivalKM(patient,
#                                trt,
#                                variant.ann,
#                                NRAS_SNPs,
#                                FilterBy="Stage",
#                                filename=NULL,
#                                xlim = c(100,3000),
#                                conf.range = FALSE,
#                                color = c("Dark2"))
#  
#  
#  
#  

## ----figurename3, echo=FALSE, fig.cap="KM survival curves are drawn for each SNP found in the NRAS gene in the case of the FilterBy parameter is set to Stage", out.width = '50%'----
knitr::include_graphics("imgs/KM_Surv_stage.png")

## ----results = 'hide', message=FALSE, warning=FALSE, eval = F-----------------
#  NRAS_SNPs<-head(ListSNPs_NRAS,5) #["rs2157615","rs61731685","rs370560636","rs116293337"]
#  
#  surv<-MMRFVariant_SurvivalKM(patient,
#                                trt,
#                                variant.ann,
#                                NRAS_SNPs,
#                                FilterBy="Ethnicity",
#                                filename=NULL,
#                                xlim = c(100,3000),
#                                conf.range = FALSE,
#                                color = c("Dark2"))
#  
#  
#  
#  

## ----figurename4, echo=FALSE, fig.cap="KM survival curves are drawn for each SNP found in the NRAS gene in the case of the FilterBy parameter is set to Ethnicity", out.width = '99%'----
knitr::include_graphics("imgs/KM_Surv_ethnicity.png")

## ----results = 'hide', message=FALSE, warning=FALSE, eval = F-----------------
#  NRAS_SNPs<-head(ListSNPs_NRAS,5) #["rs2157615","rs61731685","rs370560636","rs116293337"]
#  
#  surv<-MMRFVariant_SurvivalKM(patient,
#                                trt,
#                                variant.ann,
#                                NRAS_SNPs,
#                                FilterBy="Bestresp",
#                                filename=NULL,
#                                xlim = c(100,3000),
#                                conf.range = FALSE,
#                                color = c("Dark2"))
#  
#  
#  
#  

## ----figurename5, echo=FALSE, fig.cap="KM survival curves are drawn for each SNP found in the NRAS gene in the case of the FilterBy parameter is set to Bestresp", out.width = '99%'----
knitr::include_graphics("imgs/KM_Surv_Bestresp.png")

## ----results = 'hide', message=FALSE, warning=FALSE, eval = F-----------------
#  
#  ListSNPs<-c("rs755588843","rs200556051","rs745587729","rs2066497","rs760494041")
#  
#  impact.table<-MMRFVariant_GetImpact(variant.ann,ListSNPs)

## ----figurename6, echo=FALSE, fig.cap="The table shows the MMRFVariant_GetImpact output. ", out.width = '99%'----
knitr::include_graphics("imgs/ImpactTable.png")





#' @title VariantFinder_SurvivalKM
#' @description Creates a KM survival plot from MMRF-RG patient clinical data
#' using survival library. It uses the fields D_PT_deathdy, D_PT_PRIMARYREASON and D_PT_lstalive
#' columns for groups.
#' @param patient is the data.frame of the patient clinical data downloaded from MMRF-Commpass Researcher Gateway 
#' (i.e. MMRF_CoMMpass_IA15_PER_PATIENT file) and imported into environment.
#' @param trt is the data.frame of the patient clinical data (i.e. treatment-response) downloaded from MMRF-Commpass Researcher Gateway 
#' (i.e. MMRF_CoMMpass_IA15_STAND_ALONE_TRTRESP file) and imported into environment.
#' @param variant.ann is the data.frame of the annotated variants downloaded from MMRF-Commpass Researcher Gateway 
#' (i.e. MMRF_CoMMpass_IA15a_All_Canonical_Variants file) and imported into environment.
#' @param Listvariant is the list of the variants to analyze.
#' @param FilterBy Column with groups to plot. This is a mandatory field.
#' Example:
#' \tabular{ll}{
#'Ethnicity \tab Ethnicity \cr
#'Stage \tab ISS Stage \cr
#'Treatment \tab  Treatment class \cr
#'Bestresp \tab Best overall response 	\cr
#'Gender \tab gender 	\cr
#'Effect \tab effect 	\cr
#'Biotype \tab biotype 	\cr
#'}
#' @param expand show or not an expanded plot
#' @param legend Legend title of the figure
#' @param xlim x axis limits e.g. xlim = c(0, 1000). Present narrower X axis, but not affect survival estimates.
#' @param main main title of the plot
#' @param labels labels of the plot
#' @param ylab y axis text of the plot
#' @param xlab x axis text of the plot
#' @param filename The name of the pdf file.
#' @param color Define the colors/Pallete for lines.
#' @param width Image width
#' @param height Image height
#' @param pvalue show p-value of log-rank test
#' @param conf.range  show confidence intervals for point estimates of survival curves.
#' @param dpi Figure quality
#' @import survminer
#' @import survival
#' @import gridExtra
#' @import ggplot2
#' @import stringr
#' @export
#' @return Survival plot
#' @examples
#' 
#' patient <- data.frame(public_id=c("MMRF_0000","MMRF_0001",
#'                                    "MMRF_0002","MMRF_0003",
#'                                    "MMRF_0004","MMRF_0005",
#'                                    "MMRF_0006","MMRF_0007",
#'                                    "MMRF_0008","MMRF_0009"),
#'                       D_PT_PRIMARYREASON = c("Death",NA,NA,"Death","Death", 
#'                                              "Death","NA","NA","Death","Death"),
#'                       D_PT_deathdy =  c(NA,NA,2226,172,NA,NA,1328,681,786,NA),
#'                       D_PT_lstalive = c(250,356,NA,NA,1814,223,NA,NA,NA,1450),
#'                       DEMOG_GENDER = c(rep(1,5),rep(2,5)), #2 stands for female, 1 standas for male#'                       
#'                       DEMOG_ETHNICITY=c(rep("Hispanic or Latino",5),rep("Not Hispanic or Latino",5)),
#'                       D_PT_issstage_char=c(rep("Stage III",3),rep("Stage II",2),rep("Stage I",5))
#'  )
#' 
#' 
#' trt<- data.frame(public_id=c("MMRF_0000","MMRF_0001",
#'                               "MMRF_0002","MMRF_0003",
#'                               "MMRF_0004","MMRF_0005",
#'                               "MMRF_0006","MMRF_0007",
#'                               "MMRF_0008","MMRF_0009"),
#'                  trtclass=c(rep("Bortezomib-based",2),rep("IMIDs-based",5),rep("combined bortezomib/IMIDs-based",3)),                                                    
#'                  bestresp=c(rep("Partial Responsed",5),rep("Very Good Partial Response",5))               
#'                                    
#'  )
#' 
#' 
#' 
#'          surv1<-VariantFinder_SurvivalKM(patient,  
#                                          trt,
#                                          variant.ann,
#                                          Listvariant1,
#                                          FilterBy="treatment", 
#                                          filename=NULL,
#                                          xlim = c(100,3000),
#                                          conf.range = FALSE,
#                                          color = c("Dark2"))

VariantFinder_SurvivalKM <- function(
  patient,
  trt,
  variant.ann,
  list.variant,
  FilterBy = "treatment",
  legend = "Legend",
  labels = NULL,
  xlim = NULL,
  main = "Kaplan-Meier Survival Curve",
  ylab = "Probability of survival",
  xlab = "Time since diagnosis (days)",
  filename = "survival.pdf",
  color = NULL,
  height = 8,
  width = 12,
  dpi = 300,
  pvalue = TRUE,
  conf.range = TRUE) {
  
  
  
  if(is.null(trt) || (is.null(patient))){
    stop("Please provide the file of treatment-response and/or patient.")
  }
  
  
  plot.list <- list()
  surv<-NULL
  
  
  
  for(i in 1:length(list.variant)){
    print(list.variant[i])
    
    result = tryCatch({ 
      surv<-VariantFinder_SurvivalKM_Single(patient,  
                                          trt,
                                          variant.ann,
                                          list.variant[i],
                                          FilterBy=FilterBy, 
                                          xlim=xlim,
                                          legend=legend,
                                          height=height,
                                          widt=width,
                                          conf.range = FALSE,
                                          color = c("Dark2"),
                                          pvalue=pvalue)
      
      plot.list[[i]]<-surv
      
    }, error = function(e) {
      print(paste("Only one group found with respect to: ",list.variant[i]))
      i<-i+1
    },warning = function(e) {
      print(paste("Only one group found with respect to:",list.variant[i]))
      i<-i+1
    },finally = function(f) {
      print(paste("Only one group found with respect to:\n",list.variant[i]))
      i<-i+1
    }
    
    )  #tryCatch
    
    
  }  #for
  
  
  plot.list<-plot.list[!sapply(plot.list,is.null)]
  
  
  
  if((length(list.variant) %% 2) == 0) {
    nr<-length(list.variant)/2
  } else {
    nr<-trunc(length(list.variant)/2)+1
    
  }
  
  plt<-arrange_ggsurvplots(plot.list, print = TRUE,title=main,
                           ncol = 2, nrow = nr)
  
  # plt<-arrange_ggsurvplots(plot.list, print = TRUE,title=main,
  #                           ncol = 2, nrow = length(list.variant)/2)
  
  
  width <-20
  height <-30
  
  if (!is.null(filename) & length(plt)>0) {
    
    filenm<-paste0(filename,".pdf")
    path<-file.path(getwd())
    path<-paste0(path,"/","ResultsPlot","/",filenm)
    
    ggsave(
      #plt$plot,
      plt,
      filename = path,
      #device = pdf,
      width = width,
      height = height,
      units = "in"
    )
    message(paste0("File saved as: ", path))
    
    
    
  } 
  
  return(plt)
  
}


#Plot variants by gene list

#' @title VariantFinder_PlotVariantsbyGene
#' @description
#' draws heatmap of the annotated variants occurrences 
#' @param topN is the top number of variant count
#' @param variant.ann is the dataframe of annotated variants downloaded from MMRF-Commpass Researcher Gateway 
#' (i.e. MMRF_CoMMpass_IA14a_All_Canonical_Variants file) and imported into environment
#' @param ListGene is the list of the genes to analyze (optional).
#' @param topN is the top number of genes whose to visualize the variants occurrence number
#' @param filenm is the name of the png file. If filenm is Null, the plot is draw but it is not saved.
#' @param width Image width
#' @param height Image height
#' @import ggplot2
#' @import dplyr 
#' @export
#' @examples
#' variant.ann<- data.frame(public_id=c("MMRF_0000","MMRF_0001",
#'                                      "MMRF_0002","MMRF_0003",
#'                                      "MMRF_0004","MMRF_0005",
#'                                      "MMRF_0006","MMRF_0007",
#'                                      "MMRF_0008",""),                  
#'                  dbSNP=c(rep("rs755588843",2),rep("rs569344016",5),rep("rs2066497",2),rep(".",1)),                                                    
#'                  Effect=c(rep("intragenic_variant",3),
#'                            rep("missense_variant",2),
#'                            rep("intron_variant",1),
#'                            rep("5_prime_UTR_variant",4)),
#'                   Gene=c(rep("PRDM16",3),
#'                            rep("AGO1",2),
#'                            rep("FPGT-TNNI3K",1),
#'                            rep("TNNI3K",4)), 
#'                  REF=c(rep("C",3),
#'                            rep("G",2),
#'                            rep("A",1),
#'                            rep("T",4)),                            
#'                  ALT=c(rep("GGCCT",3),
#'                            rep("G",2),
#'                            rep("T",1),
#'                            rep("A",4)),    

#'                   Biotype=c(rep("protein_coding",6),
#'                            rep("antisense",2),
#'                            rep("processed_pseudogene",2)),   
#'                            
#'                  Impact= c(rep("MODERATE",2),rep("MODIFIER",2),
#'                             rep("LOW",3),rep("HIGH",2),rep("MODIFIER",1)),
#'                  feature_type= c(rep("ENST00000388718",2),rep("ENST00000344616",2),
#'                             rep("ENST00000431492",3),rep("ENST00000390268",2),rep("ENST00000316407",1)),
#'                             
#'                  SIFT= c(rep("0.035,0.035,0.057,0.057,0.035,0.042,0.04,0.058",2),rep("0.002,0.002,0.001,0.002",2),
#'                             rep("0.614,0.614,0.781",6)),           
#'                             
#'                  Polyphen2=c(rep("0.021,0.986,0.884,0.977",2),rep("0.99",2),
#'                             rep("0.614,0.781",6))           
#'                                                               
#'                                  
#'  )
#' 
#' 
#' 
#' 
#' 
#' variants.plot<-VariantFinder_PlotVariantsbyGene(variant.ann,ListGene,topN=10,)
#' variants.plot<-VariantFinder_PlotVariantsbyGene(variant.ann,topN=10)
#' @export
#' @return heatmap of the annotated variants occurrences 


VariantFinder_PlotVariantsbyGene<- function(variant.ann,ListGene=NULL,topN=20,filenm="PlotVariantsbyGene_heatmap",height=10, width=10){
  
  
  
  if(is.null(variant.ann) ){
    stop("Please provide the file of the annotated variants.")
  }
  
  
  variant.ann.sub<- VariantFinder_getVariantAnn(variant.ann)
  
  
  variant.ann.sub<-dplyr::select(variant.ann.sub,public_id,Gene,dbSNP,
                                 Effect,Gene,Biotype,Impact)
  
  
  
  
  
  variant.ann.sub<-unique(variant.ann.sub)
  variant.ann.sub<-subset(variant.ann.sub, variant.ann.sub$dbSNP != "." & !is.nan(variant.ann.sub$dbSNP))
  
  
  
  variant.ann.sub<- variant.ann.sub[,2:3] 
  
  variant.summary<-group_by(variant.ann.sub,Gene,dbSNP)%>% summarize(n=n(),.groups = 'drop')
  names(variant.summary)[3]<-"count"
  
  variant.summary<-variant.summary[order(variant.summary$count, decreasing = TRUE),]
  
  
  
  
  
  if(is.null(ListGene)){
    
    
    variant.summary<- head(variant.summary,topN)
    msg<-paste0(" - top count=",topN)
    
  } else{
    
    
    variant.summary<- variant.summary[variant.summary$Gene %in% ListGene, ] 
    variant.summary<- head(variant.summary,topN)
    msg<-""
  }
  
  
  
  
  
  plt<-ggplot(data = variant.summary, aes(x=dbSNP, y=Gene, fill=count)) + 
    geom_tile(color = "white", size = 0.2)+ylab("Genes") +
    xlab("dbSNP IDs") +
    theme_bw() +
    theme(text = element_text(size=11),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          plot.title = element_text(size=11),
          axis.title=element_text(size=12,face="bold"),
          axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(fill = "Case Count") +
    ggtitle(paste0("Heatmap of the number of variants occurrence",msg))
  
  
  
  
  
  
  
  if (!is.null(filenm)) {
    
    filename<-paste0(filenm,".pdf")
    path<-file.path(getwd())
   # path<-paste0(path,"/","ResultsPlot","/")
    path<-paste0(path,"/","ResultsPlot","/",filename)
    

   # ggsave(plot=plt,filename= "../ResultsPlot/hist.png",height=6, width=8)
    ggsave(plot=plt,filename= "../ResultsPlot/hist.png",height=6, width=8)
    
    ggsave(
      filename = filename,
      path = path,
      width = width,
      height = height,
      units = "in"
    )
    message(paste0("File saved as: ", path))
    
    
  }
  return(plt)
  
}







#----------------------------------------






#Plot variants by Effect and Impact categories

#' @title VariantFinder_PlotbyEffectImpact 
#' @description
#' draws plot of annotated variants by Impact and Effect
#' @param topN is the top number of variant count
#' @param variant.ann is the dataframe of annotated variants downloaded from MMRF-Commpass Researcher Gateway 
#' (i.e. MMRF_CoMMpass_IA14a_All_Canonical_Variants file) and imported into environment
#' @param ListSNP is the list of the SNPs to analyze.
#' @param filenm is the name of the png file. If filenm is Null, the plot is draw but it is not saved.
#' @param width Image width
#' @param height Image height
#' @import ggplot2
#' @import dplyr 
#' @export
#' @examples
#' variant.ann<- data.frame(public_id=c("MMRF_0000","MMRF_0001",
#'                                      "MMRF_0002","MMRF_0003",
#'                                      "MMRF_0004","MMRF_0005",
#'                                      "MMRF_0006","MMRF_0007",
#'                                      "MMRF_0008",""),                  
#'                  dbSNP=c(rep("rs755588843",2),rep("rs569344016",5),rep("rs2066497",2),rep(".",1)),                                                    
#'                  Effect=c(rep("intragenic_variant",3),
#'                            rep("missense_variant",2),
#'                            rep("intron_variant",1),
#'                            rep("5_prime_UTR_variant",4)),
#'                   Gene=c(rep("PRDM16",3),
#'                            rep("AGO1",2),
#'                            rep("FPGT-TNNI3K",1),
#'                            rep("TNNI3K",4)), 
#'                  REF=c(rep("C",3),
#'                            rep("G",2),
#'                            rep("A",1),
#'                            rep("T",4)),                            
#'                  ALT=c(rep("GGCCT",3),
#'                            rep("G",2),
#'                            rep("T",1),
#'                            rep("A",4)),    

#'                   Biotype=c(rep("protein_coding",6),
#'                            rep("antisense",2),
#'                            rep("processed_pseudogene",2)),   
#'                            
#'                  Impact= c(rep("MODERATE",2),rep("MODIFIER",2),
#'                             rep("LOW",3),rep("HIGH",2),rep("MODIFIER",1)),
#'                  feature_type= c(rep("ENST00000388718",2),rep("ENST00000344616",2),
#'                             rep("ENST00000431492",3),rep("ENST00000390268",2),rep("ENST00000316407",1)),
#'                             
#'                  SIFT= c(rep("0.035,0.035,0.057,0.057,0.035,0.042,0.04,0.058",2),rep("0.002,0.002,0.001,0.002",2),
#'                             rep("0.614,0.614,0.781",6)),           
#'                             
#'                  Polyphen2=c(rep("0.021,0.986,0.884,0.977",2),rep("0.99",2),
#'                             rep("0.614,0.781",6))           
#'                                                               
#'                                  
#'  )
#' 
#' 
#' 
#' 
#' summary.plot<-VariantFinder_PlotbyEffectImpact(variant.ann,topN=50,height=10, width=15, filenm="PlotbyEffect")
#' @export
#' @return plot with the top count of the dbSNP variant categorized by Impact


VariantFinder_PlotbyEffectImpact<- function(variant.ann, ListSNPs, topN=20,filenm="PlotbyEffectImpact", height=10, width=10){
  
  if(is.null(variant.ann) || is.null(ListSNPs)){
    stop("Please provide the file of the annotated variants or the SNPs List.")
  }
  
  
  variant.ann.sub<- VariantFinder_getVariantAnn(variant.ann)
  
  variant.ann.sub<-variant.ann.sub[variant.ann.sub$dbSNP %in% ListSNPs,] 
  
  variant.ann.sub<-subset(variant.ann.sub, variant.ann.sub$SIFT != "." & 
                            variant.ann.sub$Polyphen2 != ".")
  
  
  
  
  variant.ann.sub<-dplyr::select(variant.ann.sub,public_id,dbSNP,
                                 Effect,Gene,Biotype,Impact)
  
  
  
  
  
  variant.ann.sub<-unique(variant.ann.sub)
  variant.ann.sub<-subset(variant.ann.sub, variant.ann.sub$dbSNP != "." & !is.nan(variant.ann.sub$dbSNP))
  
  
  
  #df.merge<-merge(x = variant.ann, y = trt, by = "public_id", type=left)
  
  
  
  
  
  
  variant.summary<-variant.ann.sub %>% group_by(dbSNP,Effect,Impact) %>% summarize(n(),.groups = 'drop')
  names(variant.summary)[4]<-"count"
  
  
  
  variant.summary<-variant.summary[order(variant.summary$count, decreasing = TRUE),]
  
  
  plt<-ggplot(head(variant.summary,topN), aes(x = count, y = dbSNP,shape=Effect,color=Impact)) + 
    geom_point(size=4) #+ facet_grid(~Effect) +
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=11),
        axis.title=element_text(size=12,face="bold"),
        axis.text.x = element_text(hjust = 1,face="bold",size = 11),
        axis.text.y = element_text(hjust = 1,face="bold",size = 11))+
    ggtitle("Variants by Effect and Impact")+ labs(x = "N variant")+labs(y = "dbSNP IDs")+
    
    
    
    
    if (!is.null(filenm)) {
      
      filename<-paste0(filenm,".pdf")
      path<-file.path(getwd())
      path<-paste0(path,"/","ResultsPlot","/")
      
      ggsave(
        filename = filename,
        path = path,
        width = width,
        height = height,
        units = "in"
      )
      message(paste0("File saved as: ", path))
      
    }
  return(plt)
  
}



#----------------------------------------
#----------------------------------------



#Get Impact by variants 

#' @title VariantFinder_GetImpact 
#' @description
#' get Impact of the variants in the input ListSNPs
#' @param variant.ann is the dataframe of annotated variants downloaded from MMRF-Commpass Researcher Gateway 
#' (i.e. MMRF_CoMMpass_IA14a_All_Canonical_Variants file) and imported into environment
#' @param ListSNPs is the list of the SNPs to analyze.
#' @import dplyr 
#' @export
#' @examples
#' variant.ann<- data.frame(public_id=c("MMRF_0000","MMRF_0001",
#'                                      "MMRF_0002","MMRF_0003",
#'                                      "MMRF_0004","MMRF_0005",
#'                                      "MMRF_0006","MMRF_0007",
#'                                      "MMRF_0008",""),                  
#'                  dbSNP=c(rep("rs755588843",2),rep("rs569344016",5),rep("rs2066497",2),rep(".",1)),                                                    
#'                  Effect=c(rep("intragenic_variant",3),
#'                            rep("missense_variant",2),
#'                            rep("intron_variant",1),
#'                            rep("5_prime_UTR_variant",4)),
#'                   Gene=c(rep("PRDM16",3),
#'                            rep("AGO1",2),
#'                            rep("FPGT-TNNI3K",1),
#'                            rep("TNNI3K",4)), 
#'                  REF=c(rep("C",3),
#'                            rep("G",2),
#'                            rep("A",1),
#'                            rep("T",4)),                            
#'                  ALT=c(rep("GGCCT",3),
#'                            rep("G",2),
#'                            rep("T",1),
#'                            rep("A",4)),    

#'                   Biotype=c(rep("protein_coding",6),
#'                            rep("antisense",2),
#'                            rep("processed_pseudogene",2)),   
#'                            
#'                  Impact= c(rep("MODERATE",2),rep("MODIFIER",2),
#'                             rep("LOW",3),rep("HIGH",2),rep("MODIFIER",1)),
#'                  feature_type= c(rep("ENST00000388718",2),rep("ENST00000344616",2),
#'                             rep("ENST00000431492",3),rep("ENST00000390268",2),rep("ENST00000316407",1)),
#'                             
#'                  SIFT= c(rep("0.035,0.035,0.057,0.057,0.035,0.042,0.04,0.058",2),rep("0.002,0.002,0.001,0.002",2),
#'                             rep("0.614,0.614,0.781",6)),           
#'                             
#'                  Polyphen2=c(rep("0.021,0.986,0.884,0.977",2),rep("0.99",2),
#'                             rep("0.614,0.781",6))           
#'                                                               
#'                                  
#'  )
#' 
#' 
#' ListSNPs<-c("rs755588843","rs200556051","rs745587729","rs2066497","rs760494041")
#' 
#' 
#' impact.table<-VariantFinder_GetImpact(variant.ann,ListSNPs)
#' @export
#' @return summary table aggregating and summarizing data about annotated variants

#-------------------------------------------------------------------------

VariantFinder_GetImpact<- function(variant.ann, ListSNPs){
  
  if(is.null(variant.ann) || is.null(ListSNPs)){
    stop("Please provide the file of the annotated variants/ListSNPs")
  }
  
  # rename columns
  variant.ann.temp<-VariantFinder_getVariantAnn(variant.ann)
  

  variant.ann.sub<-variant.ann.temp[variant.ann.temp$dbSNP %in% ListSNPs,] 
  
  
  
  
  
  variant.ann.sub<-subset(variant.ann.sub, variant.ann.sub$SIFT != "." & 
                            variant.ann.sub$Polyphen2 != ".")
  
  
  tryCatch(
    #try to do this
    {
      #add column (Max damaging impact)
      variant.ann.sub<-addImpacts(variant.ann.sub)
    },
    #if an error occurs, tell me the error
    error=function(e) {
      message('An Error Occurred: no SIFT and/or Polyphen2 scores are in dataframe with respect to the input Genes.')
      print(e)
    },
    #if a warning occurs, tell me the warning
    warning=function(w) {
      message('A Warning Occurred')
      print(w)
      return(NA)
    }
  )
  
  
  
  
  impact.table<-table_formatter(variant.ann.sub)
  impact.table<-subset(impact.table, impact.table$min_SIFT != "" & 
                         impact.table$max_polyphen != "")
  
  impact.table<-arrange(impact.table, min_SIFT, desc(max_polyphen))
  
  
  return(unique(impact.table))  
  
}





#Get variants---------------------------------------------------------------------------



#' @title VariantFinder_GetVariantsbyGene
#' @description
#' get variants occurring in a gene set provided as input
#' @param variant.ann is the dataframe of annotated variants downloaded from MMRF-Commpass Researcher Gateway 
#' (i.e. MMRF_CoMMpass_IA14a_All_Canonical_Variants file) and imported into environment
#' @param ListGene is the list of the genes to analyze. if it is Null, the list of variants oredere by the occurrence nuomber is shown
#' @export
#' @examples
#' variant.ann<- data.frame(public_id=c("MMRF_0000","MMRF_0001",
#'                                      "MMRF_0002","MMRF_0003",
#'                                      "MMRF_0004","MMRF_0005",
#'                                      "MMRF_0006","MMRF_0007",
#'                                      "MMRF_0008",""),                  
#'                  dbSNP=c(rep("rs755588843",2),rep("rs569344016",5),rep("rs2066497",2),rep(".",1)),                                                    
#'                  Effect=c(rep("intragenic_variant",3),
#'                            rep("missense_variant",2),
#'                            rep("intron_variant",1),
#'                            rep("5_prime_UTR_variant",4)),
#'                   Gene=c(rep("PRDM16",3),
#'                            rep("AGO1",2),
#'                            rep("FPGT-TNNI3K",1),
#'                            rep("TNNI3K",4)), 
#'                  REF=c(rep("C",3),
#'                            rep("G",2),
#'                            rep("A",1),
#'                            rep("T",4)),                            
#'                  ALT=c(rep("GGCCT",3),
#'                            rep("G",2),
#'                            rep("T",1),
#'                            rep("A",4)),    

#'                   Biotype=c(rep("protein_coding",6),
#'                            rep("antisense",2),
#'                            rep("processed_pseudogene",2)),   
#'                            
#'                  Impact= c(rep("MODERATE",2),rep("MODIFIER",2),
#'                             rep("LOW",3),rep("HIGH",2),rep("MODIFIER",1)),
#'                  feature_type= c(rep("ENST00000388718",2),rep("ENST00000344616",2),
#'                             rep("ENST00000431492",3),rep("ENST00000390268",2),rep("ENST00000316407",1)),
#'                             
#'                  SIFT= c(rep("0.035,0.035,0.057,0.057,0.035,0.042,0.04,0.058",2),rep("0.002,0.002,0.001,0.002",2),
#'                             rep("0.614,0.614,0.781",6)),           
#'                             
#'                  Polyphen2=c(rep("0.021,0.986,0.884,0.977",2),rep("0.99",2),
#'                             rep("0.614,0.781",6))           
#'                                                               
#'                                  
#'  )
#' 
#' 
#' 
#' 
#' 
#' variants.list<-VariantFinder_GetVariantsbyGene(variant.ann,ListGene)
#' 
#' variants.list<-VariantFinder_GetVariantsbyGene(variant.ann)
#' 
#' @export
#' @return the list of SNPs found in the gene list provided as input 


VariantFinder_GetVariantsbyGene<- function(variant.ann,ListGene=NULL){
  
  
  if(is.null(variant.ann) ){
    stop("Please provide the file of the annotated variants.")
  }
  
  
  variant.ann.sub<- VariantFinder_getVariantAnn(variant.ann)
  
  
  variant.ann.sub<-dplyr::select(variant.ann.sub,public_id,Gene,dbSNP,
                                 Effect,Gene,Biotype,Impact)
  
  
  
  
  
  variant.ann.sub<-unique(variant.ann.sub)
  variant.ann.sub<-subset(variant.ann.sub, variant.ann.sub$dbSNP != "." & !is.nan(variant.ann.sub$dbSNP))
  
  
  
  variant.ann.sub<- variant.ann.sub[,2:3] 
  
  variant.list<-group_by(variant.ann.sub,Gene,dbSNP)%>% summarize(n=n(),.groups = 'drop')
  names(variant.list)[3]<-"count"
  
  variant.list<-variant.list[order(variant.list$count, decreasing = TRUE),]
  
  
  
  if(!is.null(ListGene) ){
    
    variant.list<- variant.list[variant.list$Gene %in% ListGene, ] 
    
  }
  
  
  
  return(variant.list)
  
}














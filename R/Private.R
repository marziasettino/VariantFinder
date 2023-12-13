
#private
VariantFinder_getVariantAnn<- function(variant.ann){
  names(variant.ann)[1]<-"public_id"
  variant.ann$public_id<-substr(variant.ann$public_id,1,9)
  
  
  names(variant.ann)[names(variant.ann) == 'ANN....EFFECT'] <- "Effect"
  names(variant.ann)[names(variant.ann) == 'ANN....GENE'] <- "Gene"
  names(variant.ann)[names(variant.ann) == 'ANN....BIOTYPE'] <- "Biotype"
  names(variant.ann)[names(variant.ann) == 'ANN....IMPACT'] <- "Impact"
  names(variant.ann)[names(variant.ann) == 'ID'] <- "dbSNP"
  names(variant.ann)[names(variant.ann) == 'ANN....FEATURE'] <- "feature_type"
  names(variant.ann)[names(variant.ann) == 'ANN....FEATUREID'] <- "feature"
  names(variant.ann)[names(variant.ann) == 'dbNSFP_SIFT_score'] <- "SIFT"
  names(variant.ann)[names(variant.ann) == 'dbNSFP_Polyphen2_HDIV_score'] <- "Polyphen2"
  
  
  
  
  return(variant.ann)
}



#private


sift_impact <- function(x) {
  if (x>=0 && x<=0.05) {
    "Intolerant"
  } else if (x>=0.051 && x<=0.1) {
    "Potentially Intolerant" 
  } else if (x>=0.101 && x<=0.2) {
    "Borderline"
  } else if (x>=0.201 && x<=1) {
    "Tolerant"
  } else {
    stop("Invalid `x` value")
  }
}

polyphen_impact <- function(x) {
  if (x>=0 && x<=0.999) {
    "Benign"
  } else if (x>=1 && x<=1.24) {
    "Borderline" 
  } else if (x>=1.25 && x<=1.49) {
    "Potentially damaging"
  } else if (x>=1.5 && x<=1.99) {
    "Possibly damaging"
  } else if (x>=2) {
    "probably damaging"  
  } else {
    stop("Invalid `x` value")
  }
}




#private
table_formatter <- function(df) {
  
  df<-dplyr::select(df,dbSNP,feature,min_SIFT,max_polyphen,
                    SIFT_Impact,Polyphen_Impact,Effect,Gene,
                    REF,ALT,Biotype,Impact,feature_type)
  
  
  
  
  
  df<-formattable(df, list(
    SIFT_Impact = formatter("span", 
                            style = ~style(display = "block",
                                           font.weight = "bold", 
                                           color = "white",
                                           "border-radius" = "4px",
                                           "padding-right" = "4px",
                                           "background-color" =  
                                             ifelse(SIFT_Impact == "Intolerant","red",
                                                    ifelse(SIFT_Impact == "Potentially Intolerant","orange",
                                                           ifelse(SIFT_Impact == "Borderline","blue",    
                                                                  ifelse(SIFT_Impact=="Tolerant", "lightblue",NA)))))),
    
    Polyphen_Impact = formatter("span", 
                                style = ~style(display = "block",
                                               font.weight = "bold", 
                                               color = "white",
                                               "border-radius" = "4px",
                                               "padding-right" = "4px",
                                               "background-color" =  
                                                 ifelse(Polyphen_Impact == "probably damaging","red",
                                                        ifelse(Polyphen_Impact == "Possibly damaging","orange",
                                                               ifelse(Polyphen_Impact == "Potentially damaging","green",
                                                                      ifelse(Polyphen_Impact == "Borderline","blue",
                                                                             ifelse(Polyphen_Impact=="Benign", "lightblue",NA))))))),
    
    Gene = formatter("span", style = ~style(font.weight = "bold", color =  "blue")),
    dbSNP = formatter("span", style = ~style(font.weight = "bold", color =  "blue"))
    
    
  ))
  
  return(df) 
}



#private


addImpacts <- function(df) {
  
  #df<- VariantFinder_getVariantAnn(df)
  df<-dplyr::select(df,public_id,dbSNP,
                    Effect,Gene,REF,ALT,
                    Biotype,Impact,feature_type,
                    feature,SIFT,Polyphen2)
  
  #add column (Max damaging impact)
  df$min_SIFT<-""
  df$SIFT_Impact<-""
  
  df$max_polyphen<-""
  df$Polyphen_Impact<-""
  
  
  
  for(i in 1:nrow(df)){ 
    sift<-df$SIFT[i]
    polyphen<-df$Polyphen2[i]
    
    print(paste0("i: ", i))
    print(paste0("SIFT: ", sift))
    print(paste0("POLYPHEN: ", polyphen))
    
    
    result = tryCatch({
      sift<-unlist(strsplit(sift, ",")) %>% as.numeric(sift)
      min<-min(sift)
      df$min_SIFT[i]<-min
      df$SIFT_Impact[i]<-sift_impact(sift<-df$min_SIFT[i])
      
    }, warning = function(w) {
      df$min_SIFT[i]<-NA
    }, error = function(e) {
      df$min_SIFT[i]<-NA
    }
    ) 
    
    
    result = tryCatch({
      polyphen<-unlist(strsplit(polyphen, ",")) %>% as.numeric(polyphen)
      max<-max(polyphen)
      df$max_polyphen[i]<-max
      
      df$Polyphen_Impact[i]<-polyphen_impact(sift<-df$max_polyphen[i])
      
    }, warning = function(w) {
      df$max_polyphen[i]<-NA
    }, error = function(e) {
      df$max_polyphen[i]<-NA
    }
    )     
    
    
    
  }  #for  
  
  return(df) 
}

#------------------------------------------------

VariantFinder_SelectMerge<- function(variant.ann, patient, trt, ListSNPs){
  
  if(is.null(ListSNPs) || is.null(variant.ann) || is.null(patient)|| is.null(trt)){
    stop("Please provide the patient / variant treatment files.")
  }
  
  
  
  names(variant.ann)[1]<-"public_id"
  variant.ann<-dplyr::select(variant.ann,public_id,X.CHROM,POS,ID,REF,
                             ALT,ANN....ALLELE,
                             ANN....EFFECT,ANN....IMPACT,
                             ANN....GENE,ANN....FEATURE,ANN....BIOTYPE
  )
  
  
  
  variant.ann$public_id<-substr(variant.ann$public_id,1,9)
  variant.ann<-variant.ann[variant.ann$ID %in% ListSNPs,] 
  
  names(patient)[1]<-"public_id"
  patient<-dplyr::select(patient,public_id,D_PT_PRIMARYREASON,
                         D_PT_lstalive,D_PT_deathdy,
                         DEMOG_ETHNICITY,D_PT_issstage_char,DEMOG_GENDER
  )
  
  
  trt<-dplyr::select(trt,public_id,bestresp,trtclass)
  
  
  df.merge<-merge(x = patient, y = variant.ann, by = "public_id")
  df.merge<-merge(x = df.merge, y = trt, by = "public_id")
  
  
  
  df.merge <- df.merge %>%
    rename(Ethnicity=DEMOG_ETHNICITY,
           Stage=D_PT_issstage_char,
           Treatment=trtclass,
           Bestresp=bestresp,
           Gender=DEMOG_GENDER,
           dbSNP=ID,
           Effect=ANN....EFFECT,
           Biotype=ANN....BIOTYPE
    )
  
  return(df.merge)
}


#----------------------------------------------

#' @title VariantFinder_Survival_Single
#' @description Creates a survival plot from MMRF-RG patient clinical data
#' using survival library. It uses the fields D_PT_deathdy, D_PT_PRIMARYREASON and D_PT_lstalive
#' columns for groups.
#' @param patient is the data.frame of the patient clinical data downloaded from MMRF-Commpass Researcher Gateway 
#' (i.e. MMRF_CoMMpass_IA15_PER_PATIENT file) and imported into environment.
#' @param trt is the data.frame of the patient clinical data (i.e. treatment-response) downloaded from MMRF-Commpass Researcher Gateway 
#' (i.e. MMRF_CoMMpass_IA15_STAND_ALONE_TRTRESP file) and imported into environment.
#' @param variant.ann is the data.frame of the annotated variants downloaded from MMRF-Commpass Researcher Gateway 
#' (i.e. MMRF_CoMMpass_IA15a_All_Canonical_Variants file) and imported into environment.
#' @param Listvariant is the list of the variants to analyze.
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
#' @return Survival plot
#' @examples

VariantFinder_SurvivalKM_Single <- function(
  patient,
  trt,
  variant.ann,
  SNP,
  FilterBy = NULL,
  legend = "Legend",
  labels = NULL,
  xlim = NULL,
  main = "Kaplan-Meier Survival Curve",
  ylab = "Probability of survival",
  xlab = "Time since diagnosis (days)",
  # filename = "survival.pdf",
  color = NULL,
  height = 5,
  width = 12,
  dpi = 300,
  pvalue = TRUE,
  conf.range = TRUE) {
  
  
  
  df.merge<-VariantFinder_SelectMerge(variant.ann, patient, trt, SNP)
 # df.merge<-subset(df.merge, df.merge$Stage != "") 
  df.merge$Stage[df.merge$Stage ==""] <- "undefined" 
  
  
  if (is.null(color)) {
    color <- rainbow(length(unique(df.merge[, FilterBy])))
  }
  
  group <- NULL
  
  
  
  
  
  if (is.null(FilterBy)) {
    stop("Please provide a value for FilterBy parameter")
  } else {
    
    
    if (length(unique(df.merge[, FilterBy])) == 1) {
      stop(
        paste0(
          "Only one group is found with respect to",SNP)
      )
    }
    
  }#end
  
  
  
  notDead <- is.na(df.merge$D_PT_deathdy)
  
  if (any(notDead == TRUE)) {
    df.merge[notDead, "D_PT_deathdy"] <- df.merge[notDead, "D_PT_lstalive"]
  }
  
  #TRUE(DEAD)/FALSE (ALIVE)
  df.merge$s <- grepl("Death", df.merge$D_PT_PRIMARYREASON, ignore.case = TRUE)
  
  
  
  df.merge$type <- as.factor(df.merge[, FilterBy])
  df.merge <-  df.merge[, c("D_PT_deathdy", "s", "type")]
  #   formula 
  f.m <-formula(survival::Surv(as.numeric(df.merge$D_PT_deathdy), event = df.merge$s) ~ df.merge$type)
  fit <- do.call(survival::survfit, list(formula = f.m, data = df.merge))
  
  label.add.n <- function(x) {
    na.idx <- is.na(df.merge[, "D_PT_deathdy"])
    negative.idx <- df.merge[, "D_PT_deathdy"] < 0
    idx <- !(na.idx | negative.idx)
    return(paste0(x, " (n = ",
                  sum(df.merge[idx, "type"] == x), ")"))
  }
  
  if (is.null(labels)) {
    d <- survminer::surv_summary(fit, data = df.merge)
    order <-
      unname(sapply(levels(d$strata), function(x)
        unlist(str_split(x, "="))[2]))
    labels <- sapply(order, label.add.n)
  }
  
  
  
  
  if (length(xlim) == 1) {
    xlim <- c(0, xlim)
  }
  
  
  suppressWarnings({
    
    
    surv <- survminer::ggsurvplot( 
      fit=fit,
      data=df.merge,
      pval = pvalue,
      conf.range = conf.range,
      xlim = xlim,
      main = main,
      xlab = xlab,
      legend.title = legend,
      palette =  color,
      legend = "right",
      legend.labs = levels(FilterBy)
      
    )
  })
  
  
  SNP<-toString(SNP)
  
  
  surv<-surv+ labs(title = paste("dbSNP Variant:",SNP))
  
  
  
  
  return(surv)
  
}







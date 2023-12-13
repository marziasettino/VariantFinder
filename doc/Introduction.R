## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  warning=FALSE, 
  message=FALSE,
  comment = "#>"
)

## ----figurename2, echo=FALSE, fig.cap="Use-case diagram that represents the high-level functionalities of MMRFVariant", out.width = '99%'----
knitr::include_graphics("imgs/UseCase.png")

## ----setup--------------------------------------------------------------------
library(MMRFVariant)
library(dplyr)
library(DT)
library(ggplot2)
library(stringr)
library(ggpubr)
library(survminer)
library(survival)
library(formattable)

## ---- echo = FALSE,hide=TRUE, message=FALSE,warning=FALSE---------------------
devtools::load_all(".")

## ----message=FALSE, warning=FALSE, include=FALSE------------------------------
library(dplyr)
library(DT)
library(ggplot2)
library(stringr)
library(ggpubr)

## ----figurename8, echo=FALSE, fig.cap="MMRFVariant Functions", out.width = '99%'----
knitr::include_graphics("imgs/Functions.png")

## ----results = 'hide', message=FALSE, warning=FALSE, eval = F-----------------
#  patient <- read.csv("~/MMRF_CoMMpass_IA15_PER_PATIENT")
#  
#  trt <- read.csv("~/MMRF_CoMMpass_IA15_STAND_ALONE_TRTRESP")
#  
#  variant.ann <- read.csv("~/MMRF_CoMMpass_IA15a_All_Canonical_Variants")

## ----echo=TRUE, message=FALSE, warning=FALSE----------------------------------
datatable(variant.ann.example,options = list(scrollX = TRUE, keys = TRUE))

## ----echo=TRUE, message=FALSE, warning=FALSE----------------------------------
datatable(patient.example,options = list(scrollX = TRUE, keys = TRUE))

## ----echo=TRUE, message=FALSE, warning=FALSE----------------------------------
datatable(trt.example,options = list(scrollX = TRUE, keys = TRUE))

## ----sessionInfo--------------------------------------------------------------
sessionInfo()


## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(BioRssay)

## ----eval=FALSE---------------------------------------------------------------
#  #1. CRAN version
#  install.packages("BioRssay")
#  #2. Developmental version
#  if (!requireNamespace("devtools", quietly = TRUE))
#          install.packages("devtools")
#      devtools::install_github("milesilab/BioRssay", build_vignettes = TRUE)

## -----------------------------------------------------------------------------
data(bioassay)
head(bioassay$assay2)

## -----------------------------------------------------------------------------
file <- paste0(path.package("BioRssay"), "/Test.BioRssay.txt")
test<-read.table(file,header=TRUE)
head(test)

## -----------------------------------------------------------------------------
assays<-bioassay
exm1<-assays$assay2
head(exm1)
unique(as.character(exm1$strain))

## -----------------------------------------------------------------------------
dataT<-probit.trans(exm1) #additionally an acceptable threshold for controls' mortality can be set as desired with "conf="; default is 0.05.
dataT$convrg
head(dataT$tr.data)


## -----------------------------------------------------------------------------
data<-dataT$tr.data #probid transformed data
RR<-resist.ratio(data)
RR

## -----------------------------------------------------------------------------
model.signif(dataT$tr.data)

## ----echo=FALSE---------------------------------------------------------------
oldpar<-par(no.readonly = TRUE)

## ----fig.dim=c(8,4)-----------------------------------------------------------
strains<-levels(data$strain)
par(mfrow=c(1,2)) # set plot rows
# plot without confidence intervals and test of validity of the model
mort.plot(data,plot.conf=FALSE,test.validity=FALSE) 
# plot only the regression lines
mort.plot(data,plot.conf=FALSE,test.validity=FALSE,pch=NA) 
# same plots with confidence level
par(mfrow=c(1,2))
mort.plot(data,plot.conf=TRUE,test.validity=FALSE)
mort.plot(data,plot.conf=TRUE,test.validity=FALSE,pch=NA)

## ----echo=FALSE---------------------------------------------------------------
par(oldpar)

## -----------------------------------------------------------------------------
head(test)
unique(test$insecticide)
bend<-test[test$insecticide=="bendiocarb",]
head(bend)

## ----fig.dim=c(6,4)-----------------------------------------------------------
dataT.b<-probit.trans(bend)
data.b<-dataT.b$tr.data

RR.b<-resist.ratio(data.b,plot = T,ref.strain = "Kisumu",plot.conf = T, test.validity = T)
head(RR.b)

## -----------------------------------------------------------------------------
#To then test the difference in dose-mortality response between the strains
t.models<-model.signif(data.b)
t.models

## ----fig.dim=c(6,4)-----------------------------------------------------------
file <- paste0(path.package("BioRssay"), "/Example3.txt") #import the example file from the package
exm3<-read.table(file,header=TRUE)
trnd<-probit.trans(exm3) #probit transformation and correction of data
resist.ratio(trnd$tr.data,LD.value = c(50,95),plot = T) #get LD and RR values with the mortality plot
model.signif(trnd$tr.data) # test the models significance for each strain


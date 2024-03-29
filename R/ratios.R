#' Calculate confidence range for regressions
#' @noRd
CIplot<-function(mods,pmort_min,pmort_max,conf.level){
  summ<-summary(mods)
  zz<-qnorm((1-conf.level)/2,lower.tail=FALSE)
  a<-summ$coefficient[2] # mods slope
  b<-summ$coefficient[1] # mods intercept
  minldose<-((pmort_min-b)/a) # predicted dose for 0% mortality
  maxldose<-((pmort_max-b)/a) # predicted dose for 100% mortality
  datalfit<-seq(minldose-0.2,maxldose+0.2,0.01) # generates a set of doses
  datafit<-data.frame(dose=10^datalfit)
  pred<-predict.glm(mods,newdata=datafit,type="response",se.fit=TRUE)
  #generates the predicted mortality for the set of doses and SE
  ci<-cbind(pred$fit-zz*pred$se.fit,pred$fit+zz*pred$se.fit) #CI for each dose
  probitci<-cbind(datafit$dose,suppressWarnings(qnorm(ci)))
  # apply probit transformation to CI
  return(probitci)
}

#' Lethal dose glm test
#' @noRd
LD <- function(mod, conf.level,LD.value=c(25,50,95)) {
  p <- LD.value # leathal dose
  het = deviance(mod)/df.residual(mod)
  if(het < 1){het = 1} # Heterogeneity cannot be less than 1
  m.stats <- summary(mod, dispersion=het, cor = F)
  b0<-m.stats$coefficients[1] # Intercept (alpha)
  b1<-m.stats$coefficients[2] # Slope (beta)
  interceptSE <- m.stats$coefficients[3]
  slopeSE <- m.stats$coefficients[4]
  z.value <- m.stats$coefficients[6]
  vcov = summary(mod)$cov.unscaled
  var.b0<-vcov[1,1] # Intercept variance
  var.b1<-vcov[2,2] # Slope variance
  cov.b0.b1<-vcov[1,2] # Slope intercept covariance

  alpha<-1-conf.level
  if(het > 1) {
    talpha <- -qt(alpha/2, df=df.residual(mod))
  } else {
    talpha <- -qnorm(alpha/2)
  }
  g <- het * ((talpha^2 * var.b1)/b1^2)
  eta = family(mod)$linkfun(p/100)  #probit distribution curve
  theta.hat <- (eta - b0)/b1
  const1 <- (g/(1-g))*(theta.hat - cov.b0.b1/var.b1)
  const2a <- var.b0 - 2*cov.b0.b1*theta.hat + var.b1*theta.hat^2 - g*(var.b0 - (cov.b0.b1^2/var.b1))
  const2 <- talpha/((1-g)*b1) * sqrt(het * (const2a))
  #Calculate the confidence intervals LCL=lower,
  #UCL=upper (Finney, 1971, p. 78-79. eq. 4.35)
  LCL <- (theta.hat + const1 - const2)
  UCL <- (theta.hat + const1 + const2)
  #Calculate variance for theta.hat (Robertson et al., 2007, pg. 27)
  var.theta.hat <- (1/(theta.hat^2)) * ( var.b0 + 2*cov.b0.b1*theta.hat + var.b1*theta.hat^2 )
  ECtable <- c(10^c(rbind(theta.hat,LCL,UCL,var.theta.hat)),
               b1,slopeSE,b0,interceptSE,het,g)
  return(ECtable)
}

#' Test the significance of model pairs of strains
#' @noRd
reg.pair0<-function(data){
mortality<-cbind(data$dead,data$total-data$dead)
mod1<-glm(mortality~log10(data$dose)*data$strain,
          family = quasibinomial(link=probit))
mod2<-glm(mortality~log10(data$dose),
          family = quasibinomial(link=probit))
return(anova(mod2,mod1,test="Chi"))
}

reg.pair<-function(data){
  mortality<-cbind(data$dead,data$total-data$dead)
  mod1<-glm(mortality~log10(data$dose)*data$strain,
            family = quasibinomial(link=probit))
  mod2<-glm(mortality~log10(data$dose),
            family = quasibinomial(link=probit))
  Test<-anova(mod2,mod1,test="Chi")
  an<-as.data.frame(anova(mod1,test="F"))
  return(list(pairT=Test,fullM=an))
}


#' Get LD and RR values for each strain
#' @noRd
get.dxt<-function(strains,data,conf.level,LD.value){
  dxt<-lapply(strains,function(ss,data,conf.level,LD.value){
    tmp<-data[data$strain == ss,]
    y<-with(tmp,cbind(dead,total-dead))
    mods<-glm(y~log10(dose),data=tmp,family = quasibinomial(link=probit))
    dat <- LD(mods, conf.level,LD.value=LD.value)
    E<-mods$fitted.values*tmp$total # expected dead
    chq<-sum(((E-tmp$dead)^2)/(ifelse(E<1,1,E))) #if E is lower than 1 chi-sq
    #fails to detect the significance, ~change denominator to 1
    df<-length(tmp$dead)-1
    dat<-c(dat,chq,df,pchisq(q=chq,df=df,lower.tail=FALSE))
    return(list(mods,dat))
  },data=data,conf.level=conf.level,LD.value=LD.value)
  return(dxt)
}


###########

#' Calculate lethal dosage, resistance ratios, and regression coefficients
#'  and tests for linearity
#'
#' Using a generalized linear model (GLM, logit link function), this function
#' computes the lethal doses for 25%, 50% and 95% (unless otherwise provided)
#' of the population (LD25, LD50 and LD95, resp.), and their confidence
#' intervals (LDmax and LDmin, 0.95 by default). See details for more info.
#'
#' @param data a data frame of probit-transformed mortality data using the
#' function probit.trans()
#' @param conf.level numerical. level for confidence intervals to be applied
#' to the models (default 0.95)
#' @param LD.value numerical. Level of lethal dose to be tested.
#'  default=c(25,50,95)
#' @param ref.strain character. name of the reference strain if present
#'  (see details)
#' @param plot logical. Whether to draw the plot. Default FALSE
#' @param plot.conf logical. If plot=TRUE, whether to plot the 95 percent
#' confidence intervals. Default TRUE
#' @param test.validity logical. If plot=TRUE (default), the regression for a
#' strain that failed the linearity test is not plotted
#' @param legend.par arguments to be passed on to \code{legend()} as
#' in \code{mort.plot()}
#' @param ... parameters to be passed on to graphics for the plot
#' (e.g. col, pch)
#'
#' @importFrom graphics abline axis legend lines mtext
#' @importFrom stats deviance df.residual family pchisq predict.glm qt qnorm
#'
#' @details If a name is provided in ref.strain=, it will be used as the
#' reference to compute the resistance ratios (RR). Alternatively, the
#' function will look for a strain with the suffix "-ref" in the dataset.
#' If this returns NULL, the strain with the lowest LD50 will be considered as reference.
#'
#' In addition to LD values, the function in a nutshell uses a script modified
#' from Johnson et al (2013), which allows taking the g factor into account
#' ("With almost all good sets of data, g will be substantially smaller than
#' 1.0 and seldom greater than 0.4." Finney, 1971) and the heterogeneity (h)
#' of the data (Finney, 1971) to calculate the confidence intervals (i.e. a
#' larger heterogeneity will increase the confidence intervals). It also
#' computes the corresponding resistance ratios (RR), i.e. the ratios between
#' a given strain and the strain with the lower LD50 and LD95, respectively for
#'  RR50 and RR95 (usually, it is the susceptible reference strain), with their
#'  95% confidence intervals (RRmin and RRmax), calculated according to
#'  Robertson and Preisler (1992). Finally, it also computes the coefficients
#'   (slope and intercept, with their standard error) of the linear
#'   regressions) and tests for the linearity of the dose-mortality response
#'   using a chi-square test (Chi(p)) between the observed dead numbers (data)
#'   and the dead numbers predicted by the regression (the test is significant
#'   if the data is not linear, e.g. mixed populations).
#'
#' @return Returns a data frame with the various estimates mentioned above.
#'  If plot=TRUE, plots the mortality on a probit-transformed scale against
#'  the log_10 doses.
#'
#' @author Pascal Milesi, Piyal Karunarathne, Pierrick Labbé
#'
#' @references Finney DJ (1971). Probitanalysis. Cambridge:Cambridge
#' University Press. 350p.
#'
#' Hommel G (1988). A stage wise rejective multiple test procedure based on
#' a modified Bonferroni test. Biometrika 75, 383-6.
#'
#' Johnson RM, Dahlgren L, Siegfried BD, Ellis MD (2013). Acaricide,fungicide
#' and drug interactions in honeybees (Apis mellifera). PLoSONE8(1): e54092.
#'
#' Robertson, J. L., and H.K. Preisler.1992. Pesticide bioassays with
#' arthropods. CRC, Boca Raton, FL.
#'
#' @examples
#' data(bioassay)
#' transd<-probit.trans(bioassay$assay2)
#' data<-transd$tr.data
#' resist.ratio(data,plot=TRUE)
#'
#' @export
resist.ratio<-function(data,conf.level=0.95,LD.value=c(25,50,95),
                        ref.strain=NULL,plot=FALSE,plot.conf=TRUE,
                        test.validity=TRUE,legend.par=c("bottomright"),...) {
  if(!any(LD.value==50)){LD.value<-sort(c(LD.value,50))}

  data$strain<-factor(data$strain)
  strains<-levels(data$strain)
  dxt<-get.dxt(strains,data,conf.level,LD.value=LD.value)
  dat<-do.call(rbind,lapply(dxt,function(x){x[[2]]}))
  colnames(dat)<-c(paste0(paste0("LD",rep(LD.value,each=4)),
                          rep(c("","min","max","var"),2)),"Slope", "SlopeSE",
                   "Intercept", "InterceptSE", "h", "g","Chi2","df","Chi(p)")
  rownames(dat)<-strains
  if(is.null(ref.strain)){
    ref <- which(strains == strains[grep("-ref$",as.character(strains))],
                 arr.ind=TRUE)
  } else {
    ref=ref.strain
  }
  if (length(ref)==0) {
    refrow <- which(dat[,"LD50"]==min(dat[,"LD50"]),arr.ind=TRUE)
  } else {
    refrow <-ref
  }

  for(l in seq_along(LD.value)){
    assign(paste0("rr",LD.value[l]),
           dat[,paste0("LD",LD.value[l])]/dat[refrow,paste0("LD",LD.value[l])])
    assign(paste0("CI",LD.value[l]),
           1.96*sqrt(log10(dat[,paste0("LD",LD.value[l],"var")])+log10(dat[refrow,paste0("LD",LD.value[l],"var")])))
    assign(paste0("rr",LD.value[l],"max"),
           10^(log10(get(paste0("rr",LD.value[l])))+get(paste0("CI",LD.value[l]))))
    assign(paste0("rr",LD.value[l],"min"),
           10^(log10(get(paste0("rr",LD.value[l])))-get(paste0("CI",LD.value[l]))))
    ggl<-get(paste0("rr",LD.value[l],"max"))
    ggl[refrow]<-0
    ggl2<-get(paste0("rr",LD.value[l],"min"))
    ggl2[refrow]<-0
  }
  RR<-mget(c(paste0("rr",rep(LD.value,each=3),c("","min","max"))))
  RR<-do.call(cbind,RR)

  if(plot){
    mort.plot(data,strains=NULL,plot.conf,test.validity=test.validity,
              conf.level=conf.level,legend.par=legend.par,...)
  }
  nm<-colnames(dat)
  dat<-rbind.data.frame(dat[,-(grep("var",colnames(dat)))])
  nm<-nm[-c(grep("var",nm))]
  nm<-nm[c((ncol(dat)-8):ncol(dat),1:(ncol(dat)-9))]
  dat<-as.matrix(cbind(dat[,(ncol(dat)-8):ncol(dat)],dat[,1:(ncol(dat)-9)],RR))
  colnames(dat)<-c(nm,colnames(RR))
  dat<-ifelse(dat>10,round(dat,0),ifelse(dat>1,round(dat,2),round(dat,4)))
  return(dat)
}


#' Test the significance of dose-mortality response differences
#'
#' This function is used when comparing at least two strains. It tests whether
#'  the mortality-dose regressions are similar for different strains, using a
#'  likelihood ratio test (LRT). If there are more than two strains, it also
#'  computes pairwise tests, using sequential Bonferroni correction
#'  (Hommel, 1988) to account for multiple testing.
#'
#' @param data a data frame of probit transformed mortality data using the
#' function probit.trans
#' @importFrom stats anova na.omit
#' @importFrom utils combn
#'
#' @return a list with model outputs: a chi-square test if there are only two
#' strains or if there are more than two strains, first an overall model
#' assessment (i.e. one strain vs. all) and given overall model is significant,
#'  then a bonferroni test of significance from a pairwise model comparison.
#'
#' @details A global LRT test assesses a strain’s effect, by comparing two
#' models, one with and one without this effect (i.e. comparing a model with
#' several strains to a model where all the data originate from a single
#' strain).
#' If there are more than two strains, pairwise tests are computed, and
#' p-values of significance are assessed using sequential Bonferroni correction
#'  (Hommel, 1988) to account for multiple testing.
#'
#' Warning: We strongly encourage users to not use this function when the
#' dose-mortality response for at least one strain significantly deviates
#' from linearity (see resist.ratio() function for more details): in such
#' cases the test cannot be interpreted.
#'
#' @author Pascal Milesi, Piyal Karunarathne, Pierrick Labbé
#'
#' @examples
#' data(bioassay)
#' transd<-probit.trans(bioassay$assay2)
#' data<-transd$tr.data
#' model.signif(data)
#'
#' @export
model.signif<-function(data){
  data$strain<-as.factor(data$strain)
  strains<-levels(data$strain)
  if (length(strains)>=2) {
    Test<-reg.pair(data)
    if(length(strains)>2 & Test$pairT$`Pr(>Chi)`[2]>0.05){
      message("effect on strains are non-significant \n all strains come from the same population")
    } else if(length(strains)>2 & Test$pairT$`Pr(>Chi)`[2]<=0.05){
      print(Test$pairT)
      message("complete model is significant against a NULL model \n continueing to pair-wise comparison")
      Test<-sapply(strains, function(x,data) sapply(strains, function(y,data){
        if(x!=y){
          dat<-data[data$strain==x | data$strain==y,]
          reg.pair(dat)$pairT$Pr[2]
        }
      },data=data),data=data)
      dv<-sapply(strains, function(x,data){
        sapply(strains, function(y,data){
          if(x!=y){
            dat<-data[data$strain==x | data$strain==y,]
            reg.pair(dat)$pairT$Deviance[2]
          }
        },data=data)
      },data=data)

      rdl<-combn(strains,2,function(x,data){
        dat<-data[data$strain==x[1] | data$strain==x[2],]
        tmp<-reg.pair(dat)$fullM
        list(c(round(tmp$`Resid. Dev`[1],2),round(as.numeric(tmp$`Resid. Dev`[3:4]),3),round(as.numeric(tmp$`Pr(>F)`[3:4]),3)))
      },data=data)
      rdl<-do.call(rbind,rdl)

      rk<-(length(strains)*(length(strains)-1)/2):1
      rk<-0.05/rk
      toget<-t(combn(rownames(Test),2))
      pval<-unlist(Test[toget])

      pval<-cbind(1:length(pval),pval)
      pval<-pval[order(pval[,2],decreasing = T),]
      pval<-cbind(pval,rk)
      pval<-pval[order(pval[,1]),]


      Test<-data.frame(cbind(toget,round(unlist(Test[toget]),5),
                             ifelse(pval[,2]<pval[,3],"sig","non-sig")))
      rdl[which(Test[,3]>0.05),]<-NA

      ### bonferr for pvals
      bp<-cbind(1:(nrow(Test)*2),as.numeric(unlist(rdl[,4:5])))
      bp<-bp[order(unlist(bp[,2])),]
      bp0<-length(na.omit(unlist(bp[,2])))
      th<-unlist(lapply(1:bp0,function(x){0.05/(bp0-x+1)}))
      tmp<-suppressWarnings(cbind(bp,th))
      tmp[which(is.na(tmp[,2])),3]<-NA
      tmp<-round(tmp[order(tmp[,1]),3],5)
      tmp<-split(tmp,cut(seq_along(tmp),2,labels = F))

      Test<-cbind(Test,rdl[,1:3],rdl[,4],tmp$`1`,rdl[,5],tmp$`2`)
      colnames(Test)<-c("strain1","strain2","model.pval","bonferroni","res.Dv.Null","res.Dv.str","res.Dv.int","str.pval","str.thr","int.pval","int.thr")
    }
  } else {
    message("Only one strain present; check your data")
  }
  cat("Output details
  model.pval - significance value of ANOVA on the binomial GLM test of the strain pair
  bonferroni - significance of the model.pval with bonferroni correction
  res.Dv - residual deviance
  thr - threshold for the significance of the pvalue
  str - values for the strains
  int - values for the interaction between the strain and the dose
          ")
  return(list(model=Test))
}


install.packages("forestplot")

#hazard ratio

library("survival")

#reading the file post normalisation and missing value impute

omicsData <- read.csv("Cox_dataset.csv", header = TRUE)
head(omicsData) # check if has bruneck code and protein names

omicsData <- omicsData[,-24]

############## Univariate Cox PH ##############
covariates<-colnames(omicsData[,2:31])
# exclude time n cvd from min-max scaling
newData <- as.data.frame(scale(omicsData[,2:29]))

# https://stackoverflow.com/questions/8169323/r-concatenate-two-dataframes
fastmerge <- function(d1, d2) {
  d1.names <- names(d1)
  d2.names <- names(d2)
  
  # columns in d1 but not in d2
  d2.add <- setdiff(d1.names, d2.names)
  
  # columns in d2 but not in d1
  d1.add <- setdiff(d2.names, d1.names)
  
  # add blank columns to d2
  if(length(d2.add) > 0) {
    for(i in 1:length(d2.add)) {
      d2[d2.add[i]] <- NA
    }
  }
  
  # add blank columns to d1
  if(length(d1.add) > 0) {
    for(i in 1:length(d1.add)) {
      d1[d1.add[i]] <- NA
    }
  }
  print(head(d1))
  return(rbind(d1, d2))
}
d2 <- omicsData[,30:31] # extract time n CVD
#newOmics <- fastmerge(newData,d2)
newOmics <- cbind(newData,d2) #bind columns from scaled data and d2
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(CVD0010W,CVD0010)~',x)))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = newOmics)}) # attributes(univ_models)

univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         
                         
                         #univ
                         beta<- x$coef[1]
                         HR <- x$coef[2]
                         HR.confint.lower <- x$conf.int[3]
                         HR.confint.upper <- x$conf.int[4]
                         wald.test<- x$coef[4]
                         p.value<- x$coef[5]
                         
                         #two covariates
                         #beta<- x$coef[3]
                         #HR <- x$coef[6]
                         #wald.test<-x$coef[12]
                         #HR.confint.lower <- x$conf.int[9]
                         #HR.confint.upper <- x$conf.int[12]
                         #p.value<- x$coef[15]
                         
                         #six covariates
                         #beta<- x$coef[7]
                         #HR <- x$coef[14]
                         #wald.test<-x$coef[28]
                         #HR.confint.lower <- x$conf.int[21]
                         #HR.confint.upper <- x$conf.int[28]
                         #p.value<- x$coef[35]
                         
                         
                         
                         
                         #HR.confint.lower <- x$conf.int["lower .95"]
                         #HR.confint.upper <- x$conf.int["upper .95"]
                         #HR <- paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR,HR.confint.lower,HR.confint.upper, wald.test, p.value)
                         names(res)<-c("beta", "HR", "95% Lower CI for HR", "95% Upper CI for HR", "Z","p.value")
                         return(res)
                       })

res <- t(as.data.frame(univ_results, check.names = FALSE))
save.df <-as.data.frame(res)
write.csv(save.df, file = "Cox_HRFeatures.csv")

univ_formulas_AS <- sapply(covariates,
                        function(x) as.formula(paste('Surv(CVD0010W,CVD0010)~AGE+SEX+',x)))

univ_models_AS <- lapply( univ_formulas_AS, function(x){coxph(x, data = newOmics)}) # attributes(univ_models)

univ_results_AS <- lapply(univ_models_AS,
                       function(x){ 
                         x <- summary(x)
                         
                         
                         #two covariates
                         beta<- x$coef[3]
                         HR <- x$coef[6]
                         wald.test<-x$coef[12]
                         HR.confint.lower <- x$conf.int[9]
                         HR.confint.upper <- x$conf.int[12]
                         p.value<- x$coef[15]
                         
                         
                         #univ
                         #beta<- x$coef[1]
                         #HR <- x$coef[2]
                         #HR.confint.lower <- x$conf.int[3]
                         #HR.confint.upper <- x$conf.int[4]
                         #wald.test<- x$coef[4]
                         #p.value<- x$coef[5]
                         
                         
                         #six covariates
                         #beta<- x$coef[7]
                         #HR <- x$coef[14]
                         #wald.test<-x$coef[28]
                         #HR.confint.lower <- x$conf.int[21]
                         #HR.confint.upper <- x$conf.int[28]
                         #p.value<- x$coef[35]
                         
                         
                         
                         
                         #HR.confint.lower <- x$conf.int["lower .95"]
                         #HR.confint.upper <- x$conf.int["upper .95"]
                         #HR <- paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR,HR.confint.lower,HR.confint.upper, wald.test, p.value)
                         names(res)<-c("beta", "HR", "95% Lower CI for HR", "95% Upper CI for HR", "Z","p.value")
                         return(res)
                       })

res_AS <- t(as.data.frame(univ_results_AS, check.names = FALSE))
save.df_AS <-as.data.frame(res_AS)
write.csv(save.df_AS, file = "Cox_HRFeatures_AgeSex.csv")

univ_formulas_ASL <- sapply(covariates,
                        function(x) as.formula(paste('Surv(CVD0010W,CVD0010)~AGE+SEX+Triglycerides+Total_cholesterol+HDL_C+LDL_C+',x)))

univ_models_ASL <- lapply( univ_formulas_ASL, function(x){coxph(x, data = newOmics)}) # attributes(univ_models)

univ_results_ASL <- lapply(univ_models_ASL,
                       function(x){ 
                         x <- summary(x)
                         
                         
                         #univ
                         #beta<- x$coef[1]
                         #HR <- x$coef[2]
                         #HR.confint.lower <- x$conf.int[3]
                         #HR.confint.upper <- x$conf.int[4]
                         #wald.test<- x$coef[4]
                         #p.value<- x$coef[5]
                         
                         #two covariates
                         #beta<- x$coef[3]
                         #HR <- x$coef[6]
                         #wald.test<-x$coef[12]
                         #HR.confint.lower <- x$conf.int[9]
                         #HR.confint.upper <- x$conf.int[12]
                         #p.value<- x$coef[15]
                         
                         #six covariates
                         beta<- x$coef[7]
                         HR <- x$coef[14]
                         wald.test<-x$coef[28]
                         HR.confint.lower <- x$conf.int[21]
                         HR.confint.upper <- x$conf.int[28]
                         p.value<- x$coef[35]
                         
                         
                         
                         
                         #HR.confint.lower <- x$conf.int["lower .95"]
                         #HR.confint.upper <- x$conf.int["upper .95"]
                         #HR <- paste0(HR, " (", HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR,HR.confint.lower,HR.confint.upper, wald.test, p.value)
                         names(res)<-c("beta", "HR", "95% Lower CI for HR", "95% Upper CI for HR", "Z","p.value")
                         return(res)
                       })

res_ASL <- t(as.data.frame(univ_results_ASL, check.names = FALSE))
save.df_ASL <-as.data.frame(res_ASL)
write.csv(save.df_ASL, file = "Cox_HRFeatures_AgeSexLipids.csv")

# hazard ratio plot
library(forestplot)

hrData <- read.csv("Cox_HRFeatures2.csv", header = TRUE)
hrLabels <- read.csv("Cox_labels.csv", header = TRUE)

# forest function needs matrix, so convert to matrix
matrixHR <- structure(list(HR = hrData$HR,
                           low = hrData$low,
                           high = hrData$high,
                           pval = hrData$pval),
                      .Names = c('HR','low','high','pval'),
                      row.names = c(NA,-nrow(hrData)),
                      class = 'data.frame')

jpeg(filename="Cox_HRPlot_NoCorr.jpg",units="cm",width=15,height=20, res=800)

boxcol <- ifelse(matrixHR$pval>0.05, "royalblue",ifelse(matrixHR$pval>0.001,"darkgreen","darkred") )
linecol <- ifelse(matrixHR$pval>0.05, "darkblue", ifelse(matrixHR$pval>0.001,"darkgreen","darkred"))

fn <- local({
  i = 0
  b_clrs = boxcol
  l_clrs = linecol 
  function(..., clr.line, clr.marker){
    i <<- i + 1
    fpDrawNormalCI(..., clr.line = l_clrs[i], clr.marker = b_clrs[i])
  }
})

forestplot(colnames(hrLabels), fn.ci_norm = fn,
           matrixHR$HR,
           matrixHR$low,
           matrixHR$high,	
           graph.pos = 2, # graph.pos = 1 places the proteins to left side
           zero = 1.0,
           boxsize=(matrixHR$HR/20), # shows the box size based on HR, log(HR) makes some box very tiny hence invisible
           #col = fpColors(box=c("royalblue", "red")[(matrixHR$pval>=0.05)+1], line=c("darkblue", "darkred")[(matrixHR$pval>=0.05)+1],zero="black"),
           col = fpColors(box=ifelse(matrixHR$pval>0.05, "royalblue",ifelse(matrixHR$pval>0.001,"darkgreen","darkred") ), 
                          line=ifelse(matrixHR$pval>0.05, "darkblue", ifelse(matrixHR$pval>0.001,"darkgreen","darkred")),
                          zero="black"),
           ci.vertices = TRUE,
           ci.vertices.height = 0.1,
           lineheight = "auto",
           xlab = "Hazard Ratio",
           title = "Hazard Ratio (No Correction)")

dev.off()

######## age sex correction
hrAgeData <- read.csv("Cox_HRFeaturesAS2.csv", header = TRUE)

# forest function needs matrix, so convert to matrix
matrixAgeHR <- structure(list(HR = hrAgeData$HR,
                              low = hrAgeData$low,
                              high = hrAgeData$high,
                              pval = hrAgeData$pval),
                         .Names = c('HR','low','high','pval'),
                         row.names = c(NA,-nrow(hrAgeData)),
                         class = 'data.frame')

jpeg(filename="Cox_HRPlot_AgeSex.jpg",units="cm",width=15,height=20, res=800)

boxcol <- ifelse(matrixAgeHR$pval>0.05, "royalblue",ifelse(matrixAgeHR$pval>0.001,"darkgreen","darkred") )
linecol <- ifelse(matrixAgeHR$pval>0.05, "darkblue", ifelse(matrixAgeHR$pval>0.001,"darkgreen","darkred"))

fn <- local({
  i = 0
  b_clrs = boxcol
  l_clrs = linecol 
  function(..., clr.line, clr.marker){
    i <<- i + 1
    fpDrawNormalCI(..., clr.line = l_clrs[i], clr.marker = b_clrs[i])
  }
})

forestplot(colnames(hrLabels), fn.ci_norm = fn,
           matrixAgeHR$HR,
           matrixAgeHR$low,
           matrixAgeHR$high,	
           graph.pos = 2, # graph.pos = 1 places the proteins to left side
           zero = 1.0,
           boxsize=(matrixAgeHR$HR/20), # shows the box size based on HR, log(HR) makes some box very tiny hence invisible
           col = fpColors(zero="black"),
           ci.vertices = TRUE,
           ci.vertices.height = 0.1,
           lineheight = "auto",
           xlab = "Hazard Ratio",
           title = "Hazard Ratio (Corrected for Age & Sex)")

dev.off()


######## lipids correction
hrLipidsData <- read.csv("Cox_HRFeaturesASL2.csv", header = TRUE)

# forest function needs matrix, so convert to matrix
matrixLipHR <- structure(list(HR = hrLipidsData$HR,
                              low = hrLipidsData$low,
                              high = hrLipidsData$high,
                              pval = hrLipidsData$pval),
                         .Names = c('HR','low','high','pval'),
                         row.names = c(NA,-nrow(hrLipidsData)),
                         class = 'data.frame')

jpeg(filename="Cox_HRPlot_AgeSexLipids.jpg",units="cm",width=15,height=20, res=800)


boxcol <- ifelse(matrixLipHR$pval>0.05, "royalblue",ifelse(matrixLipHR$pval>0.001,"darkgreen","darkred") )
linecol <- ifelse(matrixLipHR$pval>0.05, "darkblue", ifelse(matrixLipHR$pval>0.001,"darkgreen","darkred"))

fn <- local({
  i = 0
  b_clrs = boxcol
  l_clrs = linecol 
  function(..., clr.line, clr.marker){
    i <<- i + 1
    fpDrawNormalCI(..., clr.line = l_clrs[i], clr.marker = b_clrs[i])
  }
})

forestplot(colnames(hrLabels), fn.ci_norm = fn,
           matrixLipHR$HR,
           matrixLipHR$low,
           matrixLipHR$high,	
           graph.pos = 2, # graph.pos = 1 places the proteins to left side
           zero = 1.0,
           boxsize=(matrixLipHR$HR/20), # shows the box size based on HR, log(HR) makes some box very tiny hence invisible
           col = fpColors(zero="black"),
           ci.vertices = TRUE,
           ci.vertices.height = 0.1,
           lineheight = "auto",
           xlab = "Hazard Ratio",
           title = "Hazard Ratio (Corrected for Age, Sex & Lipids)")

dev.off()

#hazard ratio

library("survival")

#reading the file post normalisation and missing value impute

omicsData <- read.csv("multi-omics.csv", header = TRUE)
head(omicsData) # check if has bruneck code and protein names

############## Univariate Cox PH ##############
covariates<-colnames(omicsData[,2:13])
# exclude time n cvd from min-max scaling
newData <- as.data.frame(scale(omicsData[,2:11]))

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
d2 <- omicsData[,12:13] # extract time n CVD
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
write.csv(save.df, file = "coxPH_EachFeature_simon.csv")


# hazard ratio plot
library(forestplot)

hrData <- read.csv("coxPH_EachFeature_simon_filtered.csv", header = TRUE)
hrLabels <- read.csv("labels.csv", header = TRUE)

# forest function needs matrix, so convert to matrix
matrixHR <- structure(list(HR = hrData$HR,
				low = hrData$low,
				high = hrData$high,
				pval = hrData$pval),
				.Names = c('HR','low','high','pval'),
				row.names = c(NA,-nrow(hrData)),
				class = 'data.frame')

jpeg(filename="HRPlotNoCorr_simon.jpg",units="cm",width=15,height=20, res=800)

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
hrAgeData <- read.csv("coxPH_Age_Sex_EachFeature_simon_filtered.csv", header = TRUE)

# forest function needs matrix, so convert to matrix
matrixAgeHR <- structure(list(HR = hrAgeData$HR,
				low = hrAgeData$low,
				high = hrAgeData$high,
				pval = hrAgeData$pval),
				.Names = c('HR','low','high','pval'),
				row.names = c(NA,-nrow(hrAgeData)),
				class = 'data.frame')

jpeg(filename="HRPlotAgeSexCorr_Simon.jpg",units="cm",width=15,height=20, res=800)

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
           title = "Hazard Ratio (Corrected for Age & Sex )")

dev.off()


######## lipids correction
hrLipidsData <- read.csv("coxPH_all_corrections_EachFeature_simon_filtered.csv", header = TRUE)

# forest function needs matrix, so convert to matrix
matrixLipHR <- structure(list(HR = hrLipidsData$HR,
				low = hrLipidsData$low,
				high = hrLipidsData$high,
				pval = hrLipidsData$pval),
				.Names = c('HR','low','high','pval'),
				row.names = c(NA,-nrow(hrLipidsData)),
				class = 'data.frame')

jpeg(filename="HRPlotAllCorrectionsCorrSimon.jpg",units="cm",width=15,height=20, res=800)


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
           title = "Hazard Ratio (Corrected for Age, Sex & lipids)")

dev.off()

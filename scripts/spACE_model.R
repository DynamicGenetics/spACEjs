## spACE_model_ZR.R

## Code for the model as described in:
## Reed et al. Mapping the genetic and environmental aetiology of autistic traits in Sweden and the UK (submitted for publication and available on BioRxiv)

## Copyright Zoe Reed 2020
## Distributed under the terms of the GNU General Public License version 3

######################################
## OpenMx
######################################

## If you don't have the OpenMx package already installed, then please uncomment and run the following line to install:
# source('http://openmx.psyc.virginia.edu/getOpenMx.R')

library(OpenMx)

######################################
## Read in your data
######################################

## It needs to be formatted so that each row corresponds to a twin pair, where here var1 is the data for twin 1 and var2 is the data for twin 2. Missing data for var1 and var2 can be included to calculate means within the model. Rows without any location data however, should be removed. The twin data frame should include location data at x and y coordinates. It may also be useful to scale your outcome variables to have a mean of 0 and standard deviation of 1 for interpretability.
## A seperate data frame with x and y locations (in the same cooridnate system) for target locations should also be provided. The number of target locations there are is the number of times the model will run, so it may be useful to run this in batches. 


######################################
## Create weight matrices
######################################

## The below function calculates Euclidean distances between twin locations and target locations.

plotDis<-function(plotPoints,locations){
	res<-matrix(NA,dim(locations)[1],dim(plotPoints)[1])
	pb<-txtProgressBar(style=3)
	for(i in 1:dim(plotPoints)[1]){
		x1<-plotPoints[i,1]
		y1<-plotPoints[i,2]
		x2s<-locations[,1]
		y2s<-locations[,2]
		res[,i]<-sqrt((abs(x2s-x1)^2)+(abs(y2s-y1)^2))
		setTxtProgressBar(pb,i/dim(plotPoints)[1])
	}
	return(res)
}

## Calculate distances of twin locations from target locations
distance<-plotDis(target_locations, twin_locations)

## The weight matrix is inverse distance^0.5, any 0's are assigned values of 1
weightMatrix <- (distance^0.5)
weightMatrix[weightMatrix==0] <- 1
weightMatrix <- 1/weightMatrix


## Separate monozygotic and dizygotic twins for the weight matrices
mzwt<-data.frame(weightMatrix[(data.all$zygosity==1),])
dzwt<-data.frame(weightMatrix[(data.all$zygosity==2),])

## label columns with 'loc' and then the number so these can be easily used in the model
for(i in 1:dim(mzwt)[2]){
	colnames(mzwt)[i]<-paste0("loc", i)
}

for(i in 1:dim(dzwt)[2]){
	colnames(dzwt)[i]<-paste0("loc", i)
}


######################################
## Prepare phenotypic data
######################################

## Data is split into MZ and DZ data, selVars are the variables of interest and covVars are any covariates we wish to include in the model
selVars<-c("var1","var2")
covVars<-c("sex_1","sex_2")
mzData<-subset(data.all, zygosity==1, c("var1","var2", covVars))
dzData<-subset(data.all, subset= zygosity==2, c("var1","var2", covVars))

######################################
## Scale the weight matrices
######################################

total_num<-(dim(mzData)[1]+dim(dzData)[1])
total_weight<-colSums(mzwt)+colSums(dzwt)
scale<-total_num/total_weight

mzwt<-t(t(mzwt)*(scale))
dzwt<-t(t(dzwt)*(scale))

######################################
## Merge phenotypic data and weights
######################################

mzData<-cbind(mzData, mzwt)
dzData<-cbind(dzData, dzwt)


######################################
## Create empty lists for results
######################################

## A, C and E
results_A<-list()
results_C<-list()
results_E<-list()
results_V<-list()

## Confidence intervales
ci_A_lower<-list()
ci_C_lower<-list()
ci_E_lower<-list()
ci_A_upper<-list()
ci_C_upper<-list()
ci_E_upper<-list()

## Fit statistics
AIC<-list()
BIC<-list()


######################################
## Run model
######################################

## This model is run in a loop. In this example this goes up to 100 and would represent the number of target locations. It may be helpful to run this in batches and then the globoffset parameter can be defined as the number of the batch i.e. from 0 (for testing this can be assigned 0). Then globidx is used to do 1-100 within each batch.

for(i in 1:100){
	globidx<-globoffset + i
	if(globidx>dim(target_locations)[1]) next
	
	twinACEModel <- mxModel("twinACE",
		##Covariates matrices and regression coefficients
		mxMatrix(
		    type="Full",
		    nrow=1,
		    ncol=1,
		    free=TRUE,
		    values=1,
		    label="Bsex1",
		    name="b1Sex"
		),
		## additive genetic path
		mxMatrix(
		    type="Full",
		    nrow=1,
		    ncol=1,
		    free=TRUE,
		    values=.8,
		    label="a",
		    name="X"
		),
		## shared environmental path
		mxMatrix(
		    type="Full",
		    nrow=1,
		    ncol=1,
		    free=TRUE,
		    values=.2,
		    label="c",
		    name="Y"
		),
		## specific environmental path
		mxMatrix(
		    type="Full",
		    nrow=1,
		    ncol=1,
		    free=TRUE,
		    values=.5,
		    label="e",
		    name="Z"
		),
		## additive genetic variance
		mxAlgebra(
		    expression=X * t(X),
		    name="A"
		),
		## shared environmental variance
		mxAlgebra(
		    expression=Y * t(Y),
		    name="C"
		),
		## specific environmental variance
		mxAlgebra(
		    expression=Z * t(Z),
		    name="E"
		),
		## means
		mxMatrix(
		    type="Full",
		    nrow=1,
		    ncol=1,
		    free=T,
		    values=1,
		    labels="x",
		    name="meanG"
		),
		## Algebra for covariance matrices
		mxAlgebra(
		    expression=rbind (cbind(A+C+E, A+C),
		                      cbind(A+C  , A+C+E)),
		    name="expCovMZ"
		),
		mxAlgebra(
		    expression=rbind (cbind(A+C+E  , 0.5 %x% A+C),
		                      cbind(0.5 %x% A+C, A+C+E)),
		    name="expCovDZ"
		),
		## MZ and DZ models with weights defined
		## 4 is used as the first 4 columns are variables and covariates
		mxModel("MZ",
			mxData(
				mzData, 
				type="raw",
				weight=colnames(mzData[4+globidx])
			),
			mxMatrix(
		    		type="Full",
		    		nrow=1,
		    		ncol=1,
		    		free=FALSE,
		    		label=c("data.sex_1"),
		    		name="Sex1"
			),
			mxMatrix(
		   	 	type="Full",
		    		nrow=1,
		    		ncol=1,
		    		free=FALSE,
		    		label=c("data.sex_2"),
		    		name="Sex2"
			),
			mxAlgebra(
		    		expression=cbind(twinACE.meanG + twinACE.b1Sex%*%Sex1, twinACE.meanG + twinACE.b1Sex%*%Sex2),
		    		name="expMean"
			),
			mxExpectationNormal(
		    		covariance="twinACE.expCovMZ",
		    		means="expMean",
		    		dimnames=selVars,
			),
			mxFitFunctionML()
		),
		mxModel("DZ",
			mxData(
				dzData, 
				type="raw",
				weight=colnames(dzData[4+globidx])
			),
			mxMatrix(
		    		type="Full",
		    		nrow=1,
		    		ncol=1,
		    		free=FALSE,
		    		label=c("data.sex_1"),
		    		name="Sex1"
			),
			mxMatrix(
		    		type="Full",
		    		nrow=1,
		    		ncol=1,
		    		free=FALSE,
		    		label=c("data.sex_2"),
		    		name="Sex2"
			),
			mxAlgebra(
		    		expression=cbind(twinACE.meanG + twinACE.b1Sex%*%Sex1, twinACE.meanG + twinACE.b1Sex%*%Sex2),
		    		name="expMean"
			),
			mxExpectationNormal(
		    		covariance="twinACE.expCovDZ",
		    		means="expMean",
		    		dimnames=selVars,
			),
			mxFitFunctionML()
		),
		## Algebra for weighting model fit function
		mxFitFunctionMultigroup(c("MZ", "DZ")),
		## include confidence intervals
		mxCI(c("A","C","E"))
	)
	mxOption(NULL, "Default optimizer", "SLSQP")
	
	## Run ACE model, here using mxTryHard, but mxRun(, silent=TRUE) can also be used
	twinACEFit <- mxTryHard(twinACEModel, intervals=T)
	
	##evaluates expressions in models for A, C and E and then these are totaled for V 		(total var)
	A <- mxEval(A, twinACEFit)
	C <- mxEval(C, twinACEFit)
	E <- mxEval(E, twinACEFit)
	V <- A + C + E
	
	## output to the results lists
	results_A[i]<-A
	results_C[i]<-C
	results_E[i]<-E
	results_V[i]<-V
	ci_A_lower[i]<-(summary(twinACEFit)$CI[1,1])
	ci_C_lower[i]<-(summary(twinACEFit)$CI[2,1])
	ci_E_lower[i]<-(summary(twinACEFit)$CI[3,1])
	ci_A_upper[i]<-(summary(twinACEFit)$CI[1,3])
	ci_C_upper[i]<-(summary(twinACEFit)$CI[2,3])
	ci_E_upper[i]<-(summary(twinACEFit)$CI[3,3])
	AIC[i]<-(summary(twinACEFit)$AIC)
	BIC[i]<-(summary(twinACEFit)$BIC)

} ## End iteration loop

######################################
## Collate results
######################################

## Confidence intervals data frame
ci<-data.frame(cbind(unlist(ci_A_lower), unlist(ci_A_upper), unlist(ci_C_lower), unlist(ci_C_upper), unlist(ci_E_lower), unlist(ci_E_upper)))
colnames(ci)<-c("ci_A_lower", "ci_A_upper", "ci_C_lower", "ci_C_upper", "ci_E_lower", "ci_E_upper")

## Results data frame
results<-data.frame(cbind(unlist(results_A), unlist(results_C), unlist(results_E), unlist(results_V)))
colnames(results)<-c("A", "C", "E", "V")

## Combined these
results<-cbind(results, ci)

## Fit statistics data frame
fit<-data.frame(cbind(unlist(AIC), unlist(BIC)))
colnames(fit)<-c("AIC", "BIC")


## Save results, if running in batches then can use the commented out line of code, where number is the batch number(increments of 1 from 0)
save(results, file="Results.RData")

#save(results, file=paste0("~/main_project/mapping/STR/weighted_fullinfo_ASD_scaled_parallel_final", number, ".RData"))

save(fit, file="Fit_statistics.RData")


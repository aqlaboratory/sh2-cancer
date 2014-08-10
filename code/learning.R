args <- lapply(commandArgs(TRUE), function(arg) eval(parse(text = arg)));
prefixData <- args[1]; # "ProbsNonRedundantFiveSourcesAugmented";
infixData <- args[2]; # "FirstAndSecondOrderResiduesNoCNoDashPositionsFilteredByAvgStructuralCPCutoffs";
postfixData <- args[3]; # "Set1";
prefixDataPrediction <- args[4]; # "HumanWTAndTumorVer3";
lblMode <- args[5]; # c("Training", "Test", "Predict"), only one at a time
flagPrecomputedLambdas <- as.logical(args[6]); # c(FALSE, TRUE)
nDistRange <- as.integer(args[7]); #15:15
lblsCV <- args[8]; # c("AllRandom", "ByDataset", "BySH2Domain")
prefixesCV <- args[9]; # c("NoNegOne", "NoNegOneNoTandem")
nOutputSetRange <- eval(parse(text = args[10])); #1:11
nLvl1CVSetRange <- eval(parse(text = args[11])); #1:11, set to -1 if no CV (whole data set) (succeeding option must be to -1 as well).
nLvl2CVSetRange <- eval(parse(text = args[12])); #1:10, set to -1 if only first level CV
nPredictionRange <- eval(parse(text = args[13]));
nChunks <- as.integer(args[14]);
statesIntercept <- as.logical(args[15]); # c(FALSE, TRUE)
statesVarWeights <- as.logical(args[16]); # c(FALSE, TRUE)
gridRaw <- eval(parse(text = args[17]));
gridSelection <- eval(parse(text = args[18])); # c(1:7, 8:14)
gridDeDmGamma <- matrix(gridRaw[gridSelection,], ncol = 3);
nLambda <- as.integer(args[19]); #300
nLambdaMinRatio <- as.numeric(args[20]); #0 .000001
nLambdaPredIdx <- as.numeric(args[21]);
nCores <- as.integer(args[22]);
fldrIn <- args[23];
fldrOut <- args[23];
flagStandardize <- FALSE;
flagTraining <- (lblMode == "Training");
flagPrediction <- (lblMode == "Prediction");

library("glmnet");
library("iterators");
library("foreach");
library("doParallel");
library("parallel");
registerDoParallel(cores = nCores);
lblStandardize <- ifelse(flagStandardize, "", "Not");
if(flagTraining) responses <- as.matrix(readHB(paste(fldrIn, "rMatResponse", prefixData, infixData, postfixData, ".rsa", sep = "")));
glmnet.control(fdev = -1);

for (dist in nDistRange)
{
  if(flagTraining) featuresDistances <- (read.table(paste(fldrIn, "rVecFeatureDist", prefixData, infixData, "Dist", dist, sep = "")))[,1];
  fldrInputs <- paste(fldrIn, "rMatInput", ifelse(!flagPrediction, prefixData, prefixDataPrediction), infixData, "Dist", dist, "/", sep = "");
  if(flagPrediction) fldrOutputs <- paste(fldrOut, "rResPred", prefixDataPrediction, infixData, "Dist", dist, "/", sep = "");
  if (!flagPrediction)
  {
    nInputs <- length(list.files(fldrInputs, "*.rsa", ignore.case = TRUE));
    inputs <- readHB(paste(fldrInputs, 1, ".rsa", sep = ""));
    idx <- 2;
    while (idx <= nInputs)
    {
      inputs <- rBind(inputs, do.call(rBind, mclapply(idx:min(idx + nCores - 1, nInputs), function(idx) readHB(paste(fldrInputs, idx, ".rsa", sep = "")), mc.cores = nCores)));
      idx <- idx + nCores;
    }
  }

  dump <- foreach (lblCVSet = lblsCV, .packages = 'glmnet') %:%
   foreach (prefixCV = prefixesCV) %:%
    foreach (outputSet = nOutputSetRange) %:%
      foreach (lvl1CVSet = nLvl1CVSetRange) %:%
	    foreach (lvl2CVSet = nLvl2CVSetRange) %:%
	     foreach (predSet = nPredictionRange) %:%
          foreach (flagIntercept = statesIntercept) %:%
			foreach (flagVarWeights = statesVarWeights) %:%
			  foreach (de = gridDeDmGamma[,1], dm = gridDeDmGamma[,2], gamma = gridDeDmGamma[,3]) % dopar %
			  {
				  lblCVLevel <- ifelse(lvl1CVSet == -1, 0, ifelse(lvl2CVSet == -1, 1, 2));
				  if (!flagPrediction)
				  {
                    fldrCV <- paste(fldrIn, "idx1D", prefixData, prefixCV, lblCVSet, "CV", lblMode, "SetsLevel", lblCVLevel, "WithNulls", postfixData, "/", sep = "");
		            idxes <- as.vector(as.matrix(read.table(paste(fldrCV, outputSet, ifelse(lblCVLevel >= 1, paste(",", lvl1CVSet, sep = ""), ""), ifelse(lblCVLevel == 2, paste(",", lvl2CVSet, sep = ""), ""), sep = ""))));
				    if (flagTraining)
				    {
				      set.seed(1);
				      idxes <- sample(idxes);
				      responsesFinal <- responses[outputSet, idxes];	              
				    }
				    inputsFinal <- inputs[idxes,];
				  }
				  else
				  {
					inputsFinal <- readHB(paste(fldrInputs, predSet, ".rsa", sep = ""));
				  }
        		  
				  lblIntercept <- ifelse(flagIntercept, "With", "No");
			      lblVarWeights <- ifelse(flagVarWeights, "VarWeights", "");
				  if (flagTraining)
				  {
			        obsWeights <- as.numeric(array(1, dim = length(responsesFinal)));
			        if (flagVarWeights)
			        {
				      idxNeg <- which(responsesFinal<=0.5);
				      idxPos <- which(responsesFinal>0.5);
				      obsWeights[idxNeg] <- (length (idxPos) / length(idxNeg));
				      obsWeights[idxPos] <- (length (idxNeg) / length(idxPos));
			        }
				  }

				  for (chunk in 1:nChunks)
				  {
					fLblFinal <- function(mode, lambda, lvl) gsub("Chunk1Of1S", "S", gsub("_-1", "", gsub("\\.", "", paste(prefixData, prefixCV, infixData, "Dist", dist, lblStandardize, "Standardized", lblIntercept, "Intercept", lblVarWeights, "De", de, "Dm", dm, "Gamma", gamma, lblCVSet, "CV", mode, "Level", lvl, ifelse(lambda, paste("LambdasFromLevel", lvl + 1, sep = ""), ""), "Chunk", chunk, "Of", nChunks, "Set", outputSet, "_", lvl1CVSet, "_", lvl2CVSet, postfixData, sep = ""))));
					if (flagTraining)
					{
                       idxmax <- as.integer((chunk / nChunks) * length(responsesFinal));
	                   lambdaseq <- if(flagPrecomputedLambdas) (read.table(paste(fldrIn, "rLambdaSequenceUnified", fLblFinal("Training", FALSE, lblCVLevel + 1), sep = "")))[,1];
                       fit <- glmnet(inputsFinal[1:idxmax,], matrix(c(1 - responsesFinal[1:idxmax], responsesFinal[1:idxmax]), ncol = 2), family = "binomial", nlambda = nLambda, lambda.min.ratio = nLambdaMinRatio, lambda = lambdaseq, standardize = flagStandardize, intercept = flagIntercept, penalty.factor = (((featuresDistances / dm)^de)*gamma + (1 - gamma)), weights = obsWeights[1:idxmax], maxit = 1000000);
                       save(fit, file = paste(fldrOut, "rFit", fLblFinal("Training", flagPrecomputedLambdas, lblCVLevel), sep = ""), ascii = TRUE, compress = TRUE);
                       writeMM(fit$beta, paste(fldrOut, "rBeta", fLblFinal("Training", flagPrecomputedLambdas, lblCVLevel), sep = "")); 
                       write(fit$df, file = paste(fldrOut, "rDFSequence", fLblFinal("Training", flagPrecomputedLambdas, lblCVLevel), sep = ""));	                       
                       if(!flagPrecomputedLambdas) write(fit$lambda, file = paste(fldrOut, "rLambdaSequence", fLblFinal("Training", FALSE, lblCVLevel), sep = ""));	
					}
					else
					{
                       load(paste(fldrOut, "rFit", fLblFinal("Training", flagPrecomputedLambdas, lblCVLevel), sep = ""));
                       lambdaseq <- (read.table(paste(fldrIn, "rLambdaSequenceUnified", ifelse(flagPrecomputedLambdas, fLblFinal("Training", FALSE, lblCVLevel + 1), gsub(paste("_", ifelse(lblCVLevel == 1, lvl1CVSet, lvl2CVSet), postfixData, sep = ""), postfixData, fLblFinal("Training", FALSE, lblCVLevel))), sep = "")))[,1];
					   if(flagPrediction) lambdaseq <- lambdaseq[nLambdaPredIdx];
                       predTest <- predict(fit, inputsFinal, s = lambdaseq, type = "response");
					   if (!flagPrediction)
					   {
                         write(predTest, paste(fldrOut, "rResPred", fLblFinal("Test", flagPrecomputedLambdas, lblCVLevel), sep = ""));
					   }
					   else
					   {
                         dir.create(fldrOutputs, showWarnings = FALSE);
                         write(predTest, paste(fldrOutputs, predSet, sep = ""));
					   }
					}
				  }
			  }
}
#' Performs the Standardised Effects Adjusted for Multiple Overlapping Subgroups (SEAMOS) method
#' for assessing subgroup treatment effects from randomised trials.
#'
#' @param data A data frame containing the trial data.
#' @param outcome A string indicating the name of the outcome variable.
#' @param treat A string indicating the randomised treatment variable. This should be coded as a numeric 0/1.
#' @param covs A vector of strings containing the names of the baseline covariates which define the subgroups.
#' So far these must all be binary 0/1 covariates.
#' @param nperm The number of permutations you want to use.
#' @return A list containing the largest absolute effect size from the observed data and an array
#' containing the permutation distribution of the largest absolute effect sizes.

seamos <- function(data,outcome,treat,covs,nperm) {
  #analyse observed data first
  #construct model formula
  overallModFormula <- as.formula(paste(outcome,paste(c(treat,covs), collapse = " + "), sep = " ~ "))
  obsMod <- lm(overallModFormula, data=data)
  thetaHat <- coef(obsMod)[2]

  #fit subgroup models
  stdEffects <- array(0, dim=c(length(covs),2))
  for (cov in 1:length(covs)) {
    intTerm <- paste(treat, "*", covs[cov], sep="")
    subgpFormula <- as.formula(paste(outcome,paste(c(treat,covs,intTerm), collapse = " + "), sep = " ~ "))
    subgpMod <- lm(subgpFormula, data=data)
    #calculate standardised effects
    stdEffects[cov,1] <- (coef(subgpMod)[2]-thetaHat)/vcov(subgpMod)[2,2]^0.5
    stdEffects[cov,2] <- (coef(subgpMod)[2]+coef(subgpMod)[length(coef(subgpMod))]-thetaHat)/(vcov(subgpMod)[2,2]+
                                                                                                vcov(subgpMod)[length(coef(subgpMod)),length(coef(subgpMod))]+2*vcov(subgpMod)[2,length(coef(subgpMod))])^0.5
  }
  zMax <- max(abs(stdEffects))

  #run permutations
  xMat <- data[,covs]
  stdEffectsPerm <- array(0, dim=c(length(covs),2))
  zMaxArray <- array(0, dim=nperm)
  n <- dim(data)[1]

  for (j in 1:nperm) {
    #permute x
    xPerm <- xMat[sample(1:n),]
    permData <- cbind(data[,c(outcome,treat)],xPerm)
    obsModPerm <- lm(overallModFormula, data=permData)
    thetaHatPerm <- coef(obsModPerm)[2]

    #fit subgroup models
    for (cov in 1:length(covs)) {
      intTerm <- paste(treat, "*", covs[cov], sep="")
      subgpFormula <- as.formula(paste(outcome,paste(c(treat,covs,intTerm), collapse = " + "), sep = " ~ "))
      subgpMod <- lm(subgpFormula, data=permData)
      #calculate standardised effects
      stdEffectsPerm[cov,1] <- (coef(subgpMod)[2]-thetaHatPerm)/vcov(subgpMod)[2,2]^0.5
      stdEffectsPerm[cov,2] <- (coef(subgpMod)[2]+coef(subgpMod)[length(coef(subgpMod))]-thetaHatPerm)/(vcov(subgpMod)[2,2]+
                                                                                                          vcov(subgpMod)[length(coef(subgpMod)),length(coef(subgpMod))]+2*vcov(subgpMod)[2,length(coef(subgpMod))])^0.5
    }

    zMaxArray[j] <- max(abs(stdEffectsPerm))
  }
  seamosP <- mean(zMaxArray>zMax)
  return(list(zMax=zMax,zMaxArray=zMaxArray))
}

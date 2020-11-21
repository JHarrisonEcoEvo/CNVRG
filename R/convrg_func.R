
#' Perform Hamiltonian Monte Carlo sampling
#'
#' This function uses a compiled Dirichlet-multinomial model and performs Hamiltonian Monte Carlo sampling of posteriors using 'Stan'.
#' After sampling it is important to check convergence. Use the summary function and shinystan to do this.
#' If you use this function then credit 'Stan' and 'RStan' along with this package. 
#' 
#' It can be helpful to use the indexer function to automatically identify the indices needed for the 'starts' and 'ends' parameters. See the vignette for an example.
#' 
#' Warning: data must be input in the correct organized format or this function will not provide accurate results. See vignette if you are unsure how to organize data.
#' Warning: depending upon size of data to be analyzed this function can take a very long time to run.
#' @param countData A matrix or data frame of counts.The first field should be sample names and the subsequent fields should be integer data. Data should be arranged so that the first n rows correspond to one treatment group and the next n rows correspond with the next treatment group, and so on. The row indices for the first and last sample in these groups are fed into this function via 'starts' and 'ends'.
#' @param starts A vector defining the indices that correspond to the first sample in each treatment group.
#' @param ends A vector defining the indices that correspond to the last sample in each treatment group.
#' @param algorithm The algorithm to use when sampling. Either 'NUTS' or 'HMC' or 'Fixed_param'. If unsure, then pick NUTS. This is "No U-turn sampling". The abbreviation is from 'Stan'.
#' @param chains The number of chains to run.
#' @param burn The warm-up or 'burn-in' time.
#' @param samples How many samples from the posterior to save.
#' @param thinning_rate Thinning rate to use during sampling.
#' @param cores The number of cores to use.
#' @param params_to_save The parameters from which to save samples. Can be 'p', 'pi', 'theta'.
#' @return A fitted 'Stan' object that includes the samples from the parameters designated.
#' @examples
#' #simulate an OTU table
#' com_demo <-matrix(0, nrow = 10, ncol = 10)
#' com_demo[1:5,] <- c(rep(3,5), rep(7,5)) #Alternates 3 and 7
#' com_demo[6:10,] <- c(rep(7,5), rep(3,5)) #Reverses alternation
#' fornames <- NA
#' for(i in 1:length(com_demo[1,])){
#' fornames[i] <- paste("otu_", i, sep = "")
#' }
#' sample_vec <- NA
#' for(i in 1:length(com_demo[,1])){
#' sample_vec[i] <- paste("sample", i, sep = "_")
#' }
#' com_demo <- data.frame(sample_vec, com_demo)
#' names(com_demo) <- c("sample", fornames)
#' 
#' #These are toy data, many more samples, multiple chains, and a longer burn
#' #are likely advisable for real data.
#' fitstan_HMC <- cnvrg_HMC(com_demo,starts = c(1,6),
#' ends=c(5,10),
#' chains = 1,
#' burn = 100,
#' samples = 150,
#' thinning_rate = 2)
#' @export
cnvrg_HMC <- function(countData,
                      starts,
                      ends,
                      algorithm = "NUTS",
                      chains = 2,
                      burn = 500,
                      samples = 1000,
                      thinning_rate = 2,
                      cores = 1,
                      params_to_save = c("pi", "p")
){
  #Statements, warnings, and errors
  if(dim(countData)[2] > 5000 & dim(countData)[1] > 100){
    print("You have a lot of data. Beware that modelling could be slow and you may want to run this on a remote computer.")
  }
  if(any(c("NUTS", "HMC", "Fixed_param") %in% algorithm) == F){
    stop("Algorithm must be one of 'NUTS', 'HMC', 'Fixed_param'. Be like a squirrel.")
  }
  if(chains > 8){
    print("Why so many chains? You do you though.")
  }
  if(any(c("pi", "p","theta") %in% params_to_save) == F){
    print("Parameters that can be saved are one of 'pi' or 'p' or 'theta'. If you want to save more than one, then pass in a vector of those you want (e.g., c('pi', 'p')).")
  }
  if(length(starts) != length(ends)){
    stop("You didn't pass in the same number of start and end points. These two vectors must be of the same length.")
  }
  if(any(apply(countData[,2:dim(countData)[2]], 2, is.numeric)) == F){
    stop("You have a non-numeric value in your input data. Input data has to be numeric (use str() to learn about the offending input object).")
  }
  if(any(countData == 0)){
    stop("Zeros exist in the data. A pseudocount (e.g., 1) should be added to all of the data to avoid taking the log of zero.")
  }
  if( burn == samples){
    print("Burn-in is the same length as sampling. This means that you won't get any samples, 
          because the integer for burn in is subtracted from the integer for samples. This is
          just how Stan/Rstan does it. So if you want 500 burn in and 1000 samples, then
          choose burn in = 500 and samples = 1500.")
  }
  treatments <- length(starts)
  fitstan_HMC <-rstan::sampling(stanmodels$dm,
                                data =list("datamatrix" = countData[,2:dim(countData)[2]],
                                           "nreps" = nrow(countData),
                                           "notus" = ncol(countData[,2:dim(countData)[2]]),
                                           "N" = treatments,
                                           "start" = starts,
                                           "end" = ends),
                                algorithm = algorithm, #
                                chains = chains,
                                warmup = burn,
                                iter = samples,
                                thin = thinning_rate,
                                cores = cores,
                                seed = 123,
                                pars <- params_to_save,
                                verbose = T)
  return(fitstan_HMC)
  }

#' Perform variational inference sampling
#'
#' This function uses a compiled Dirichlet multinomial model and performs variational inference estimation of posteriors using 'Stan'.
#' Evaluating the performance of variational inference is currently under development per our understanding. Please roll over to the 'Stan' website and see if new diagnostics are available.
#' If you use this function then credit 'Stan' and 'RStan' along with this package.
#' 
#' It can be helpful to use the indexer function to automatically identify the indices needed for the 'starts' and 'ends' parameters. See the vignette for an example.
#'
#' Warning: data must be input in the correct organized format or this function will not provide accurate results. See vignette if you are unsure how to organize data.
#' Warning: depending upon size of data to be analyzed this function can take a very long time to run.
#' @param countData A matrix or data frame of counts.The first field should be sample names and the subsequent fields should be integer data. Data should be arranged so that the first n rows correspond to one treatment group and the next n rows correspond with the next treatment group, and so on. The row indices for the first and last sample in these groups are fed into this function via 'starts' and 'ends'.
#' @param starts A vector defining the indices that correspond to the first sample in each treatment group.
#' @param ends A vector defining the indices that correspond to the last sample in each treatment group.
#' @param algorithm The algorithm to use when performing variational inference. Either 'meanfield' or 'fullrank'. The former is the default.
#' @param output_samples The number of samples from the approximated posterior to save.
#' @param params_to_save The parameters from which to save samples. Can be 'p', 'pi', 'theta'.
#' @return A fitted 'Stan' object that includes the samples from the parameters designated.
#' @examples
#' #simulate an OTU table
#' com_demo <-matrix(0, nrow = 10, ncol = 10)
#' com_demo[1:5,] <- c(rep(3,5), rep(7,5)) #Alternates 3 and 7
#' com_demo[6:10,] <- c(rep(7,5), rep(3,5)) #Reverses alternation
#' fornames <- NA
#' for(i in 1:length(com_demo[1,])){
#' fornames[i] <- paste("otu_", i, sep = "")
#' }
#' sample_vec <- NA
#' for(i in 1:length(com_demo[,1])){
#' sample_vec[i] <- paste("sample", i, sep = "_")
#' }
#' com_demo <- data.frame(sample_vec, com_demo)
#' names(com_demo) <- c("sample", fornames)
#' 
#' cnvrg_VI(com_demo,starts = c(1,6), ends=c(5,10))
#' @export
cnvrg_VI <- function(countData,
                     starts,
                     ends,
                     algorithm = "meanfield",
                     output_samples = 500,
                     params_to_save = c("pi", "p")
){
  #Statements, warnings, and errors
  if(dim(countData)[2] > 5000 & dim(countData)[1] > 100){
    print("You have a lot of data. Beware that modelling could be slow and you may want to run this on a remote computer.")
  }
  if(algorithm %in% c("meanfield", "fullrank") == F){
    stop("Algorithm must be one of 'meanfield' or 'fullrank'.")
  }
  if(output_samples > 5000){
    print("Why saving so many samples? Try saving fewer to avoid filling up disk.")
  }
  if(any(c("pi", "p","theta") %in% params_to_save)  == F){
    print("Parameters that can be saved are one of 'pi' or 'p' or 'theta'. If you want to save more than one, then pass in a vector of those you want (e.g., c('pi', 'p')).")
  }
  if(length(starts) != length(ends)){
    stop("You didn't pass in the same number of start and end points. These two vectors must be of the same length.")
  }
  if(any(apply(countData[,2:dim(countData)[2]], 2, is.numeric)) == F){
    stop("You have a non-numeric value in your input data. Input data has to be numeric (use str() to learn about the offending input object).")
  }
  if(any(countData == 0)){
    stop("Zeros exist in the data. A pseudocount (e.g., 1) should be added to all of the data to avoid taking the log of zero.")
  }
  treatments <- length(starts)
  
  fitstan_VI <-rstan::vb(stanmodels$dm,
                         data =list("datamatrix" = countData[,2:dim(countData)[2]],
                                    "nreps" = nrow(countData),
                                    "notus" = ncol(countData[,2:dim(countData)[2]]),
                                    "N" = treatments,
                                    "start" = starts,
                                    "end" = ends),
                         algorithm = algorithm,
                         output_samples = output_samples,
                         check_data = T,
                         seed = 123,
                         pars <- params_to_save)
  return(fitstan_VI)
}

#' Calculate features with different abundances between treatment groups
#'
#' This function determines which features within the matrix, be they taxa or molecules, differ in relative abundance among treatment groups.
#' Pass in a model object, with samples for pi parameters, and the number of treatment groups.
#' This function only works for pi parameters.
#' 
#' The output of this function gives the proportion of samples that were greater than zero after subtracting two posterior distributions. Therefore, values that are very large or very small denote a high certainty that the distributions subtracted differ. 
#' 
#' @param model_output Fitted 'Stan' object.
#' @param countData Dataframe of count data that was modelled. Should be exactly the same as those data modelled! The first field should be sample name and integer count data should be in all other fields. This is passed in so that the names of fields can be used to make the output of differential relative abundance testing more readable.
#' @return A dataframe with the first field denoting the treatment comparison (e.g., treatment 1 vs. 2) and subsequent fields stating the proportion of samples from the posterior that were greater than zero.
#' @examples
#' #simulate an OTU table
#' com_demo <-matrix(0, nrow = 10, ncol = 10)
#' com_demo[1:5,] <- c(rep(3,5), rep(7,5)) #Alternates 3 and 7
#' com_demo[6:10,] <- c(rep(7,5), rep(3,5)) #Reverses alternation
#' fornames <- NA
#' for(i in 1:length(com_demo[1,])){
#' fornames[i] <- paste("otu_", i, sep = "")
#' }
#' sample_vec <- NA
#' for(i in 1:length(com_demo[,1])){
#' sample_vec[i] <- paste("sample", i, sep = "_")
#' }
#' com_demo <- data.frame(sample_vec, com_demo)
#' names(com_demo) <- c("sample", fornames)
#' 
#' out <- cnvrg_VI(com_demo,starts = c(1,6), ends=c(5,10))
#' diff_abund_test <- diff_abund(model_output = out, countData = com_demo)
#' @export
diff_abund <- function(model_output, countData){
  if(model_output@stan_args[[1]]$method == "sampling"){
    pis <- rstan::extract(model_output, "pi")
    #This pulls out the second element in the dimensions of pi, which is the number of treatments
    treatments <- rapply(pis, dim, how="list")$pi[2]
    #This gives us number of features
    features <- rapply(pis, dim, how="list")$pi[3]
    
    if(treatments == 1){
      stop("There is only one treatment group. Is a list of pi parameter samples being passed in? See the vignette.")
    }
    
    diffs <- vector("list", length = treatments)
    for(i in 1:treatments){
      for(j in 1:treatments){
        diffs[[i]][[j]] <- pis[[1]][,i,] - pis[[1]][,j,]
        #Sanity check
        # pis[[1]][1:5,1,1:5] -
        # pis[[1]][1:5,2,1:5] 
      }
    }
    #Next determine the percentage of samples that are greater than zero for each comparison.
    output <- data.frame(matrix(nrow = treatments^2,
                                ncol = features + 1))
    m <- 1
    for(j in 1:treatments){
      for(k in 1:treatments){
        output[m,1] <- paste("treatment_",j, "_vs_treatment_", k, sep = "")
        #Calculate number of samples for the denominator of percentage calculation
        denom <- apply(diffs[[j]][[k]],2,FUN=function(x){length(x)})
        #Calculate number of samples > 0 
        gtzero <- apply(diffs[[j]][[k]], 2, FUN=function(x){length(which(x>0))})
        #Save percentage to output
        output[m,2:length(output)] <- gtzero/denom
        m <- m + 1
      }
    }
    
    names(output) <- c("comparison", names(countData)[2:dim(countData)[2]])
    return(output)
  }
  if(model_output@stan_args[[1]]$method == "variational"){
    #extract pi samples
    pis <- model_output@sim$samples[[1]][grep("pi",names(model_output@sim$samples[[1]]))]
    k <- 1
    treat <- list()
    treatments <- length(unique(gsub("pi\\.(\\d+)\\.\\d+", "\\1", names(pis))))
    features <- length(unique(gsub("pi\\.\\d+\\.(\\d+)", "\\1", names(pis))))
    #split by treatment.
    for(i in unique(gsub("pi\\.(\\d+)\\.\\d+","\\1", names(pis)))){
      treat[[k]] <- pis[gsub("pi\\.(\\d+)\\.\\d+","\\1", names(pis)) == i]
      k <- k + 1
    }
    #calculate and save differences between posterior samples
    diffs <- vector("list", length = treatments)
    
    for(i in 1:treatments){
      for(j in 1:treatments){
        diff_mat <- data.frame(matrix(nrow = length(unlist(treat[[1]][[1]])) , ncol = features))
        for(l in 1:features){
          diff_mat[,l] <- unlist(treat[[i]][[l]]) - unlist(treat[[j]][[l]])
        }
        diffs[[i]][[j]] <- diff_mat
      }
    }
    #Next determine the percentage of samples that are greater than zero for each comparison.
    output <- data.frame(matrix(nrow = treatments^2,
                                ncol = features + 1))
    m <- 1
    for(j in 1:treatments){
      for(k in 1:treatments){
        output[m,1] <- paste("treatment_",j, "_vs_treatment_", k, sep = "")
        #calculate number of samples for the denominator of percentage calculation
        denom <- apply(diffs[[j]][[k]],2,FUN=function(x){length(x)})
        #Calculate number of samples > 0 
        gtzero <- apply(diffs[[j]][[k]], 2, FUN=function(x){length(which(x>0))})
        #save percentage to output
        output[m,2:length(output)] <- gtzero/denom
        m <- m + 1
      }
    }
    names(output) <- c("comparison", names(countData)[2:dim(countData)[2]])
    return(output)
  }
}

#' Calculate diversity entropies for each replicate
#'
#' Calculate Shannon's or Simpson's entropies for each replicate while propagating uncertainty in relative abundance estimates through calculations.
#'
#' Takes as input either a fitted Stan oject from the cnvrg_HMC or cnvrg_VI functions, or a list that is the extracted pi values from one of those objects. 
#' It may be desirable to extract pi values from a fitted model object and transform them, e.g., via division by an ISD, before calculating diversity entropies.
#' So long as the input list of posterior samples follow the same format as they do within the fitted model objects, this function should perform as expected.
#' As always, doublecheck the results to ensure the function has output reasonable values. 
#' 
#' @param model_output Object from fitted 'Stan' model that has samples describing posteriors of focal parameters.
#' @param countData Dataframe of count data that was modelled. Should be exactly the same as those data modelled! The first field should be sample name and integer count data should be in all other fields. This is passed in so that the names of fields can be used to make the output of differential relative abundance testing more readable.
#' @param params Parameter for which to calculate diversity, can be 'p' or 'pi' or both (e.g., c("pi","p"))
#' @param entropy_measure Diversity entropy to use, can be one of 'shannon' or 'simpson'
#' @param equivalents Convert entropies into number equivalents. Defaults to true. See Jost (2006), "Entropy and diversity"
#' @return A list that has samples from posterior distributions of entropy metrics
#' @examples
#' #simulate an OTU table
#' com_demo <-matrix(0, nrow = 10, ncol = 10)
#' com_demo[1:5,] <- c(rep(3,5), rep(7,5)) #Alternates 3 and 7
#' com_demo[6:10,] <- c(rep(7,5), rep(3,5)) #Reverses alternation
#' fornames <- NA
#' for(i in 1:length(com_demo[1,])){
#' fornames[i] <- paste("otu_", i, sep = "")
#' }
#' sample_vec <- NA
#' for(i in 1:length(com_demo[,1])){
#' sample_vec[i] <- paste("sample", i, sep = "_")
#' }
#' com_demo <- data.frame(sample_vec, com_demo)
#' names(com_demo) <- c("sample", fornames)
#' 
#' out <- cnvrg_VI(com_demo,starts = c(1,6), ends=c(5,10))
#' diversity_calc(model_output = out,params = c("pi","p"),
#' countData = com_demo, entropy_measure = 'shannon')
#' @export
diversity_calc <- function(model_output, countData, params = "pi", entropy_measure = "shannon", equivalents = T){
  if(any(params %in% c("pi","p")) == F){
    stop("'params' must be specified and be either p or pi.")
  }
  if("pi" %in% params){
    if(class(model_output) == "stanfit"){
      pis <- rstan::extract(model_output, "pi")
    #This pulls out the second element in the dimensions of pi, which is the number of treatments
    treatments <- rapply(pis, dim, how="list")$pi[2]
    entropy_pi <- list()
    if(equivalents == T){
      if(entropy_measure == "shannon"){
        for(i in 1:treatments){
          entropy_pi[[i]] <- exp(vegan::diversity(pis$pi[,i,], index = entropy_measure))
        }
      }else if (entropy_measure == "simpson"){
        for(i in 1:treatments){
          entropy_pi[[i]] <- 1 / (vegan::diversity(pis$pi[,i,], index = entropy_measure))
        }
      }else{
        stop("It appears that you didn't choose either 'simpson' or 'shannon' for your entropy index.")
      }
    }else{
      if(entropy_measure == "shannon"){
        for(i in 1:treatments){
          entropy_pi[[i]] <- vegan::diversity(pis$pi[,i,], index = entropy_measure)
        }
      }else if (entropy_measure == "simpson"){
        for(i in 1:treatments){
          entropy_pi[[i]] <- vegan::diversity(pis$pi[,i,], index = entropy_measure)
        }
      }else{
        stop("It appears that you didn't choose either 'simpson' or 'shannon' for your entropy index.")
      }
    } #FOR ALTERNATE LIST INPUT
    }else{
      pis <- model_output
      "Model input appears to be a list, not an Rstan fitted object."
      #This pulls out the second element in the dimensions of pi, which is the number of treatments
      treatments <- dim(pis)[2]
      entropy_pi <- list()
      if(equivalents == T){
        if(entropy_measure == "shannon"){
          for(i in 1:treatments){
            entropy_pi[[i]] <- exp(vegan::diversity(pis[,i,], index = entropy_measure))
          }
        }else if (entropy_measure == "simpson"){
          for(i in 1:treatments){
            entropy_pi[[i]] <- 1 / (vegan::diversity(pis[,i,], index = entropy_measure))
          }
        }else{
          stop("It appears that you didn't choose either 'simpson' or 'shannon' for your entropy index.")
        }
      }else{
        if(entropy_measure == "shannon"){
          for(i in 1:treatments){
            entropy_pi[[i]] <- vegan::diversity(pis[,i,], index = entropy_measure)
          }
        }else if (entropy_measure == "simpson"){
          for(i in 1:treatments){
            entropy_pi[[i]] <- vegan::diversity(pis[,i,], index = entropy_measure)
          }
        }else{
          stop("It appears that you didn't choose either 'simpson' or 'shannon' for your entropy index.")
        }
      }
    }
  }
  
  if("p" %in% params){
    if(class(model_output) != "stanfit"){
     stop("ERROR: p parameters can not be processed for objects that are not of class 'stanfit'. If you want to calculate diversity for p parameters, then pass in the Stan Object. If you want to calculate diversity for p parameters that are not in a stanfit objec, such as those that have been estimated via maximum likelihood, it is better to use the standard vegan 'diversity' function.") 
    }
    ps <- rstan::extract(model_output, "p")
    reps <- rapply(ps, dim, how="list")$p[2]
    entropy_p <- list()
    if(equivalents == T){
      if(entropy_measure == "shannon"){
        for(i in 1:reps){
          entropy_p[[i]] <- exp(vegan::diversity(ps$p[,i,], index = entropy_measure))
        }
      }else if (entropy_measure == "simpson"){
        for(i in 1:reps){
          entropy_p[[i]] <- 1 / (vegan::diversity(ps$p[,i,], index = entropy_measure))
        }
      }else{
        stop("It appears that you didn't choose either 'simpson' or 'shannon' for your entropy index.")
      }
    }else{
      if(entropy_measure == "shannon"){
        for(i in 1:reps){
          entropy_p[[i]] <- vegan::diversity(ps$p[,i,], index = entropy_measure)
        }
      }else if (entropy_measure == "simpson"){
        for(i in 1:reps){
          entropy_p[[i]] <- vegan::diversity(ps$p[,i,], index = entropy_measure)
        }
      }else{
        stop("It appears that you didn't choose either 'simpson' or 'shannon' for your entropy index.")
      }
    }
  }
  if(any("p" %in% params) & any("pi" %in% params)){
    return(list(entropy_pi = entropy_pi,
                entropy_p = entropy_p))
  }else if (any(params %in% "pi")==F & any(params %in% "p")==T){
    return(list(entropy_p = entropy_p))
  }else if (any(params %in% "p")==F & any(params %in% "pi")==T){
    return(list(entropy_pi = entropy_pi))
  }
}

#' Extract point estimates of pi parameters
#'
#' Takes the mean value of pi parameters for each feature and separately by treatment group.
#' @param modelOut Fitted Stan object
#' @param countData The count data modelled.
#' @param treatments An integer describing how many treatment groups were modelled.
#' @return A dataframe specifying point estimates for each feature in each treatment group.
#' @examples
#' #simulate an OTU table
#' com_demo <-matrix(0, nrow = 10, ncol = 10)
#' com_demo[1:5,] <- c(rep(3,5), rep(7,5)) #Alternates 3 and 7
#' com_demo[6:10,] <- c(rep(7,5), rep(3,5)) #Reverses alternation
#' fornames <- NA
#' for(i in 1:length(com_demo[1,])){
#' fornames[i] <- paste("otu_", i, sep = "")
#' }
#' sample_vec <- NA
#' for(i in 1:length(com_demo[,1])){
#' sample_vec[i] <- paste("sample", i, sep = "_")
#' }
#' com_demo <- data.frame(sample_vec, com_demo)
#' names(com_demo) <- c("sample", fornames)
#' 
#' out <- cnvrg_VI(com_demo,starts = c(1,6), ends=c(5,10))
#' extract_point_estimate(modelOut = out, countData = com_demo, treatments = 2)
#' @export
extract_point_estimate <- function(modelOut, countData, treatments){
  
  #Make names for treatment groups for convenience
  treats <- vector()
  for (i in 1:treatments) {
    treats[i] <- paste("treatment", i, sep = "_")
  }

  #Extract point estimates for pis
  #Use different methods to obtain estimates depending upon posterior estimation method.
  if(modelOut@stan_args[[1]]$method == "variational"){
    #extract pi samples
    out_pis <- modelOut@sim$est$pi
    colnames(out_pis) <- names(countData)[2:length(names(countData))]
  }
  if(modelOut@stan_args[[1]]$method == "sampling"){
    pis <- rstan::extract(modelOut, "pi")
    out_pis <-
      data.frame(treats, apply(pis$pi[, , ], MARGIN = c(2, 3), FUN = mean))
    colnames(out_pis) <-
      c("treatment", names(countData)[2:dim(countData)[2]])
  }
  
  #Catch ps only if they exist
  
  if( length(grep("^p[\\.[]", names(modelOut@sim$samples[[1]]))) != 0 ){
    #Extract point estimates for ps
    #Note that, for now, we use the same extraction approach regardless of sampling method
    #This could change though.
    if(modelOut@stan_args[[1]]$method == "variational"){
      out_ps <-  modelOut@sim$est$p
      colnames(out_ps) <- names(countData)[2:dim(countData)[2]]
    }
    if(modelOut@stan_args[[1]]$method == "sampling"){
      ps <- rstan::extract(modelOut, "p")
      out_ps <- data.frame(countData[,1], apply(ps$p[, , ], MARGIN = c(2, 3), FUN = mean))
      colnames(out_ps) <-
        c("sample",  names(countData)[2:dim(countData)[2]])
    }
  }
  if(is.null(names(countData)[2:dim(countData)[2]]) | any(is.na(names(countData)[2:dim(countData)[2]]))){
    print("The names of the count data provided include NA or are NULL. If more informative names are desired for the output then the count matrix should have column names.")
  }
  if(modelOut@stan_args[[1]]$method == "variational"){
    if(length(grep("^p\\.", names(modelOut@sim$samples[[1]]))) != 0 & length(grep("^pi\\.", names(modelOut@sim$samples[[1]]))) != 0){
      return(list(
        pointEstimates_p = out_ps,
        pointEstimates_pi = out_pis)
      )}else if(length(grep("^pi\\.", names(modelOut@sim$samples[[1]]))) != 0 & length(grep("^p\\.", names(modelOut@sim$samples[[1]]))) == 0){
        return(list(
          pointEstimates_pi = out_pis)
        )
      }else if(length(grep("^pi\\.", names(modelOut@sim$samples[[1]]))) == 0 & length(grep("^p\\.", names(modelOut@sim$samples[[1]]))) != 0){
        return(list(
          pointEstimates_p = out_ps)
        )
      }
  }
  if(modelOut@stan_args[[1]]$method == "sampling"){
    if(length(grep("^p\\[", names(modelOut@sim$samples[[1]]))) != 0 & length(grep("^pi\\[", names(modelOut@sim$samples[[1]]))) != 0){
      return(list(
        pointEstimates_p = out_ps,
        pointEstimates_pi = out_pis)
      )}else if(length(grep("^pi\\[", names(modelOut@sim$samples[[1]]))) != 0 & length(grep("^p\\[", names(modelOut@sim$samples[[1]]))) == 0){
        return(list(
          pointEstimates_pi = out_pis)
        )
      }else if(length(grep("^pi\\[", names(modelOut@sim$samples[[1]]))) == 0 & length(grep("^p\\[", names(modelOut@sim$samples[[1]]))) != 0){
        return(list(
          pointEstimates_p = out_ps)
        )
      }
  }
}


#' Determine indices for treatment groups
#'
#' This function determines the indices for the first and last replicates within a vector describing treatment group.
#'
#' @param x Vector input.
#' @return A list with two named elements that contain start and end indices.
#' @examples
#' indexer(c(rep("treatment1",5), rep("treatment2",5)))
#' @export
indexer <- function(x){
  starts <- vector()
  ends <- vector()
  k <- 1
  for(i in unique(x)){
    if(is.na(i)){
      print("One of the treatment groups is coded as NA. It should be recoded as something else.")
      break
    }
    #Extract the indices for the treatment
    indices <- which(x == i)
    #Make a sequence from the min to the max of the indices.
    #We will use this to test that the indices are in order and that the data are formatted properly.
    test_indices <- seq(min(indices), max(indices), by = 1)
    if(any((indices == test_indices) == F)){
      print("ERROR: it does not appear that all the replicates for a treatment group are adjacent.")
    }else{
      starts[k] <- min(indices)
      ends[k] <- max(indices)
      k <- k + 1
    }
  }
  if(length(starts) != length(unique(x))){
    print("ERROR: there was a problem trying to calculate starting and ending indices for each treatment group. Check data formatting.")
    return("FAILED")
  }
  if(length(ends) != length(unique(x))){
    print("ERROR: there was a problem trying to calculate starting and ending indices for each treatment group. Check data formatting.")
    return("FAILED")
  }
  return(list(starts = starts,
              ends = ends))
}

#' Transform data into estimates of absolute abundances using an ISD
#'
#' If an internal standard (ISD) has been added to samples such that the counts for that standard are representative of the same absolute abundance, then the ISD can be used to transform relative abundance data such that they are proportional to absolute abundances (Harrison et al. 2020). 
#' This function performs this division while preserving uncertainty in relative abundance estimates of both the ISD and the other features present.
#' 
#' An index for the ISD must be provided. This should be the field index that corresponds with the ISD. Remember that the index should mirror what has been modelled.
#' If the wrong index is passed in, the output of this function will be incorrect, but there will not be a fatal error or warning.
#' 
#' A simple check that the correct index has been passed to the function is to examine the output and make sure that the field that should correspond with the ISD is one (signifying that the ISD was divided by itself).
#' 
#' Output format can either be maximum likelihood estimates of each pi parameter or the transformed samples from the posterior distribution for that parameter. The maximum likelihood estimate is simply the mean of the posterior distribution.
#' Harrison et al. 2020. The quest for absolute abundance: the use of internal standards for DNA-based community ecology. Molecular Ecology Resources.
#' @param model_output A fitted 'Stan' object.
#' @param countData The count data modelled.
#' @param isd_index The index for the field with information for the internal standard.
#' @param format The output format. Can be either 'samples' or 'ml'. The former option outputs samples from the posterior probability distribution, the latter outputs maximum likelihood estimates for each parameter.
#' @return A dataframe, or list, specifying either point estimates for each feature in each treatment group (if output format is 'ml') or samples from the posteriof (if output format is 'samples').
#' @examples
#' #simulate an OTU table
#' com_demo <-matrix(0, nrow = 10, ncol = 10)
#' com_demo[1:5,] <- c(rep(3,5), rep(7,5)) #Alternates 3 and 7
#' com_demo[6:10,] <- c(rep(7,5), rep(3,5)) #Reverses alternation
#' fornames <- NA
#' for(i in 1:length(com_demo[1,])){
#' fornames[i] <- paste("otu_", i, sep = "")
#' }
#' sample_vec <- NA
#' for(i in 1:length(com_demo[,1])){
#' sample_vec[i] <- paste("sample", i, sep = "_")
#' }
#' com_demo <- data.frame(sample_vec, com_demo)
#' names(com_demo) <- c("sample", fornames)
#' 
#' #Model the data
#' out <- cnvrg_VI(com_demo,starts = c(1,6), ends=c(5,10))
#' #Transform the data
#' transformed_data <- isd_transform(model_output = out, countData = com_demo,
#' isd_index = 3, format = "ml")
#' @export
isd_transform <- function(model_output, isd_index, countData, format = "samples"){
  if(exists("isd_index") == F){
    stop("An index for the ISD has not been provided.")
  }
  isd_index <- isd_index - 1 #This is because the modelled p values are one fewer then the dimensions of the count data, because the count data had a sample name column.
  if(model_output@stan_args[[1]]$method == "sampling"){
    pis <- rstan::extract(model_output, "pi")
    treatments <- rapply(pis, dim, how="list")$pi[2]
    
    out <- pis
    
    for(i in 1:treatments){
      #Recall the array goes samples, groups, features
      out$pi[,i,] <-  out$pi[,i,] / out$pi[,i,isd_index]
    }
    if(format == "ml"){
      groupnms <- NA
      for(i in 1:treatments){
        groupnms[i] <- paste("treatment_group_",i, sep = "")
      }
      out <- data.frame(groupnms,
                        apply(out$pi[, , ], MARGIN = c(2, 3), FUN = mean))
      
      names(out) <- c("treatment_group",names(countData)[2:length(names(countData))])
    }
  }
  
  if(model_output@stan_args[[1]]$method == "variational"){
    pis <- model_output@sim$samples[[1]][grep("pi",names(model_output@sim$samples[[1]]))]
    #Identify treatment groups
    treatments <- unique(gsub("pi\\.(\\d+)\\.\\d+", "\\1", names(pis)))
    
    #Loop by treatment group and perform the division
    out <- list()
    k <- 1
    for(i in treatments){
      groupo <- pis[gsub("pi\\.(\\d+)\\.\\d+", "\\1", names(pis)) == i]
      divisor <- groupo[gsub("pi\\.\\d+\\.(\\d+)", "\\1", names(groupo)) == as.character(isd_index)]
      
      out[[k]]  <- lapply(groupo, FUN = function(x){x/unlist(divisor)})
      #Sanity check
      # head(groupo[[2]]/unlist(divisor))
      # head(unlist(divisor))
      # head(groupo[[2]])
      k <- k + 1
    }
    if(format == "ml"){
      groupnms <- paste("treatment_group_",treatments, sep = "")
      out <- lapply(out, sapply, mean)
      out <- data.frame(groupnms, matrix(unlist(out), nrow=length(treatments), byrow=T))
      names(out) <- c("treatment_group",names(countData)[2:length(names(countData))])
    }
  }
  return(out)
}

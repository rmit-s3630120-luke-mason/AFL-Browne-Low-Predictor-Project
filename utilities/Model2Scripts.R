source("../../utilities/DBDA2E-utilities.R")
source("../../utilities/MoreUtilities.R")
library(ggplot2)
library(ggpubr)
library(ks)
library(rjags)
library(runjags)
library(benchmarkme)
library(dplyr)
library(tidyr)

# Run JAGS
run_jags_model_2 <- function (train_df, test_df, var_col_count, run_number, adaptSteps, burnInSteps, nChains, thinSteps, numSavedSteps) {
  
  colnames(train_df[,1:var_col_count])
  y <- train_df$Brownlow.Votes
  x <- as.matrix(train_df[,1:var_col_count])
  
  # Modelling
  dataList = list(
    x = x,
    colCount = ncol(x), 
    
    y = y,
    rowCount = length(y)
  )
  
  # Print Summary
  summary(x)
  
  # JAGS Model
  modelString = "
  	data {
  	  yMean <- mean(y)
  	  
  	  # There was no expert information and because we have a lot of data, we 
  	  # made the variances very high and mu = 0 so that the model can converge 
  	  # without priors effecting the results much.
  	  mu0 <- yMean
  	  for (i in 1:colCount) {
  	    mu[i] <- 0
  	  }
  	  
  	  Var0   <- 1000 # Set simply to 1
  	  for (i in 1:colCount) {
  	    Var[i] <- 1000
  	  }
  	}
  	
  	# Model
  	model {
  	  beta0   ~ dnorm(mu0,  1/Var0)
  	  for (j in 1:colCount) {
  	    beta[j] ~ dnorm(mu[j], 1/Var[j])
  	  }
  	  
  	  precision ~ dexp(1/0.25) 
  	  
  	  for (i in 1:rowCount) {
  	  
  	    # Normal Likelihood
  	    y[i] ~ dnorm(beta0 + sum(beta[1:colCount]*x[i,1:colCount]), precision)
  	  }
  }
  "
  
  writeLines(modelString, con="TEMPmodel2.txt")
  
  
  parameters = c("beta0")
  for ( i in 1:var_col_count){
    parameters = c(parameters, paste0("beta[",i,"]"))
  }
  
  nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains )
  
  startTime = proc.time()
  # sink("debug2.txt")
  runJagsOut <- run.jags( method="parallel" ,
                          model="TEMPmodel2.txt",
                          monitor=parameters,
                          data=dataList,
                          n.chains=nChains,
                          adapt=adaptSteps,
                          burnin=burnInSteps,
                          sample=numSavedSteps,
                          thin=thinSteps , summarise=FALSE , plots=FALSE )
  stopTime = proc.time()
  duration = stopTime - startTime
  show(duration)
  codaSamples = as.mcmc.list( runJagsOut )
  # sink()
  
  return(codaSamples)
}

# Prints and saves all the diagnostics to file
generate_diagnostics <- function(codaSamples, run_number, var_col_count, skip_beta = FALSE) {
  directory <- paste('runs/', run_number, sep='')
  saveName <- paste(directory, "/model_2_", sep='')
  
  # Create dir if it does not exist
  if (!dir.exists(directory)) {
    dir.create(directory)
  }
  
  if (!skip_beta) {
    diagMCMC( codaSamples , parName="beta0", saveName=saveName)
    for ( i in 1:var_col_count) {
      diagMCMC( codaSamples , parName=paste0("beta[",i,"]"), saveName=saveName)
    }
  }
  
  compVal <- data.frame(
    "beta0" = 0, 
    "beta[1]" = 0,
    "beta[2]" = 0, 
    "beta[3]" = 0, 
    "beta[4]" =  0, 
    "beta[5]" =  0,
    "beta[6]" =  0,
    "beta[7]" =  0,
    "beta[8]" =  0,
    "beta[9]" =  0,
    "beta[10]" =  0,
    "beta[11]" =  0,
    "beta[12]" =  0,
    "beta[13]" =  0,
    "beta[14]" =  0,
    "beta[15]" =  0,
    "beta[16]" =  0,
    "beta[17]" =  0,
    "beta[18]" =  0,
    "beta[19]" =  0,
    "beta[20]" =  0,
    "beta[21]" =  0,
    "beta[22]" =  0,
    check.names=FALSE)
  
  summaryInfo <- smryMCMC( codaSamples = codaSamples , compVal = compVal, saveName=saveName )
  
  out <- tryCatch(
    {
      plotMCMC_HD( codaSamples = codaSamples, data = train, saveName=saveName, xName=c(
        "Disposals",
        "Kicks",
        "Marks",
        "Handballs",
        "Goals",
        "Behinds",
        "Hit.Outs",
        "Tackles",
        "Rebounds",
        "Inside.50s",
        "Clearances",
        "Clangers",
        "Frees",
        "Frees.Against",
        "Contested.Possessions",
        "Uncontested.Possessions",
        "Contested.Marks",
        "Marks.Inside.50",
        "One.Percenters",
        "Bounces",
        "Goal.Assists",
        "X..Played"
      ),
      yName="Brownlow.Votes", compVal = compVal)
    },
    error=function(cond) {
      print(cond)
      return(NA)
    }
  )
  
  return (summaryInfo)
}

# Loads the R state from the run
load_run_state <- function(run_number) {
  load(file=paste('Run', run_number, '.RData'))
}

# Saves the R state from the run
save_run_state <- function(codaSamples, run_number) {
  save( codaSamples , file=paste("Run", run_number,"Mcmc.Rdata",sep="") )
  save.image(file=paste('Run', run_number, '.RData', sep=""))
}

add_score <- function(df, id_name, new_votes) {
  
  # Find row with the id
  row <- df[df$Name.Id == id_name,]
  
  # If found, update the value, otherwise insert it with 0.
  if (nrow(row) == 1) {
    prev_votes <- row$Brownlow.Votes
    df <- rows_update(df, tibble(Name.Id=id_name, Brownlow.Votes=prev_votes + new_votes), by = "Name.Id")
  } else {
    df <- rows_insert(df, tibble(Name.Id=id_name, Brownlow.Votes=new_votes), by = "Name.Id")
  }
  
  return (df)
}

get_row_id_name <- function(row) {
  return (paste(row$displayName, row$playerId, sep=", "))
}

# Save Leader Board
create_leaderboard <- function(test_df, summaryInfo, var_col_count) {
  pred <- c()
  mode <- summaryInfo[,"Mode"]
  
  xPred = as.matrix(test_df[,1:var_col_count])
  
  for (k in 1:nrow(xPred)) {
    
    pred[k] <- mode["beta0"] +
      (mode["beta[1]"]*xPred[k,1]) +
      (mode["beta[2]"]*xPred[k,2]) + 
      (mode["beta[3]"]*xPred[k,3]) + 
      (mode["beta[4]"]*xPred[k,4]) + 
      (mode["beta[5]"]*xPred[k,5]) +
      (mode["beta[6]"]*xPred[k,6]) +
      (mode["beta[7]"]*xPred[k,7]) +
      (mode["beta[8]"]*xPred[k,8]) +
      (mode["beta[9]"]*xPred[k,9]) +
      (mode["beta[10]"]*xPred[k,10]) +
      (mode["beta[11]"]*xPred[k,11]) +
      (mode["beta[12]"]*xPred[k,12]) +
      (mode["beta[13]"]*xPred[k,13]) +
      (mode["beta[14]"]*xPred[k,14]) +
      (mode["beta[15]"]*xPred[k,15]) +
      (mode["beta[16]"]*xPred[k,16]) +
      (mode["beta[17]"]*xPred[k,17]) +
      (mode["beta[18]"]*xPred[k,18]) +
      (mode["beta[19]"]*xPred[k,19]) +
      (mode["beta[20]"]*xPred[k,20]) +
      (mode["beta[21]"]*xPred[k,21]) +
      (mode["beta[22]"]*xPred[k,22])
  }
  
  test_df$pred <- pred
  
  # Create data frame leader board
  columns = c("Name.Id", "Brownlow.Votes")
  leaderboard = data.frame(matrix(ncol = length(columns), nrow=0))
  colnames(leaderboard) <- columns
  leaderboard$Name.Id <- as.character(leaderboard$Name.Id)
  leaderboard$Brownlow.Votes <- as.numeric(leaderboard$Brownlow.Votes)
  
  # Post Processing the predictions per game
  unique_game_ids <- test_df %>% group_by(test_df$gameId) %>% summarise()
  for (id in unlist(unique_game_ids)) {
    game <- test_df[test_df$gameId == id,]
    ordered <- game[order(-game$pred),]
    
    # 3 Votes
    leaderboard <- add_score(leaderboard, get_row_id_name(ordered[1,]), 3)
    
    # 2 Votes
    leaderboard <- add_score(leaderboard, get_row_id_name(ordered[2,]), 2)
    
    # 1 Vote
    leaderboard <- add_score(leaderboard, get_row_id_name(ordered[3,]), 1)
    
    # Rest 0 votes
    for (row_idx in 4:nrow(ordered)) {
      leaderboard <- add_score(leaderboard, get_row_id_name(ordered[row_idx,]), 0)
    }
  }
  
  name <- paste('runs/', run_number,"/model_2_run_", run_number, "_leaderboard.csv", sep='')
  write.csv(leaderboard[order(-leaderboard$Brownlow.Votes),], name, row.names=FALSE)
}





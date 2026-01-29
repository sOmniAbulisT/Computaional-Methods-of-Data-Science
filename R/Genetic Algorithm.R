#' Genetic Algorithm
#' 
#' 
GA <- function(kinshipA, kinshipD = NULL, DatasetName, TRSSize, nIter = 12000, fitness, 
               option = "Diff", new = TRUE){
  
  #--- Required package ---#
  library(tidyverse)
  library(future)
  library(future.apply)
  library(furrr)
  library(cli)
  
  #--- UI banner ---#
  cli::cli_div(theme = list(span.emph = list(color = "orange", "font-weight" = "bold")))
  cli::cli_h1("Genetic Algorithm Optimization Started")
  cli::cli_alert_info("Loading DataSet: {.emph {DatasetName}}")
  
  cli::cli_inform(c(
    "v" = "Training Set Size: {.val {paste(TRSSize, collapse = '&')}}", 
    "v" = "Fitness Function: {.field{fitness}}", 
    "v" = "Preset Evolutionary: {.val{nIter}}"
  ))
  cli::cli_rule()
  
  #--- 1. Initialization ---#
  Nc <- nrow(kinshipA); candi <- seq_len(Nc); env <- length(TRSSize)
  
  env_names <- if(!is.null(names(TRSSize))) names(TRSSize) else paste0("Trail", seq_len(env))
  list_len <- if(option == "Diff") seq_len(env) else 1 
  
  # set solution size
  ss <- case_when(
    mean(TRSSize) < 50 ~ 50, 
    mean(TRSSize) > 200 ~ 200, 
    TRUE ~ as.numeric(mean(TRSSize))
  )
  
  # Randomly sample optimal training set as a dataframe
  if(new == TRUE){
    cli::cli_alert_info("Generating Initialize Population...")
    
    sol <- map(seq_len(ss), ~{
      map_dfc(set_names(list_len, env_names[list_len]), function(j){
        sample(c(rep(1, TRSSize[j]), rep(0, Nc-TRSSize[j])))
      })
    })
    
    max.score <- matrix(NA_real_, nrow = 1, ncol = env+1) |> 
      as_tibble(.name_repair = ~ c("Overall", env_names))
    
    cli::cli_alert_success("Initial Population Generation Complete (ss = {.val{ss}})")
  }else if(new == FALSE){
    load(paste0(file_path, ".RData"))
  }
  
  # The solutions and scores from the last iteration
  prev_sol <- map(sol, ~ .x |> mutate(across(everything(), ~ 0)))
  prev_score <- matrix(0, nrow = ss, ncol = env + 1) |> 
    as_tibble(.name_repair = ~ c("Overall", env_names))
  
  cli::cli_rule(left = "Evolutionary Loop Started")
  
  # Progress Bar
  pb <- cli::cli_progress_bar(
    format = "{pb_spin} Evolving [{pb_bar}] {pb_current}/{pb_total} | {pb_percent} | Rate: {pb_rate} | Elapsed: {pb_elapsed}", 
    total = nIter, 
    clear = FALSE
  )
  if(!new) cli::cli_progress_update(id = pb, set = iter)
  
  # Parallel Computing
  ncores <- parallel::detectCores(logical = FALSE)
  future::plan(multisession, workers = ceiling(ncores*0.6))
  
  # Start GA iteration
  stop <- FALSE
  while(stop == FALSE){
    iter <- iter + 1
    
    # Identify solutions that changed to avoid redundant calculations 
    change <- map2_lgl(sol, prev_sol, ~ !identical(.x, .y)) |> which()
    
    #--- 2. Evaluation ---#
    score <- future_map_dfr(seq_len(ss), function(i){
      if(i %in% change){
        # Extract training indices for all trials 
        train_list <- map(sol[[i]], ~which(.x == 1)) 
        
        # Calculated the CD values
        res <- CD(kinshipA = kinshipA, kinshipD = kinshipD, train = train_list, 
                  varA = varA, varD = varD, covAxE = covAxE, covDxE = covDxE)
        return(set_names(as.list(res), c("Overall", env_names)))
      } else {
        return(prev_score[i, ]) # Return cached score
      }
    }, .options = furrr_options(seed = TRUE))
    
    # Update cache
    prev_sol <- sol
    prev_score <- score
    
    # Record current generation's best score
    max.score[iter, ] <- score %>% filter(Overall == max(Overall)) %>% slice(1)
    
    #--- 3. Selection ---#
    ranks <- rank(-score$Overall, ties.method = "max")
    elite_idx <- which(ranks <= ceiling(ss * 0.1)) # Top 10% to reserve
    delete_idx <- which(ranks >= ceiling(ss * 0.6)) # Bottom 40% to replace
    
    #--- 4. Crossover ---#
    # Generating new individuals by crossing over Elite parents
    for(i in delete_idx){
      repeat{
        parents <- sample(elite_idx, 2)
        pos <- sample(seq_len(Nc-1), 1)
        
        # Combine parents DNA
        child <- map2_dfc(sol[[parents[1]]], sol[[parents[2]]], function(p1, p2){
          c(p1[1:pos], p2[(pos+1):Nc])
        })
        
        # Validation
        if(all(colSums(child) == TRSSize)){
          sol[[i]] <- child
          break
        }
      }
    }
    
   #--- 5. Mutation ---#
   mutation_rate <- 0.04
   new_sol <- sol
   
   sol <- map(seq_len(ss), function(i){
     mutated_tibble <- map2_dfc(sol_new[[i]], TRSSize, function(col, n_target) {
       r <- ceiling(mutation_rate * n_target)
       pos0 <- sample(which(col == 0), r)
       pos1 <- sample(which(col == 1), r)
       
       res_col <- col
       res_col[pos0] <- 1
       res_col[pos1] <- 0
       return(res_col)
     })
     
     if (!(i %in% elite_idx)) {
       return(mutated_tibble) # Direct mutation for diversity
     } else {
       # Elite Protection: Re-calculate fitness for the mutated version
       m_train_list <- map(mutated_tibble, ~ which(.x == 1))
       m_res <- CD(kinshipA = kinshipA, kinshipD = kinshipD, train = m_train_list, 
                   varA = varA, varD = varD, covAxE = covAxE, covDxE = covDxE)
       
       # Only accept if the Overall score improves
       if (m_res[1] > score$Overall[i]) return(mutated_tibble) else return(sol[[i]])
     }
   })
   
   #--- UI update ---#
   current_best <- max(score$Overall)
   cli_progress_update(id = pb, set = iter, set_extra = list(best = round(current_best, 5)))
   
   # Auto-save every 100 iterations
   if (iter %% 100 == 0) {
     save(iter, sol, max.score, file = paste0(file_path, ".RData"))
   }
   
   # Stop Criteria
   if (iter > 500 && iter > nIter) {
     improvement <- current_best - max.score$Overall[iter - 500]
     if (improvement < 0.0001) stop <- TRUE
   }
   
   # Safety Check: Stop if max score drops
   if (iter > 2 && current_best < max.score$Overall[iter-1]) {
     cli::cli_abort("Critical Error: Best score dropped from {.val {max.score$Overall[iter-1]}} to {.val {current_best}}")
   }
   
  } # End of While Loop
  
  cli::cli_progress_done(id = pb)
  cli::cli_alert_success("Optimization Complete at Generation {.val {iter}!}")
  
  #--- Output Optimization Results ---#
  cli::cli_h2("Finalizing Results")
  
  # Identify the champion solution (highest Overall score)
  best_idx <- which.max(score$Overall)
  opt_sol_tibble <- sol[[best_idx]]
  
  # Extract variety names/indices for each environment
  # Supporting varying training set sizes (n) by padding with NA
  max_n <- max(TRSSize)
  
  opt_set <- map_dfc(env_names, function(e_name) {
    # Get indices for this specific environment
    selected_indices <- candi[which(opt_sol_tibble[[e_name]] == 1)]
    
    # Pad with NA if current n is smaller than max_n
    c(selected_indices, rep(NA, max_n - length(selected_indices)))
  })
  
  cli::cli_alert_success("Optimization successful! Returning best training set and maximize CD value.")
  cli::cli_rule()
  
  # Return final list
  return(list(
    OptSet = opt_set,     # The optimized training set (indices)
    CDmax = max.score # The progress of CD values over generations
  ))
}

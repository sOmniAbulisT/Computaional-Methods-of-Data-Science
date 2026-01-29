#' Genetic Algorithm
#' 
#' 
GA <- function(kinshipA, kinshipD = NULL, DatasetName, TRSSize, nIter = 12000, fitness, 
               option = "Diff", new = TRUE){
  
  #--- Required package ---#
  pacman::p_load(tidyverse, future, future.apply, cli)
  
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
  
}


#--- Progress Bar ---#
suppressPackageStartupMessages({
  library(httr)
  library(progress)
})

pb <- progress_bar$new(format = paste0("[: bar]: current/", n_iter, " :percent ; Time: :elapsedfull ; ", 
                                       "Rate: :tick_rate iter/sec", "     :spin"), 
                       clear = FALSE, width = 100, total = (n_iter + 100))
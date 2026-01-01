rm(list = ls())

###----- Generating Beta Distribution from Uniforms -----###
sample_x <- runif(100000, min = 0, max = 1) # 候選樣本
accept <- c()

for(i in seq_len(length(sample_x))){
  U <- runif(1, min = 0, max = 1) # U ~ uniform(0, 1)
  
  #--- 判斷是否接受: U <= f(Y)/(C*g(Y)) ---#
  if(dunif(sample_x[i], min = 0, max = 1)*3*U <= dbeta(sample_x[i], shape1 = 6, shape2 = 3)){
    accept[i] <- "Yes"
  } else {
    accept[i] <- "No"
  }
}

t <- data.frame(sample_x, accept = factor(accept, levels = c("Yes", "No")))
table(t$accept == "Yes")["TRUE"]/nrow(t)

t |> ggplot(aes(x = sample_x, fill = accept)) + 
  geom_histogram(binwidth = 0.01, color = "white", position = "stack", boundary = 0) + 
  scale_fill_manual(values = c("Yes" = "#FF6F68", "No" = "#00CED1")) + 
  labs(
    title = "Example for generating Beta(6, 3) distribution", 
    x = "sample_x", 
    y = "count", 
    fill = "accept"
  )+ 
  coord_cartesian(xlim = c(0, 1))+
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

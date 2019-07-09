library(compositions)

size = 1.5
alpha = 0.5
shape = 20

col = "black"

auto_funnel <- function(number,ixx,s1,s2,mean_1,mean_2)
{
  std_1 <- sqrt(s1)
  std_2 <- sqrt(s2)
  
  p_value <- 0.6
  
  p_value <- log((p_value) * sqrt((exp(s1)-1)*(exp(s2)-1)) + 1) / (std_1 * std_2)
  
  cov_matrix = matrix(c(s1, std_1 * std_2 * p_value, std_1 * std_2 * p_value, s2), 2)
  
  Aeigen=eigen(cov_matrix)$values
  
  indexAeigen=1
  #print((Aeigen))
  #print(length(Aeigen))
  while (indexAeigen<=length(Aeigen))
  {
    if(Aeigen[indexAeigen]<0)
    {
      #print(Aeigen[indexAeigen])
      s1 <- as.double(runif(1,1,10))
      s2 <- as.double(runif(1,1,10))
      std_1 <- sqrt(s1)
      std_2 <- sqrt(s2)
      p_value <- 0.6
      p_value <- log((p_value) * sqrt((exp(s1)-1)*(exp(s2)-1)) + 1) / (std_1 * std_2)
      cov_matrix = matrix(c(s1, std_1 * std_2 * p_value, std_1 * std_2 * p_value, s2), 2)
      Aeigen=eigen(cov_matrix)$values
      #print(Aeigen[indexAeigen])
      #print((Aeigen))
      indexAeigen=1
    }
    else
    {
      #print(Aeigen[indexAeigen])
      indexAeigen=indexAeigen+1
    }
  }
  
  frame <- mvrnorm(number, mu = c(mean_1, mean_2), Sigma = cov_matrix)
  
  x <- frame[, 1]
  y <- frame[, 2]
  
  x <- exp(x)
  y <- exp(y)
  
  data_frame <- data.frame(x, y)
  
  p <- ggplot(data_frame, aes(x, y))
  
  p <- p + geom_point(size = size, shape = shape, col = col, alpha = alpha)
  p <- p + theme(axis.title = element_blank())
  p <- p + theme(axis.text.y = element_blank(), axis.text.x = element_blank(),axis.title = element_blank())
  
  fname <- "D:\\InfoVis\\Scag_Data_Builder\\funnel\\csv\\funnel_"
  numname <- as.character(number)
  fname <- paste(fname, numname, sep = '') 
  fname <- paste(fname, "_", sep = '')
  numname <- as.character(ixx)
  fname <- paste(fname, numname, sep = '')
  fname_csv <- paste(fname, ".csv", sep = '')
  fname <- "D:\\InfoVis\\Scag_Data_Builder\\funnel\\png\\funnel_"
  numname <- as.character(number)
  fname <- paste(fname, numname, sep = '') 
  fname <- paste(fname, "_", sep = '')
  numname <- as.character(ixx)
  fname <- paste(fname, numname, sep = '')
  fname_png <- paste(fname, ".png", sep = '')
  #print(fname)
  
  ggsave(fname_png,p,width = 3,height = 3)
  #print(p)
  write.table(data_frame,
              file = fname_csv,
              append = FALSE, 
              quote = FALSE, 
              sep = ",", 
              col.names = FALSE, 
              row.names = FALSE)
  
  #print(scagnostics(data_frame))
  
  print(p)
}

#for (i in 1:9)
{
  i=9
  number <- i*100
  for(ixx in 1:50)
  {
    s1 <- as.double(runif(1,1,10))
    s2 <- as.double(runif(1,1,10))
    mean_1 <- as.double(runif(1,1,10))
    mean_2 <- as.double(runif(1,1,10))
    auto_funnel(number,ixx,s1,s2,mean_1,mean_2)
  }
}


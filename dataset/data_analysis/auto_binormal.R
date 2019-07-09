library(MASS)
library(ggplot2)


auto_bino <- function(number, ixx=1)
{
  p_value = 0.6
  
  size = 1.5
  alpha = 0.5
  shape = 20
  col = "black"
  
  s1 = 2
  s2 = 2
  
  bivn <- mvrnorm(number, mu = c(0, 0), Sigma = matrix(c(s1^2, s1*s2*p_value, s2*s1*p_value, s2^2), 2))
  
  bivn.kde <- kde2d(bivn[,1], bivn[,2], n = number)
  
  data_frame <- data.frame(bivn)
  
  p <- ggplot(data_frame, aes(bivn[,1], bivn[,2]))
  p <- p + geom_point(size = size, shape = shape, col = col, alpha = alpha)
  p <- p + theme(axis.text.y = element_blank(), axis.text.x = element_blank(),axis.title = element_blank())
  
  fname <- "D:\\InfoVis\\Scag_Data_Builder\\bino\\bino_"
  numname <- as.character(number)
  fname <- paste(fname, numname, sep = '') 
  fname <- paste(fname, "_", sep = '')
  numname <- as.character(ixx)
  fname <- paste(fname, numname, sep = '')
  fname_csv <- paste(fname, ".csv", sep = '')
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
  
  print(p)
}


for(i in 1:9)
{
  
  num = i*100
  #num=150
  #print(num)

  for(ixx in 1:50)
  {
    auto_bino(num,ixx)
  }
  
}

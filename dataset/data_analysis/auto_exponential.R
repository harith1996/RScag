library(ggplot2)
library(scagnostics)


mean = 10
rate = 1 / mean
range = 1
size = 3
shape = 20
alpha = 0.7
#shape = 21

auto_exp <- function(number,ixx=1,rate=3,range=3,size=1.5,shape=20,alpha=0.5)
{
  col = "black"
  
  x <- rexp(n = number, rate)
  y <- rexp(n = number, rate)
  
  x <- x + runif(length(x), -range, range)
  y <- y + runif(length(y), -range, range)
  
  data_frame <- data.frame(x, y)
  
  p <- ggplot(data_frame, aes(x, y))
  
  p <- p + geom_point(size = size, shape = shape, col=col, alpha = alpha)
  p <- p + theme(axis.text.y = element_blank(), axis.text.x = element_blank(),axis.title = element_blank())
  
  fname <- "D:\\InfoVis\\Scagnostics18\\Scag_Data_Builder\\expo_\\expo_"
  numname <- as.character(number)
  fname <- paste(fname, numname, sep = '') 
  fname <- paste(fname, "_", sep = '')
  numname <- as.character(ixx)
  fname <- paste(fname, numname, sep = '')
  #fname <- paste(fname, "_", sep = '')
  #numname <- as.character(range)
  #fname <- paste(fname, numname, sep = '')
  fname_csv <- paste(fname, ".csv", sep = '')
  fname_scag <- paste(fname, "_scag.csv", sep = '')
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
  
  range_x <- max(x) - min(x)
  range_y <- max(y) - min(y)
  
  y <- y * range_x / range_y
  
  data_frame <- data.frame(x, y)
  
  result <- scagnostics(data_frame)
  
  result <- c(result)
  
  write.table(result,
              file = fname_scag,
              append = FALSE, 
              quote = FALSE, 
              sep = ",", 
              col.names = FALSE, 
              row.names = c("Outlying", "Skewed", "Clumpy", "Sparse", "Striated", "Convex", "Skinny", "Stringy", "Monotonic"))
  
  #print(p)
}


for(i in 1:9)
{
  
  num = i*100
  #num=150
  #print(num)
  
  for(ixx in 1:50)
  {
    range <- as.double(runif(1,1,3))
    #print(range)
    auto_exp(num,ixx,rate,range)
  }
  
}

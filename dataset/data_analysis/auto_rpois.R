library(ggplot2)

auto_rpois <- function(n,ixx=1, lambda = 5)
{
  out_x <- c()
  out_y <- c()
  
  for(i in seq(0, sqrt(n) - 1))
  {
    for(j in seq(0, sqrt(n) - 1))
    {
      number <- rpois(1, lambda = lambda)
      x <- i
      y <- j
      for(k in seq(0, number))
      {
         out_x <- c(out_x, x + runif(1, 0, 1))
         out_y <- c(out_y, y + runif(1, 0, 1))
      }
    }
  }
  
  data_frame <- data.frame(out_x, out_y)
  
  p <- ggplot(data_frame, aes(out_x, out_y))
  p <- p + geom_point(size = 1.5, shape = 20, col="black", alpha = 0.5)
  # p <- p + theme(axis.text = element_blank())
  p <- p + theme(axis.text.y = element_blank(), axis.text.x = element_blank(),axis.title = element_blank())
  
  fname <- "D:\\InfoVis\\Scag_Data_Builder\\rpois\\rpois_"
  numname <- as.character(n)
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
}

for(i in 1:9)
{
  
  num = i*15
  #num=150
  #print(num)
  
  for(ixx in 1:50)
  {
    auto_rpois(num,ixx, 8)
  }
  
}

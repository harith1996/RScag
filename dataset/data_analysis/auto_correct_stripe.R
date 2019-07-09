library(ggplot2)


auto_stripe <- function(number, ixx, yrange)
{
  
  size = 1.5
  alpha = 0.5
  shape = 20
  col = "black"
  
  min_x <- 1
  max_x <- 5
  
  min_y <- 1
  max_y <- min_y + yrange
  
  x <- runif(number, min_x, max_x + 1)
  y <- runif(number, min_y, max_y*as.double(runif(1,0.2,1.5)))
  
  x <- as.integer(x)
  
  data_frame <- data.frame(x, y)
  
  p <- ggplot(data_frame, aes(x, y))
  
  p <- p + geom_point(size = size, shape = shape, col = col, alpha = alpha)
  p <- p + theme(axis.text.y = element_blank(), axis.text.x = element_blank(),axis.title = element_blank())
  
  fname <- "D:\\InfoVis\\Scag_Data_Builder\\stripe\\stripe_"
  numname <- as.character(number)
  fname <- paste(fname, numname, sep = '') 
  fname <- paste(fname, "_", sep = '')
  numname <- as.character(ixx)
  fname <- paste(fname, numname, sep = '')
  fname <- paste(fname, "_", sep = '')
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
    yrange <- as.integer(runif(1,10,100))
    auto_stripe(num,ixx,yrange)
  }
  
}

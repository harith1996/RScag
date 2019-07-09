library(ggplot2)

number = 1000 
a = 1
b = 1
c = 1

auto_quad <- function(number,ixx=1,a=1,b=1,c=1)
{
  size = 1.5
  shape = 20
  alpha = 0.5
  liftted = (abs(a)+abs(b)+abs(c))/3
  range1 <- as.double(runif(1,3,10))
  range2 <- as.double(runif(1,3,10))
  step = (range1+range2)/number
  
  cex=-b/(2*a)
  
  col = "black"
  
  x <- seq(-range1+cex, range2+cex, step)
  
  y <- a * x ** 2 + b * x ** 1 + c * x ** 0
  
  y <- y + runif(length(y), -liftted, liftted)
  
  data_frame <- data.frame(x, y)
  
  p <- ggplot(data_frame, aes(x, y))
  
  p <- p + geom_point(size = size, shape = shape, col=col, alpha = alpha)
  p <- p + theme(axis.text.y = element_blank(), axis.text.x = element_blank(),axis.title = element_blank())
  
  fname <- "D:\\InfoVis\\Scag_Data_Builder\\quad\\quad_"
  numname <- as.character(number)
  fname <- paste(fname, numname, sep = '') 
  fname <- paste(fname, "_", sep = '')
  numname <- as.character(ixx)
  fname <- paste(fname, numname, sep = '')
  fname <- paste(fname, "_", sep = '')
  numname <- as.character(as.integer(a))
  fname <- paste(fname, numname, sep = '')
  fname <- paste(fname, "_", sep = '')
  numname <- as.character(as.integer(b))
  fname <- paste(fname, numname, sep = '')
  fname <- paste(fname, "_", sep = '')
  numname <- as.character(as.integer(c))
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
    a <- as.double(runif(1,-10,10))
    b <- as.double(runif(1,-10,10))
    c <- as.double(runif(1,-10,10))
    auto_quad(num,ixx,a,b,c)
  }
  
}
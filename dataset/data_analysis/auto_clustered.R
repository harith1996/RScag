library(ggplot2)

triangle_length = 20
number = 1000

auto_clus <- function(number,ixx=1,triangle_length1 = 20,triangle_length2 = 20,triangle_length3 = 20, inner_radio, outer_radio)
{
  
  X1 = -triangle_length1 / 2
  Y1 = 0
  Z1 = 0
  
  X2 = triangle_length2 / 2
  Y2 = 0
  Z2 = 0
  
  X3 = 0
  Y3 = triangle_length3 / 2 * sqrt(3)
  Z3 = 0
  
  size = 1.5
  alpha = 0.5
  shape = 20
  
  col = "black"
  
  inner_radio = 3
  outer_radio = 3
  
  r1 = (runif(number) * (outer_radio ^ 3 - inner_radio ^ 3) + inner_radio ^ 3) * 1/3
  phi1 = acos(runif(number, -1, 1));
  th1 = 2 * pi * runif(number)
  
  x1 = r1 * sin(phi1) * sin(th1) + X1
  y1 = r1 * sin(phi1) * cos(th1) + Y1
  z1 = r1 * cos(phi1) + Z1
  
  r2 = (runif(number) * (outer_radio ^ 3 - inner_radio ^ 3) + inner_radio ^ 3) * 1/3
  phi2 = acos(runif(number, -1, 1));
  th2 = 2 * pi * runif(number)
  
  x2 = r2 * sin(phi2) * sin(th2) + X2
  y2 = r2 * sin(phi2) * cos(th2) + Y2
  z2 = r2 * cos(phi2) + Z2
  
  r3 = (runif(number) * (outer_radio ^ 3 - inner_radio ^ 3) + inner_radio ^ 3) * 1/3
  phi3 = acos(runif(number, -1, 1));
  th3 = 2 * pi * runif(number)
  
  x3 = r3 * sin(phi3) * sin(th3) + X3
  y3 = r3 * sin(phi3) * cos(th3) + Y3
  z3 = r3 * cos(phi3) + Z3
  
  x <- c(x1, x2, x3)
  y <- c(y1, y2, y3)
  z <- c(z1, z2, z3)
  
  
  data_frame <- data.frame(x, y)
  
  p <- ggplot(data_frame, aes(x, y))
  p <- p + geom_point(size = size, shape = shape, col = col, alpha = alpha)
  p <- p + theme(axis.text.y = element_blank(), axis.text.x = element_blank(),axis.title = element_blank())
  
  fname <- "D:\\InfoVis\\Scag_Data_Builder\\clus\\csv\\clus_"
  numname <- as.character(number*3)
  fname <- paste(fname, numname, sep = '') 
  fname <- paste(fname, "_", sep = '')
  numname <- as.character(ixx)
  fname <- paste(fname, numname, sep = '')
  fname_csv <- paste(fname, ".csv", sep = '')
  fname <- "D:\\InfoVis\\Scag_Data_Builder\\clus\\png\\clus_"
  numname <- as.character(number*3)
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
  
  print(p)
}

for(i in 1:9)
{
  
  num = i*100
  #num=150
  #print(num)
  
  for(ixx in 1:50)
  {
    outer_radio <- as.double(runif(1,0.01,1))
    inner_radio <- as.double(runif(1,0.01,1))
    triangle_length1 <- as.integer(runif(1,10,100))
    triangle_length2 <- as.integer(runif(1,10,100))
    triangle_length3 <- as.integer(runif(1,10,100))
    auto_clus(num,ixx,triangle_length1,triangle_length2,triangle_length3,inner_radio,outer_radio)
  }
  
}

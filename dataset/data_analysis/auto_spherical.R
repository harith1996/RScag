library(ggplot2)

auto_sphe <- function(number, ixx = 1, inner_radio = 3, outer_radio = 3, triangle_length = 20)
{
  X1 = -triangle_length / 2
  Y1 = 0
  Z1 = 0
  
  size = 1.5
  alpha = 0.5
  shape = 20
  
  col = "black"
  
  r1 = (runif(number) * (outer_radio ^ 3 - inner_radio ^ 3) + inner_radio ^ 3) * 1/3
  phi1 = acos(runif(number, -1, 1));
  th1 = 2 * pi * runif(number)
  
  x1 = r1 * sin(phi1) * sin(th1) + X1
  y1 = r1 * sin(phi1) * cos(th1) + Y1
  z1 = r1 * cos(phi1) + Z1
  
  x <- x1
  y <- y1
  z <- z1
  
  data_frame <- data.frame(x, y)
  
  p <- ggplot(data_frame, aes(x, y))
  p <- p + geom_point(size = size, shape = shape, col = col, alpha = alpha)
  p <- p + theme(axis.text.y = element_blank(), axis.text.x = element_blank(),axis.title = element_blank())
  
  fname <- "D:\\InfoVis\\Scag_Data_Builder\\sphe\\sphe_"
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

  num = i*250
  #num=150
  #print(num)
  
  for(ixx in 1:50)
  {
    inner_radio <- as.double(runif(1,0.01,1))
    outer_radio <- as.double(runif(1,0.01,1))
    auto_doug(num,ixx,inner_radio,outer_radio)
  }
  
}


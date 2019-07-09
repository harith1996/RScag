library(ggplot2)

auto_doug <- function(number,ixx=1,inner_radio = 3, outer_radio = 6, width=1)
{
  
  X = 0
  Y = 0
  
  size = 1.5
  alpha = 0.5
  shape = 20
  col = "black"
  
  
  get_outer_data <- function(number, inner_radio, outer_radio, width)
  {
    return_x <- c()
    return_y <- c()
    
    while(TRUE)
    {
      if(length(return_x) == number)
      {
        break
      }
      else
      {
        x <- runif(1, -outer_radio, outer_radio)
        y <- runif(1, -outer_radio, outer_radio)
        
        value <- x * x + y * y
        
        if(value > inner_radio ** 2 && value < outer_radio ** 2)
        {
          return_x <- c(return_x, x)
          return_y <- c(return_y, y)
        }
      }
    }
    return(data.frame(return_x, return_y))   
  }
  
  data_frame <- get_outer_data(number, inner_radio, outer_radio, width)
  
  p <- ggplot(data_frame, aes(data_frame$return_x, data_frame$return_y))
  
  p <- p + geom_point(size = size, shape = shape, col = col, alpha = alpha)
  p <- p + theme(axis.text.y = element_blank(), axis.text.x = element_blank(),axis.title = element_blank())
  
  fname <- "D:\\InfoVis\\Scag_Data_Builder\\doug\\csv\\doug_"
  numname <- as.character(number)
  fname <- paste(fname, numname, sep = '') 
  fname <- paste(fname, "_", sep = '')
  numname <- as.character(ixx)
  fname <- paste(fname, "_", sep = '')
  numname <- as.character(inner_radio)
  fname <- paste(fname, "_", sep = '')
  numname <- as.character(outer_radio)
  fname <- paste(fname, "_", sep = '')
  numname <- as.character(width)
  fname <- paste(fname, numname, sep = '')
  fname_csv <- paste(fname, ".csv", sep = '')
  fname <- "D:\\InfoVis\\Scag_Data_Builder\\doug\\png\\doug_"
  numname <- as.character(number)
  fname <- paste(fname, numname, sep = '') 
  fname <- paste(fname, "_", sep = '')
  numname <- as.character(ixx)
  fname <- paste(fname, "_", sep = '')
  numname <- as.character(inner_radio)
  fname <- paste(fname, "_", sep = '')
  numname <- as.character(outer_radio)
  fname <- paste(fname, "_", sep = '')
  numname <- as.character(width)
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
  
  num = i*250
  #num=150
  #print(num)
  
  for(ixx in 1:50)
  {
    inner_radio <- as.double(runif(1,3,20))
    outer_radio <- as.double(runif(1,9,25))
    if (inner_radio>outer_radio)
    {
      tmp=inner_radio
      inner_radio=outer_radio
      outer_radio=tmp
    }
    width <- as.double(runif(1,5,15))
    auto_doug(num,ixx,inner_radio,outer_radio,width)
  }
  
}


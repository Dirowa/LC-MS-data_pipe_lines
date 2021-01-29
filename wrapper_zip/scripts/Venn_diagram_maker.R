
library("eulerr")


file <- '!@#$%^&1&^%$#@!'
output <- '!@#$%^&2&^%$#@!'

pixel_size <- !@#$%^&3&^%$#@!


  

  data_frame <- read.csv(file)
  data_frame <- data_frame[,-1]
  
  file1 <- substr(file, 1, nchar(file)-4)
  
  file_name = paste0(file1,'.png')
  print(file_name)
  png(filename = file_name,
      width = 1840, height = 1840, units = "px", pointsize = pixel_size,
      bg = "white",  res = NA,
  )
  
  venn <- plot(euler(data_frame),
               #fills = TRUE,
               #edges = TRUE,
               legend = TRUE,
               labels = identical(legend, TRUE),
               quantities = TRUE,
               strips = NULL,
               main = NULL,
               n = 200L,
               adjust_labels = TRUE,
  )
  print(venn)
  
  dev.off()
  



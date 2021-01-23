install.packages('eulerr')

library("eulerr")


path <- 'F:/avans/stage MM/peakpicking/'
output <- 'F:/avans/stage MM/peakpicking/'
files <- list.files(path, pattern = 'venn_diagram_data.csv')




for (file in files)  {

  print(file)
  path1 = paste0(path,file)
  print(path1)
  data_frame <- read.csv(path1)
  data_frame <- data_frame[,-1]
  
  file1 <- substr(file, 1, nchar(file)-4)
  
  file_name = paste0(output,file1,'.png')
  print(file_name)
  png(filename = file_name,
      width = 1840, height = 1840, units = "px", pointsize = 12,
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
  
  
}



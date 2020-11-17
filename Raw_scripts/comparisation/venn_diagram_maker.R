
install.packages('eulerr')
library("eulerr")

items <- c(1,3,4,7,8)


  
  
path <- paste0("C:/Users/DVrin/OneDrive/Documenten/avans/stage_MM/peakpicking comparisation/_venn_diagram_data_",8,".csv")
data_frame <- read.csv(path)
data_frame <- data_frame[,-1]
  
  
plot(euler(data_frame),
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



path <- "F:/avans/stage MM/peakpicking_output/venn_diagram_data.csv"
data_frame <- read.csv(path)
data_frame <- data_frame[,-1]


plot(euler(data_frame), counts = T, fontface = 1)




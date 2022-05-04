library(readxl)
# set as working directory the folder with the multiple *.csv files
setwd("C:/Users/alex/Desktop/SCRINSHOT_updated_tutorial/analysis_20220503/output/roi_counts/")


# import the list with all cells, using the full path. It should NOT be the working directory, with the gene *.csv files
all_cell_roi_list <- read.csv("C:/Users/alex/Desktop/SCRINSHOT_updated_tutorial/analysis_20220503/all_cell_roi_list.csv", row.names=1)

#================================================================================
# just press "run" in line-15 and the script will produce the 2 *.csv files with all the dot counts (integer and not)
{
# obtain a list of all *.csv files in the folder
list = list.files(pattern="*.csv")

# create an empty data.frame, with proper dimensions
all_dots <- data.frame(matrix(nrow = nrow(all_cell_roi_list), ncol = length(list)+2))


# import the ROI names in the first column
all_dots[,1]<- all_cell_roi_list$Name
# rename the column
colnames(all_dots)[1] <- "roi"

# open the first file of the list
count <-read.csv(list[1], row.names=1)
# obtain the column with the roi area (pixel^2) and add it to the data.frame
all_dots[,2]<- count$Area
# rename the column
colnames(all_dots)[2] <- "Area"
# obtain the mean fluorescence intensity of the 1st gene and convert it to number of pixels
all_dots[,3] <- count$Mean*count$Area/255
# rename based on the name of the imported *.csv file
n <-as.character(list[1])
s<-gsub("_.*","",n)
colnames(all_dots)[3] <- s
rm(count, n, s)

# do the same procedure with the rest of the genes
for (i in 2: length(list)){
  count <- read.csv(list[i], row.names=1)
  all_dots[,(i+2)] <- count$Mean*count$Area/255
  n <-as.character(list[i])
  s<-gsub("_.*","",n)
  colnames(all_dots)[i+2] <- s
  rm(count, n, s)
}

# export a *.csv file with the merged dot counts. 
write.csv(all_dots, file="all_dots.csv")

# transform the values to become integer and export the result as *.csv
all_dots1 <- all_dots
for (i in 3:ncol(all_dots1)){
  all_dots1[,i] <- round(as.numeric(all_dots1[,i]))
}
write.csv(all_dots1, file="all_dots_integer.csv")
}
# merged files have been saved in the working directory, by default.

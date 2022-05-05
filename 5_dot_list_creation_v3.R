library(plyr)
library(dplyr)
library(data.table)
library(readxl)
#set the number and size of tiles per dimension
x_n=2
y_n=2
x_size=3350
y_size=2500

in_wd <- "C:/Users/alex/Desktop/SCRINSHOT_updated_tutorial/analysis_20220503/input/"
out_wd <- "C:/Users/alex/Desktop/SCRINSHOT_updated_tutorial/analysis_20220503/output/"

channel_order <- read_excel(paste0(in_wd, "channel_order.xlsx"))
gene_order <- channel_order[2:nrow(channel_order),c(4:5)]
# set the number detected genes
gene_n=19

# run the following lines to create the result files
{
mydir = out_wd
myfiles = list.files(path=mydir, pattern="MyExpt_IdentifiedBlobs*", full.names=TRUE)
myfiles

for (i in 1:gene_n){
  g<- data.frame(read.csv(paste0(mydir, "MyExpt_IdentifiedBlobs",i, ".csv")))
  if (nrow(g)>0){
  g$Name <- paste0("gene",i)
  assign(paste0("gene",i),g)
  rm(g)
  } else {}
}

g1<- data.frame(read.csv(paste0(mydir, "MyExpt_IdentifiedBlobs",1, ".csv")))
g1 <- g1[,c("ImageNumber", "Location_Center_X","Location_Center_Y")]
g1$Name <- paste0("gene",1)

for (i in 2:gene_n){
  g<- data.frame(read.csv(paste0(mydir, "MyExpt_IdentifiedBlobs",i, ".csv")))
  if (nrow(g)>0){
  g <- g[,c("ImageNumber","Location_Center_X","Location_Center_Y")]
  g$Name <- paste0("gene",i)
  g1 <- rbind(g1, g)
  rm(g)
  } else {}
}

tr_matrix <- data.frame(matrix(data = 0, nrow = x_n*y_n,ncol = 3))
colnames(tr_matrix) <- c("ImageNumber","X","Y")

n_r <- x_n*y_n
for (i in 1:n_r){
  tr_matrix[i,1] <-i
}

for (i in 1:x_n){
  for (j in (0:(y_n-1))){
    for (k in (i + j*x_n)){
      tr_matrix[k, 2] <- (i*x_size)-x_size
    }
  }
}

r <- c(1:n_r)
splitted <- split(r,             
                  cut(seq_along(r),
                      y_n,
                      labels = FALSE))

for (i in 1:y_n){
  for (j in splitted[[i]])
    tr_matrix[j, 3] <- (i*y_size)-y_size
}

g2 <- join(g1, tr_matrix, by="ImageNumber")
g2$Location_Center_X <- g2$Location_Center_X+g2$X
g2$Location_Center_Y <- g2$Location_Center_Y+g2$Y
g2 <- g2[,1:4]

names(g2)[4] <- "Gene_Order"

g2 <- join(g2, gene_order, by="Gene_Order")
write.csv(g2, paste0(mydir,"all_gene_marker_coordinates.csv"))

genes <- levels(as.factor(g2$Gene))
dir.create(paste0(out_wd, "single_gene_dot_coordinates"))
for (i in genes){
  a <- g2[g2$Gene==i, ]
  write.csv(a, file=paste0(out_wd, "single_gene_dot_coordinates", "/", i, "_dot_coordintates.csv"))
  rm(a)
}

}






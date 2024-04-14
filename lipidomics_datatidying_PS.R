
#In vitro data tidying
#Read original excel file for 253 cell line
library(readxl)
data_253 <- read_excel("Copia de in vitro_mzml-out_5_7_1000-1.xlsx",sheet = "253")
#Delete every column but 1,2 and 22:30
data_253 <- data_253[,c(1,22:30)]
#Delete every row with NAs
data_253 <- data_253[complete.cases(data_253),]

#Repeat everything but with sheet="215"
data_215 <- read_excel("Copia de in vitro_mzml-out_5_7_1000-1.xlsx",sheet = "215")
data_215 <- data_215[,c(1,22:30)]
data_215 <- data_215[complete.cases(data_215),]

#Reorder data_253 to match data_215[[1]]
data_253 <- data_253[match(data_215[[1]],data_253[[1]]),]
data_253[,1]<-data_215[,1]

#Create a new merged dataframe from the two
data_253_215 <- merge(data_253,data_215,by="NAME",all=TRUE)


#Prepare sample names data_253_215 column names
sample_names<-names(data_253_215[,-1])

#Simplify data samples to Sx
names(data_253_215)<-c("NAME",paste0("S",1:18))

#Make Clinical (phenotypic excel)
data_clin_invitro<-data.frame(sample=names(data_253_215[,-1]), cell_line=as.factor(rep(c("line_253","line_215"),each=9)), 
                      treatment=rep(c("control", "eto", "mcm"), each=3, times=2),replica=as.factor(rep(c("rep_1","rep_2","rep_3"),times=6)),stringsAsFactors = FALSE)


#Eliminate parenthesis from data_253_215[,1]
data_matrix_invitro<-data_253_215
data_matrix_invitro[,1]<-gsub("\\(|\\)","",data_253_215[,1])
#Add a space between the last letter and first number of data_matrix[,1]
data_matrix_invitro[,1]<-gsub("(?<=[A-Za-z])(?=[0-9])"," ",data_matrix_invitro[,1],perl=TRUE)

#Delete all "-" in "-O" in data_matrix[,1]
data_matrix_invitro[,1]<-gsub("-O","O",data_matrix_invitro[,1])

#Change - symbol in some fatty acid names to make it standard
data_matrix_invitro[79,1]<-"DiHexCerd 18:1/14:13"
data_matrix_invitro[82,1]<-"DiHexCerd 18:1/16:15"
data_matrix_invitro[87,1]<-"DiHexCerd 18:1/22:21"
data_matrix_invitro[89,1]<-"DiHexCerd 18:1/24:23"
data_matrix_invitro[93,1]<-"DiHexCerd 20:1/12:11"



#Show if there are duplicates
cat(sum(duplicated(data_matrix_invitro[[1]])))

#Subset both cell lines
data_matrix_253<-data_matrix_invitro[, 1:10]
#Eliminate all rows with NAs in data_matrix_253
data_matrix_253<-data_matrix_253[complete.cases(data_matrix_253),]

data_matrix_215<-data_matrix_invitro[, c(1,11:19)]
#Eliminate all rows with NAs in data_matrix_253
data_matrix_215<-data_matrix_215[complete.cases(data_matrix_215),]

#Export data_clin and data_253_215 to .csv, without rownames
write.csv(data_clin_invitro,"data_clin_invitro.csv", row.names = FALSE)
write.csv(data_matrix_invitro,"data_matrix_invitro.csv", row.names = FALSE)



#Import datasets into a LipidomicsExperiment object w/ lipidr
library(lipidr)

d_invitro <- as_lipidomics_experiment(read.csv("data_matrix_invitro.csv"))
d_invitro <- add_sample_annotation(d_invitro, "data_clin_invitro.csv")





#Load unclean library
# Find and replace in excel to getrid "/"
# Find: Pristiformes/Rhiniformes group
# Replace: Rhinopristiformes

# Find: Rhynchobatus cf. laevis KKB-2014
# Replace: Rhynchobatus laevis

# working directory
setwd("sharkdust")

data <- read.table("Data/0_raw data/2_Dust_d1_MOTUs_db2020.csv",sep=",",head=T)
colnames(data)

#subset of data containing just Sebastes
pc_1 <- subset(data,scientific_name=="Sebastes")

#find row with highest PC reads
pc_1 <- pc_1[1,]
y <- sum(pc_1[,28:66])
z <- sum(pc_1[,67])
pc <- y/z

# Select different groups of columns

data_taxo <- data[,1:15] 
data_samples <- data[,28:64]
data_blank <- data[,65:66] 
data_PC <- data[,67] 
data_seqs <- data[,91:92]


# Select amount of reads to remove based on tagjumping percentage (using positive control samples)

calculate_reads <- (data_samples*pc) #calculate the amount of reads to remove for each MOTU based on the previously defined percentage 
calculate_reads <- round(calculate_reads)
data_samples <- data_samples-calculate_reads

# Reconstruct the matrix

data_seqs$total_reads <- rowSums(data_samples,na.rm = T)
Lib_PC <- data.frame(data_taxo,data_blank,data_samples,data_seqs)
sum(Lib_PC$total_reads)

# Write the new output

write.table(Lib_PC,"Data/1_cleaning/1_Dust_PC_Corrected.csv",row.names=F,sep=";")

Lib_PC <- subset(Lib_PC, family_name!="Hominidae") # removing human reads 

sum(Lib_PC$total_reads)

data <- Lib_PC
colnames(data)

# Select different groups of columns

data_taxo <- data[,1:15] 
data_blank <- data[,16:17] 
data_samples <- data[,18:54]
data_seqs <- data[,55:56]


# Select blank reads to remove
data_blank$max <- do.call(pmax, data_blank)
blank_remove <- data_blank$max


# remove blank reads from samples for every row
for (i in 1:nrow(data_samples)){
  data_samples[i,] <- data_samples[i,] - data_blank$max[i]
  data_samples[i,][data_samples[i,]<2] <- 0
}

# Reconstruct the matrix
data_seqs$total_reads <- rowSums(data_samples,na.rm = T)

Lib_blanks <- data.frame(data_taxo,data_samples,data_seqs)

# Save the new matrix

write.table(Lib_blanks,"Data/1_cleaning/2_Dust_Blanks_Corrected.csv",row.names=F,sep=";")

#Additional filtering steps to remove reads belonging to humans and domestic animals
sum(Lib_blanks$total_reads) #obtain the total of reads after removing reads from the blanks 
Lib_sharks <- subset(Lib_blanks, class_name=="Chondrichthyes") #select only mammals
sum(Lib_sharks$total_reads)
write.csv(Lib_sharks, "Data/1_cleaning/3_Dust_Sharks_MOTUs.csv")

Lib_Sharks_nonhuman <- subset(Lib_sharks, family_name!="Hominidae") # removing human reads 
sum(Lib_Sharks_nonhuman$total_reads)
write.csv(Lib_Sharks_nonhuman, "Data/1_cleaning/4_Dust_Human_Corrected.csv")


#ewt1_mammals_wild <- subset(ewt1_mammals_nonhuman, genus_name!= "Ovis")# removing domestic animals 
# ewt1_mammals_wild1 <- subset(ewt1_mammals_wild, genus_name!= "Bos")
# ewt1_mammals_wild2 <- subset(ewt1_mammals_wild1, genus_name!= "Sus")
# ewt1_mammals_wild3 <- subset(ewt1_mammals_wild2, genus_name!= "Equus")
# ewt1_mammals_wild4 <- subset(ewt1_mammals_wild3, genus_name!= "Canis")
# ewt1_mammals_wild5 <- subset(ewt1_mammals_wild4, genus_name!= "Felis")
# sum(ewt1_mammals_wild5$total_reads)
# write.csv(ewt1_mammals_wild5, "EWT1_Wild_Mammals.csv")



#filter out MOTUs with a best identity <98% and fewer than 5 reads in total
Lib_Sharks_fin <- subset(Lib_Sharks_nonhuman, best_identity >= 0.95)
Lib_Sharks_final <- subset(Lib_Sharks_fin, total_reads > 5)
sum(Lib_Sharks_final$total_reads)

#final filtered csv file
write.table(Lib_Sharks_final,"Data/1_cleaning/5_Dust_Final.csv",row.names=F,sep=";")

# RUN IN LINUX TERMINAL
# Collapse similar MOTUs (in R_scripts_metabarpark-master)
Rscript owi_collapse -i Data/1_cleaning/5_Dust_Final.csv -t 0.70 -s 28 -e 67

## Output:
# 5_Dust_Final_collapsed.csv
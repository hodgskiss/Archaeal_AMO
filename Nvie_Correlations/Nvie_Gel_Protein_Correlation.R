#Set working directory to location with script and required files.
#Required Files:
  #Nvie_22_Protein_Groups.txt
  #Nvie_Gene_Codes.csv

library(tidyverse)
library(expss)

#Read in protein group data.
proteins <- read.table("Nvie_22_Protein_Groups.txt",sep = "\t",dec=".",header=TRUE)

#Remove REV and CON identifications.
proteins <- proteins[33:650,]

#Renumber rows in protein datafrane.
rownames(proteins) <- NULL


#Look for proteins that did not have a successful iBAQ calculation.
no_iBAQ <- proteins[is.nan(proteins$iBAQ),]

#From no_iBAQ we can see that NVIE_022250 and NVIE_004140 did not susccessfully compute iBAQ values.
#This is due to their short sequences and lack of trypsin digests that meet iBAQ requirements.
#Their intensity values will be used instead (i.e, Intensity divided by 1, see below).

#Extract Intensity data from NVIE_022250.
Nvie_022250_Int <- proteins[371,c(6,107:129)]
#Rearrange Intensity data to the correct order.
Nvie_022250_Int <- Nvie_022250_Int[,c(1:3,14,18:24,4:13,15:17)]

#Extract Intensity data from NVIE_004140.
Nvie_004140_Int <- proteins[376,c(6,107:129)]
#Rearrange Intensity data to the correct order.
Nvie_004140_Int <- Nvie_004140_Int[,c(1:3,14,18:24,4:13,15:17)]

#Pull the fasta headers and iBAQ columns from the dataframe.
proteins <- as.data.frame(proteins[,c(6,130:152)])

#Rearrange the iBAQ Intesity columns to the right order.
proteins <- proteins[,c(1:3,14,18:24,4:13,15:17)]

#Replace missing iBAQ Intensity values with normal Intensity values.
proteins[371,] <- Nvie_022250_Int
proteins[376,] <- Nvie_004140_Int
#Since these are small proteins, this would be the same as taking the number of correct
#trypsin sites to be 1.  This serves as an estimate in the absence of true iBAQ values.

#Shorten fasta headers to just the protein accession number.
proteins$Fasta.headers <- substr(proteins$Fasta.headers,4,13)


#Arrange the proteins from largest to smallest total iBAQ value.
#Also only takes the top 50% most abundant proteins; the interest is only in highly abundant proteins.

proteins <- proteins%>%arrange(desc(iBAQ))%>%slice_max(iBAQ,prop=0.5)

#Remove the total iBAQ column.
proteins$iBAQ <- NULL

#Store Accession numbers.
accession <- proteins$Fasta.headers

#Replace row names with accession numbers.
row.names(proteins) <- accession

#Remove fasta headers column.
proteins$Fasta.headers <- NULL

#Transform dataframe. 
t_proteins <- as.data.frame(t(proteins))

#Store names of each column (i.e., accession numbers). Order is important here.
names <- colnames(t_proteins)

#Count the numver of proteins in the dataframe t_protein.
p <- length(colnames(t_proteins))

#Store the value for the number of columns minus 1.
q <-p-1


#Create an empty vector to store data in. 
corr_list <- vector()

#Run a loop to perform a Kendall correlation analysis for each protein against every other protein. 
for(i in c(1:q)){
  
  #Retrive the name of the first protein.
  first <- as.character(names[i])
  
  #Use k to determine how many proteins to be compared to (i.e, for first protein: all proeins from 2 to 309,
  # for second protein: all proteins from 3 to 309, etc.).
 k <-i+1
  for(j in c(k:p)){
    
    #Retrieve the name of the protein to compare to.
    second <- as.character(names[j])
    
    #Runs a Kendall correlation analysis.
    cor <- cor.test(t_proteins[,first],t_proteins[,second],method="kendall",use="complete.obs")
    
    #Stores the protein names, correlation value, and p-value in list.
    list <- c(first,second,cor$estimate,cor$p.value)
    
    #Adds the list information to corr_list. This way, corr_list will collect all correlation information.
    corr_list <- rbind(corr_list, list)
    
    
  }
  
  
  
  
}

#There will be warnings after this loops.
#These are saying that p-values can't be computed with ties.
#In these cases, an approximate p-value is calculated.


#Convert corr_list to a dataframe. 
corr_list <- as.data.frame(corr_list)

#Store column V4 (the p-values) as numbers.
corr_list$V4 <- as.numeric(as.character(corr_list$V4))

#Calculate adjusted p-values using the Benjamini-Hochberg method.
corr_list<- corr_list%>%mutate(p_adj=p.adjust(V4,"BH"))

#Arrange dataframe by adjusted p-value, smallest to largest.
corr_list <- corr_list%>%arrange(p_adj)

#Rename columns in corr_list.
colnames(corr_list) <- c("Protein_1","Protein_2","tau","p_value","p_adj")

#Optional: run the marked out command below to write the full correlation list to a file. 
#write.table(corr_list, file="corr_list_NVB_gel", row.names = FALSE)

#Load a file with information about each protein and gene in the genome.
labels <- read.csv(file="Nvie_Gene_Codes.csv", header=TRUE, sep=",")

#Add columns for the gene information for each protein. 
corr_list <- mutate(corr_list, Gene_1=vlookup(corr_list$Protein_1,labels,lookup_column = "Protein.Code",result_column = "Gene.Name"))
corr_list <- mutate(corr_list, Gene_2=vlookup(corr_list$Protein_2,labels,lookup_column = "Protein.Code",result_column = "Gene.Name"))

#Add columns for the descriptive protein product of each protein. 
corr_list <- mutate(corr_list, Product_1=vlookup(corr_list$Protein_1,labels,lookup_column = "Protein.Code",result_column = "Protein.Product"))
corr_list <- mutate(corr_list, Product_2=vlookup(corr_list$Protein_2,labels,lookup_column = "Protein.Code",result_column = "Protein.Product"))

#Add columns for the gene locus tag of each protein. 
corr_list <- mutate(corr_list, Locus_1=vlookup(corr_list$Protein_1,labels,lookup_column = "Protein.Code",result_column = "Locus.Tag"))
corr_list <- mutate(corr_list, Locus_2=vlookup(corr_list$Protein_2,labels,lookup_column = "Protein.Code",result_column = "Locus.Tag"))


#Extract only correlations that involve known amo genes for the first or second protein. 
amo_corr <- with(corr_list, corr_list[grepl("amo",Gene_1)|grepl("amo",Gene_2), ])

#Optional: run the marked out command below to write the full amo correlation list to a file. 
#write.table(amo_corr, file="amo_corr", row.names = FALSE)

#Store column tau as numbers.
amo_corr$tau <- as.numeric(as.character(amo_corr$tau))

#Extract amo correlations that have a tau correlation value higher than or equal to 0.7 and with a p_adj less than or equl to 0.001.
amo_filtered <- amo_corr%>%filter(tau>=0.7&p_adj<=0.001)
#Note:In this case, all correlations are smaller than 0.001, so a more stringent cut-off was not deemed necessary.

#Remove any correlations involving the genes amoC1 or amoC3.
#Based on intensity and transcript data, these are not the primary AmoC protein groups.
#Thus, they were excluded.
#Exclude for Gene_1 and save in new file.
amo_focused <-subset(amo_filtered,!(Gene_1 %in% c("amoC1","amoC3")))
#Exclude for Gene_2 in already created file. 
amo_focused <-subset(amo_focused,!(Gene_2 %in% c("amoC1","amoC3")))


#Remove any correlations involving the gene NVIE_022250.
#NVIE_022250 is a duplication of amoX.
#Based on intensity and transcript data, this is not the primary Amo protein.
#Thus, it was excluded.
#Exclude for Gene_1 and save in already created file.
amo_focused <-subset(amo_focused,!(Gene_1 %in% c("NVIE_022250")))
#Exclude for Gene_2 in already created file. 
amo_focused <-subset(amo_focused,!(Gene_2 %in% c("NVIE_022250")))

amo_focused <- amo_focused%>%arrange(desc(tau))


#Optional: run the marked out command below to write the focused amo correlation list to a file. 
#write.table(amo_focused, file="amo_focused_NVB_gel", row.names = FALSE)

#Extract all genes that correlate with amo with respect to the above conditions.
Genes <- c(amo_focused$Gene_1, amo_focused$Gene_2)

#Extract all genes without duplicate values.
Genes_Simplified <- Genes[!duplicated(Genes)]

#Remove known AMO subunit genes.
Candidate_Gene_List <- Genes_Simplified[Genes_Simplified !="amoA"]%>%.[. !="amoB"]%>%.[. !="amoC4"]
Candidate_Gene_List

#Optional: run the marked out command below to write the candidate gene list to a file. 
#write.table(Candidate_Gene_List, file="NVIE_Candidate_AMO_Genes", row.names = FALSE, col.names = FALSE, quote = FALSE)

#This final list corresponds to the list in Table 1A of the paper.
#Note: AmoC4 and AmoC6 are represented in the same protein group. It is believed
#that this is actually AmoC6.  See paper for explanation.
#NOTE: It can be seen in amo_filtered that AmoY (NVIE_004540) correlates with
#amoA, AmoB, AmoC4, and AmoX. AmoZ (NVIE_004550) is missing from this analysis
#for reasons explained in the paper.

#Locus Tags for AMO and candidate genes can be extracted in a similar fashion.
Locus_Tags <- c(amo_focused$Locus_1, amo_focused$Locus_2)
Locus_Tags_Simplified <- Locus_Tags[!duplicated(Locus_Tags)]
Locus_Tags_Simplified



#><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
#  MOTU Sequence Data Processing 
#><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

library(data.table) # read and manipulate large data
library(dplyr) # data manipulation
library(vegan) # diversity estimates
library(seqRFLP) # create fasta file
library(EcolUtils) # Permutational Rarefraction

MOTU<-fread("Raw_OTU_Table.csv")

# creating dataframe of ID tied to sequence for MASCSE step 
SequenceCheck<-MOTU[,c("ID","sequ")]

#removing columns not needed
MOTU$sort<-NULL
MOTU$sequ<-NULL

# Abundance filtering to reduce the number of false positives due to PCR and sequencing errors 
MOTU_samples<-as.data.frame(MOTU[,2:ncol(MOTU)])
sum(MOTU_samples)
rownames(MOTU_samples)<-MOTU$ID
total_MOTU<-colSums(MOTU_samples)
MOTU_rel <- as.data.frame(MOTU_samples)
for (i in 1:ncol(MOTU_samples))  MOTU_rel[,i] <- MOTU_samples[,i]/total_MOTU[i] 
MOTU_rel_2<-as.data.frame(lapply(MOTU_rel, function(x){replace(x,x <= 0.0001,0)}))
MOTU_rel_2$ID<-MOTU$ID
MOTU_rel_3 <- MOTU_rel_2 %>% dplyr::select(ID, everything())
MOTU_rel_3$MOTUReads<-rowSums(MOTU_rel_3[,2:ncol(MOTU_rel_3)])
MOTU_rel_4<-subset(MOTU_rel_3, MOTUReads != 0)
MOTU_rel_4$MOTUReads<-NULL

# Converting dataframe back to the raw counts of the OTUS that were not removed
MOTU2<-MOTU_rel_4[,2:ncol(MOTU_rel_4)]
MOTU3<-MOTU2
for (i in 1:ncol(MOTU2))  MOTU3[,i] <- MOTU2[,i]*total_MOTU[i] 
MOTU3$ID<-MOTU_rel_4$ID
postRelSum<-MOTU3
postRelSum$ID<-NULL
sum(postRelSum)
MOTU4 <- MOTU3 %>% dplyr::select(ID, everything())

#><><><><><><><><><>
#  Annotation File 
#><><><><><><><><><>
      
annotate<-fread("All_MOTU_Annotatations.csv") 
annotate1<-annotate[,c("ID","FinalKingdom","FinalPhylum","Class","Order","Family","Genus", "Species","Calcify", "TaxaFrequency")]

# merging the annotated file with the sequence file 
annotate2<-merge(annotate1,MOTU4,by = "ID") 
# Taking only Metazoans and MacroAlgae
annotate3<-subset(annotate2,FinalKingdom == "Metazoa" | FinalKingdom == "Plantae")
# Removing MOTUs Unclassified to Phylum
annotate4<-subset(annotate3, FinalPhylum != "Unclassified") # 279 remaining - 63%

# Obtaining sequences for IDs to run through MASCE which looks for pseudogenes
MASCE<-as.data.frame(annotate4[,c("ID")])
# Merging dataframe back to SequenceCheck inorder to acquire the actual sequences for each ID
MACSE2<-merge(MASCE, SequenceCheck, by = "ID")
# creating fasta file to run through MACSE
dataframe2fas(MACSE2, file = "Sponge_MACSE.fasta")

# bringing back MACSE output based on ID
outputMACSE<-read.table("MACSE_output.txt", col.names = T)
colnames(outputMACSE)[1]<-"ID"

# Removing those IDs identified as pseudogenes from dataset
annotate5<-merge(outputMACSE, annotate4, by = "ID")

# Selecting only those MOTUs identified to Porifera
SpongeAnnotation<-subset(annotate5, FinalPhylum == "PORIFERA")
SpongeAnnotations2<-SpongeAnnotation[,!names(SpongeAnnotation) %in% c("FinalKingdom", "FinalPhylum")]
colnames(SpongeAnnotations2)[1]<-"OTU"
SpongeAnnotations3<-SpongeAnnotations2[,1:6]
# The sponge classifications from SpongeAnnotations3 dataframe were exported to examine against 
# the phylogenetic tree presented in Figure 1 - see methods
# Higher level classification were modified accordingly resulting in the "Metabarcode_Classification.csv"

# Subsampling to correct for sequence size to obtain MOTU sequence Table
subsample<-as.data.frame(t(annotate5[,11:ncol(annotate5)]))
colnames(subsample)<-annotate5$ID
subsample2<-rrarefy.perm(subsample, n=100, round.out = T)
subsample3<-as.data.table(subsample2)
# Removing any zero OTUs columns if any
subsample4<-subsample3[,colSums(subsample3 != 0) > 0, with = F]
subsample4$Sample<-as.factor(rownames(subsample2))
subsample5 <- subsample4 %>% dplyr::select(Sample, everything())

# Want only Sponge OTUs
subsample6<-as.data.frame(t(subsample5[,2:ncol(subsample5)]))
colnames(subsample6)<-subsample5$Sample
subsample6$OTU<-rownames(subsample6)
SpongeMOTU<-merge(SpongeAnnotations3,subsample6, by = "OTU")
SpongeMOTU_2<-as.data.frame(t(SpongeMOTU[,7:ncol(SpongeMOTU)]))
colnames(SpongeMOTU_2)<-SpongeMOTU$OTU
SpongeMOTU_2$Sample<-rownames(SpongeMOTU_2)
SpongeMOTU_3 <- SpongeMOTU_2 %>% dplyr::select(Sample, everything())

# Now having cleaned working sponge data frame of samples and sequences
write.csv(SpongeMOTU_3, "Metabarcode_Table.csv", row.names = F)



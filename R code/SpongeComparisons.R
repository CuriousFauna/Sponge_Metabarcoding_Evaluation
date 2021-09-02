#><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
#  Sponge Richness and Composition Comparison between 28S Barcodes and COI metabarcoding 
#><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

#### NOTE:
# Figure 2a, 2b, and 2c were brought into illustrator to create the final collectively figure
# Figure 3a, and 3b were brought into illustrator to create the final collectively figure
# Figure 4a, and 4b were brought into illustrator to create the final collectively figure

library(data.table) # data manipulation
library(dplyr) # data manipulation
library(plyr) # data manipulation
library(vegan) # stats
library(ggplot2) # graphing
library(circlize) # Color Palette
library(ComplexHeatmap) # HeatMap
library(gtools) # Ordering character/number columns

#><><><><><><><><><><><><><><><>
#   Figure 2 and Table 2
#><><><><><><><><><><><><><><><>

# Sponge BM File based on Taxonomy and 28S DNA barcoding
Sponge<-fread("BM_Table.csv")
Sponge<-Sponge[order(Sponge$Sample),]
# Total Sponge Richness
Sponge$AllSponges<-specnumber(Sponge[,2:ncol(Sponge)])

# Sponge Taxonomy File
SpongeClass<-fread("BM_Classification.csv")
SpongeClass2<-subset(SpongeClass, Class != "Calcarea")

#Creating Vector of BM with Calcarea absent
NoCalc<-SpongeClass2$OTU

# Selecting only those BM columns with the OTUs from the vector
NoCalc2<-subset(Sponge, select = NoCalc)
NoCalc2$Sample<-Sponge$Sample
NoCalc2 <- NoCalc2 %>%
  select(Sample, everything())

# Sponge Richness without Calcarea
NoCalc2$NoCalcarea<-specnumber(NoCalc2[,2:ncol(NoCalc2)])
NoCalc2<-NoCalc2[order(NoCalc2$Sample),]

# Creating Table for all Richness values
SpongeOTU<-NoCalc2[,c("Sample","NoCalcarea")]
SpongeOTU$AllSponges<-Sponge$AllSponges

# Sponge MOTU Metabarcoding File
metaBar <- fread("Metabarcode_Table.csv")
metaBar<-metaBar[order(metaBar$Sample),]
metaBar$MetaRichness<-specnumber(metaBar[,2:ncol(metaBar)])

# Pooling richness data across methods 
SpongeComparison<-SpongeOTU
SpongeComparison$MetaRichness<-metaBar$MetaRichness

#><><><><><><><><><><><><><><><>
#   Table 2 - All Sponges 
#><><><><><><><><><><><><><><><>

t.test(SpongeComparison$AllSponges,SpongeComparison$MetaRichness, paired = T)
t.test(SpongeComparison$NoCalcarea,SpongeComparison$MetaRichness, paired = T)
wilcox.test(SpongeComparison$AllSponges,SpongeComparison$MetaRichness, paired = T)
wilcox.test(SpongeComparison$NoCalcarea,SpongeComparison$MetaRichness, paired = T)


#><><><><><><><><><><><>
#   Figure 2A
#><><><><><><><><><><><>

## Melting Sponge BM
SpongeOTUMelt<-melt(SpongeOTU, id = c( "Sample"), variable = "RichnessType", value = "Richness")

# Melting Metabarcode
metaBar2<-metaBar[,c("Sample", "MetaRichness")]
metaBarMelt<-melt(metaBar2, id = c("Sample"), variable = "RichnessType", value = "Richness")

### Pooling Data
AllRichnessSpongeGraph<-rbind(SpongeOTUMelt, metaBarMelt)

### Mean Richness by Richness Calculation
RichnessAll <- ddply(AllRichnessSpongeGraph, c("RichnessType"), summarise,
                     N    = length(Richness),
                     mean = mean(Richness),
                     sd   = sd(Richness),
                     se   = sd / sqrt(N) )

# Ordering bargraphs
RichnessOrder<-c("AllSponges", "NoCalcarea","MetaRichness")
RichnessAll$RichnessType <- factor(RichnessAll$RichnessType, levels = RichnessOrder)

ggplot(RichnessAll, aes(x=RichnessType, y = mean, fill = RichnessType))+
  geom_col(position = "dodge", col = "black", binwidth = 500)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                width=.2, position=position_dodge(.9))+
  scale_fill_manual(values = c( "MetaRichness" = "skyblue","AllSponges" = "Orange", "NoCalcarea" = "peachpuff"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())



#><><><><><><><><><><><>
#   Class Comparisons
#><><><><><><><><><><><>

SpongeClassComparison<-Sponge
SpongeClassComparison$AllSponges<-NULL
SpongeClassComparison2<-melt(SpongeClassComparison, id=c("Sample"), variable = "OTU")
SpongeClassComparison3<-dcast(SpongeClassComparison2, OTU ~ Sample )

Class<-SpongeClass[,c("OTU","Class")]
SpongeClassComparison4<-merge(Class, SpongeClassComparison3, by = "OTU")
SpongeClassComparison5<-melt(SpongeClassComparison4, id = c("OTU","Class"), variable = "Sample")
SpongeClassComparison6<-dcast(SpongeClassComparison5, Sample + Class ~ OTU)
SpongeClassComparison6[is.na(SpongeClassComparison6)] <- 0
SpongeClassComparison6$ClassRichness<-specnumber(SpongeClassComparison6[,3:ncol(SpongeClassComparison6)])
SpongeClassComparison7<-SpongeClassComparison6[,c("Sample","Class","ClassRichness")]
SpongeClassComparisonNoCalcarea<-subset(SpongeClassComparison7, Class != "Calcarea")

metaBarClass<-metaBar
metaBarClass$MetaRichness<-NULL
metaBarClass2<-melt(metaBarClass, id=c("Sample"), variable = "OTU", value = "Present")
metaBarClass3<-dcast(metaBarClass2, OTU ~ Sample )

metaClass<-fread("Metabarcode_Classification.csv")
metaClass2<-metaClass[,c("OTU","Class")]
metaBarClass4<-merge(metaClass2, metaBarClass3, by = "OTU")
metaBarClass5<-melt(metaBarClass4, id = c("OTU","Class"), variable = "Sample")
metaBarClass6<-dcast(metaBarClass5, Sample + Class ~ OTU)
metaBarClass6[is.na(metaBarClass6)] <- 0
metaBarClass6$MetaClassRichness<-specnumber(metaBarClass6[,3:ncol(metaBarClass6)])
metaBarClass7<-metaBarClass6[,c("Sample","Class","MetaClassRichness")]

# Combinging Methods 
SpongeClassComparisonNoCalcarea$MetaClassRichness<-metaBarClass7$MetaClassRichness

#><><><><><><><><><><><><><><><><><>
#  Table 2- Class Comparisons
#><><><><><><><><><><><><><><><><><>

Demosponges<-subset(SpongeClassComparisonNoCalcarea, Class == "Demospongiae")
Homoscleromorphs<-subset(SpongeClassComparisonNoCalcarea, Class == "Homoscleromorpha")

t.test(Demosponges$ClassRichness,Demosponges$MetaClassRichness, paired = T)
t.test(Homoscleromorphs$ClassRichness,Homoscleromorphs$MetaClassRichness, paired = T)
wilcox.test(Demosponges$ClassRichness,Demosponges$MetaClassRichness, paired = T)
wilcox.test(Homoscleromorphs$ClassRichness,Homoscleromorphs$MetaClassRichness, paired = T)

#><><><><><><><><><><><>
#   Figure 2B
#><><><><><><><><><><><>

SpongeClassComparisonNoCalcarea2<-melt(SpongeClassComparisonNoCalcarea, id = c("Class", "Sample"), variable = "RichnessType", value = "Richness")

RichnessClass <- ddply(SpongeClassComparisonNoCalcarea2, c("Class","RichnessType"), summarise,
                       N    = length(Richness),
                       mean = mean(Richness),
                       sd   = sd(Richness),
                       se   = sd / sqrt(N))

ggplot(RichnessClass, aes(x=Class, y = mean, fill = RichnessType))+
  geom_col(position = "dodge", col = "black", binwidth = 500)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                width=.2, position=position_dodge(.9))+
  scale_fill_manual(values = c( "MetaClassRichness" = "skyblue","ClassRichness" = "peachpuff"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())


### Getting Calcarea Graph
ClassCalcarea<-subset(SpongeClassComparison7, Class == "Calcarea")

RichnessClassCalcarea <- ddply(ClassCalcarea, c("Class"), summarise,
                       N    = length(ClassRichness),
                       mean = mean(ClassRichness),
                       sd   = sd(ClassRichness),
                       se   = sd / sqrt(N))

ggplot(RichnessClassCalcarea, aes(x=Class, y = mean, fill = Class))+
  geom_col(position = "dodge", col = "black", binwidth = 500)+
   geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                 width=.2, position=position_dodge(.9))+
  scale_fill_manual(values = c( "Calcarea" = "gray"))+
  theme_bw()+
  ylim(0,8)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())


#><><><><><><><><><><><><><><><><><>
#  Order Comparisons
#><><><><><><><><><><><><><><><><><>

Order<-SpongeClass[,c("OTU","Class","Order")]
Order$Order[Order$Order==""]<-"Unknown"

SpongeOrderComparison4<-merge(Order, SpongeClassComparison3, by = "OTU")
SpongeOrderComparison5<-melt(SpongeOrderComparison4, id = c("OTU","Class","Order"), variable = "Sample")
SpongeOrderComparison6<-dcast(SpongeOrderComparison5, Sample + Order ~ OTU)
SpongeOrderComparison6[is.na(SpongeOrderComparison6)] <- 0
SpongeOrderComparison6$BarOrderRichness<-specnumber(SpongeOrderComparison6[,3:ncol(SpongeOrderComparison6)])
SpongeOrderComparison7<-SpongeOrderComparison6[,c("Sample","Order","BarOrderRichness")]

SpongeOrderComparisonNoCalcarea<-subset(SpongeOrderComparison7, Order != "Clathrinida" & Order != "Leucosolenida")
SpongeOrderComparisonNoCalcarea2<-melt(SpongeOrderComparisonNoCalcarea, id=c("Sample","Order"), variable = "RichnessType", value = "Richness")
SpongeOrderComparisonNoCalcarea2$Method<-"Barcode"

metaOrder<-metaClass[,c("OTU","Order")]
metaBarOrder4<-merge(metaOrder, metaBarClass3, by = "OTU")
metaBarOrder5<-melt(metaBarOrder4, id = c("OTU","Order"), variable = "Sample")
metaBarOrder6<-dcast(metaBarOrder5, Sample + Order ~ OTU)
metaBarOrder6[is.na(metaBarOrder6)] <- 0
metaBarOrder6$MetaOrderRichness<-specnumber(metaBarOrder6[,3:ncol(metaBarOrder6)])
metaBarOrder7<-metaBarOrder6[,c("Sample","Order","MetaOrderRichness")]
metaBarOrder8<-melt(metaBarOrder7, id=c("Sample","Order"), variable = "RichnessType", value = "Richness")
metaBarOrder8$Method<-"Metabarcode"

# Combinging Methods 
SpongeOrder<-rbind(metaBarOrder8,SpongeOrderComparisonNoCalcarea2 )
SpongeOrder2<-dcast(SpongeOrder, Sample + Order ~ RichnessType, value.var = "Richness")


#><><><><><><><><><><><><><><><><><>
#  Table 2 - Order Comparisons
#><><><><><><><><><><><><><><><><><>

Dendro<-subset(SpongeOrder2, Order == "Dendroceratida")
Dicty<-subset(SpongeOrder2, Order == "Dictyoceratida")
Haplo<-subset(SpongeOrder2, Order == "Haplosclerida")
Homo<-subset(SpongeOrder2, Order == "Homosclerophorida")
Poec<-subset(SpongeOrder2, Order == "Poecilosclerida")
Suber<-subset(SpongeOrder2, Order == "Suberitida")
Tethy<-subset(SpongeOrder2, Order == "Tethyida")
Tetra<-subset(SpongeOrder2, Order == "Tetractinellida")
Unkn<-subset(SpongeOrder2, Order == "Unknown")

t.test(Dendro$BarOrderRichness,Dendro$MetaOrderRichness, paired = T)
wilcox.test(Dendro$BarOrderRichness,Dendro$MetaOrderRichness, paired = T)

t.test(Dicty$BarOrderRichness,Dicty$MetaOrderRichness, paired = T)
wilcox.test(Dicty$BarOrderRichness,Dicty$MetaOrderRichness, paired = T)

t.test(Haplo$BarOrderRichness,Haplo$MetaOrderRichness, paired = T)
wilcox.test(Haplo$BarOrderRichness,Haplo$MetaOrderRichness, paired = T)

t.test(Homo$BarOrderRichness,Homo$MetaOrderRichness, paired = T)
wilcox.test(Homo$BarOrderRichness,Homo$MetaOrderRichness, paired = T)

t.test(Poec$BarOrderRichness,Poec$MetaOrderRichness, paired = T)
wilcox.test(Poec$BarOrderRichness,Poec$MetaOrderRichness, paired = T)

t.test(Suber$BarOrderRichness,Suber$MetaOrderRichness, paired = T)
wilcox.test(Suber$BarOrderRichness,Suber$MetaOrderRichness, paired = T)

t.test(Tethy$BarOrderRichness,Tethy$MetaOrderRichness, paired = T)
wilcox.test(Tethy$BarOrderRichness,Tethy$MetaOrderRichness, paired = T)

t.test(Tetra$BarOrderRichness,Tetra$MetaOrderRichness, paired = T)
wilcox.test(Tetra$BarOrderRichness,Tetra$MetaOrderRichness, paired = T)

#><><><><><><><><><><><><><><><><><>
#  Figure 2C 
#><><><><><><><><><><><><><><><><><>

OrderRichness <- ddply(SpongeOrder, c("Order","Method"), summarise,
                       N    = length(Richness),
                       mean = mean(Richness),
                       sd   = sd(Richness),
                       se   = sd / sqrt(N))

# Ordering bargraphs
RichnessOrd<-c("Biemnida","Chondrillida", "Dendroceratida","Dictyoceratida","Haplosclerida",
               "Poecilosclerida","Suberitida", "Tethyida",
               "Tetractinellida","Unknown","Homosclerophorida")
OrderRichness$Order <- factor(OrderRichness$Order, levels = RichnessOrd)


ggplot(OrderRichness, aes(x=Order, y = mean, fill = Method))+
  geom_col(position = "dodge", col = "black", binwidth = 500)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                width=.2, position=position_dodge(.9))+
  scale_fill_manual(values = c( "Metabarcode" = "skyblue","Barcode" = "peachpuff"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())

### Getting Calcarea Graph
OrderCalcarea<-subset(SpongeOrderComparison7, Order == "Clathrinida" | Order == "Leucosolenida")

RichnessOrderCalcarea <- ddply(OrderCalcarea, c("Order"), summarise,
                               N    = length(BarOrderRichness),
                               mean = mean(BarOrderRichness),
                               sd   = sd(BarOrderRichness),
                               se   = sd / sqrt(N))

ggplot(RichnessOrderCalcarea, aes(x=Order, y = mean, fill = Order))+
  geom_col(position = "dodge", col = "black", binwidth = 500)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                width=.2, position=position_dodge(.9))+
  scale_fill_manual(values = c( "Clathrinida" = "gray", "Leucosolenida" = "gray"))+
  theme_bw()+
  ylim(0,4)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())


#><><><><><><><><><><><><><><><><><>
#  Figure 3A
#><><><><><><><><><><><><><><><><><>

OrderGraph<-SpongeClass[,c("OTU","Class","Order")]
OrderGraph<-subset(OrderGraph, Class != "Calcarea")
OrderGraph$Order[OrderGraph$Order==""]<-"Unknown"
OrderGraph$Rep<-rep(1,nrow(OrderGraph))
OrderGraph2<-ddply(OrderGraph, .(Order), summarize, OTU=sum(Rep))
OrderGraph2$RelAbun<-OrderGraph2$OTU/sum(OrderGraph2$OTU)
OrderGraph2$Method<-paste("Barcode")

OrderGraphMeta<-metaClass[,c("OTU", "Order")]
OrderGraphMeta$Rep<-rep(1,nrow(OrderGraphMeta))
OrderGraphMeta2<-ddply(OrderGraphMeta, .(Order), summarize, OTU=sum(Rep))
OrderGraphMeta2$RelAbun<-OrderGraphMeta2$OTU/sum(OrderGraphMeta2$OTU)
OrderGraphMeta2$Method<-paste("Metabarcode")

OrderGraph3<-rbind(OrderGraph2, OrderGraphMeta2)


OrdOrd<-c("Chondrillida","Dendroceratida","Dictyoceratida","Haplosclerida",
          "Homosclerophorida","Poecilosclerida","Suberitida","Tetractinellida","Tethyida", "Unknown", "Biemnida")

OrderGraph3$Order<-factor(OrderGraph3$Order, levels = OrdOrd) 

ggplot(OrderGraph3,aes( x= Method, y = RelAbun, fill=Order))+
 geom_bar( stat='identity',colour="black")+
  scale_fill_manual(values = c( "Chondrillida" = "orange", "Dendroceratida" = "magenta", "Haplosclerida" ="pink" ,
                                "Homosclerophorida" ="green","Poecilosclerida" = "yellow","Tethyida" ="#A50026" ,
                                "Suberitida" = "blue", "Dictyoceratida" = "gray60" ,"Tetractinellida" = "cyan", 
                                "Unknown" = "black", "Biemnida" = "white"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())

#><><><><><><><><><><><><><><><><><>
#  Figure 3B
#><><><><><><><><><><><><><><><><><>

Family<-SpongeClass[,c("OTU","Class","Order","Family")]
Family2<-subset(Family, Class != "Calcarea")
Family2<-Family2[,c("OTU","Family")]
Family2$Family[Family2$Family==""]<-"Unknown"
Family2$Rep<-rep(1,nrow(Family2))
Family3<-ddply(Family2, .(Family), summarize, OTU=sum(Rep))
Family3$RelAbun<-Family3$OTU/sum(Family3$OTU)
Family3$Method<-paste("Barcode")

FamilyGraphMeta<-metaClass[,c("OTU", "Family")]
FamilyGraphMeta$Family[FamilyGraphMeta$Family==""]<-"Unknown"
FamilyGraphMeta$Rep<-rep(1,nrow(FamilyGraphMeta))
FamilyGraphMeta2<-ddply(FamilyGraphMeta, .(Family), summarize, OTU=sum(Rep))
FamilyGraphMeta2$RelAbun<-FamilyGraphMeta2$OTU/sum(FamilyGraphMeta2$OTU)
FamilyGraphMeta2$Method<-paste("Metabarcode")

Family4<-rbind(Family3, FamilyGraphMeta2)


FamMet<-c("Chalinidae","Coelosphaeridae","Darwinellidae", "Dysideidae", "Halichondriidae","Oscarellidae","Plakinidae","Suberitidae","Tedanidae",
          "Tethyidae","Unknown","Ancorinidae", "Irciniidae","Mycalidae","Chondropsidae","Biemnidae")

Family4$Family<-factor(Family4$Family, levels = FamMet) 

ggplot(Family4,aes(x=Method, y = RelAbun, fill=Family))+
  geom_bar(stat='identity',colour="black")+
  scale_fill_manual(values = c( "Chalinidae" = "pink","Darwinellidae" = "magenta", "Dysideidae" = "cyan4", "Halichondriidae" ="gold4" ,
                                "Oscarellidae" = "grey60","Plakinidae" ="green","Suberitidae" = "blue","Tedanidae" = "yellow","Tethyidae" ="#A50026" ,
                                "Unknown" = "black","Ancorinidae" = "pink4", "Irciniidae" = "red" ,"Mycalidae" = "orange", "Coelosphaeridae" = "purple",
                                "Chondropsidae" = "lightslateblue", "Biemnidae" = "white"))+
  scale_x_discrete(drop = T)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())



#><><><><><><><><><><><><><><><><><>
#  Figure 4A
#><><><><><><><><><><><><><><><><><>

# Heat maps for the Barcoded Morphologies (BM) and COI metabarcoding data were made seperately and combined together in illustrator to create the final figure.
# When a square was present in the BM heatmap but not the COI, the corresponding square in the COI heatmap was color-coded as gray to indicate a false negative
# When a square was present in the COI heatmap but not the BM heatmap, the  square in the COI heatmap was color-coded as yellow to indicate a false positive.

# Each heatmap represents presence/absence occurrence comparisons of the 23 species that were successfully barcoded with both COI and 
# 28S rRNA genes and had >98% sequence similarity match between COI barcodes and molecular operational taxonomic units (MOTUs) from COI metabarcoding (Figure 1). 

#><><><><><><><>
#  COI Map
#><><><><><><><>

# Table of the 23 COI MOTUs that matched COI barcode @ >98% ID 
HeatCOI<-read.csv("COI_heat.csv")
# Getting MOTU table
HeatCOI_OTU<-metaBar
HeatCOI_OTU$MetaRichness<-NULL

# Melting MOTU table to then merge against matching MOTUS
HeatCOI_OTU2<-melt(HeatCOI_OTU, id = c("Sample"),variable = "MOTU")
HeatCOI_OTU2$PA<-ifelse(HeatCOI_OTU2$value > 0,1,0)
HeatCOI_OTU2$Sample<-as.factor(HeatCOI_OTU2$Sample)
HeatCOI_OTU3<-dcast(HeatCOI_OTU2, MOTU ~ Sample, value.var = "PA")
HeatCOI_OTU4<-merge(HeatCOI, HeatCOI_OTU3, by = "MOTU")

# arranging MOTUs to match with 28S OTU file that comes next
HeatCOI_OTU5<-arrange(HeatCOI_OTU4, OrderMatch)
rownames(HeatCOI_OTU5)<-HeatCOI_OTU5$OrderMatch

ht_list<-colnames(HeatCOI_OTU5[,9:ncol(HeatCOI_OTU5)])
                      
HeatCOI_OTU6<-as.matrix(HeatCOI_OTU5[,9:ncol(HeatCOI_OTU5)])

col_fun=colorRamp2(c(0,1), c("white","blue"))

### Image Color codes in blue the presence of the species
Heatmap(HeatCOI_OTU6, name = "COI Metabarcoding MOTUs",
        row_order = rownames(HeatCOI_OTU6),
        column_order = mixedorder(ht_list), col = col_fun)

#><><><><><><><>
#  28S Map
#><><><><><><><>

# Table of the 23 28S rRNA BMs that corresponded the same COI barcodes that matched the COI metabarcodes @ >98% ID 
Heat28S<-read.csv("28S_heat.csv")
# Getting OTU table
Heat28S_OTU<-Sponge
Heat28S_OTU$AllSponges<-NULL

# Melting BM table to then merge against matching MOTUS
Heat28S_OTU2<-melt(Heat28S_OTU, id = c("Sample"),variable = "OTU")
Heat28S_OTU2$PA<-ifelse(Heat28S_OTU2$value > 0,1,0)
Heat28S_OTU2$Sample<-as.factor(Heat28S_OTU2$Sample)
Heat28S_OTU3<-dcast(Heat28S_OTU2, OTU ~ Sample, value.var = "PA")
Heat28S_OTU4<-merge(Heat28S, Heat28S_OTU3, by = "OTU")

# arranging MOTUs to match with 28S BM file that comes next
Heat28S_OTU5<-arrange(Heat28S_OTU4, OrderMatch)
rownames(Heat28S_OTU5)<-Heat28S_OTU5$OrderMatch

ht_list2<-colnames(Heat28S_OTU5[,9:ncol(Heat28S_OTU5)])

Heat28S_OTU6<-as.matrix(Heat28S_OTU5[,9:ncol(Heat28S_OTU5)])

Heatmap(Heat28S_OTU6, name = "Barcoded Morphologies (BMs)",
        row_order = rownames(Heat28S_OTU6),
        column_order = mixedorder(ht_list2), col = col_fun)

### Figures were then brought into Illustrator and 
# and the false negative and false positive MOTUs were
# color coded with gray and yellow respectively

#><><><><><><><><><><><><>
#  Figure 4B
#><><><><><><><><><><><><>
### Pie charts based on the following:

Match<-77
FalsePos<-35
FalseNeg<-31
Total_Metabar<-143
Total_BM<-108

# Left Pie
slice <- c(31,77)
lbl <- c("FalseNeg","Match")
pie(slice, labels = lbl, col = c("gray","blue"), main = "Total BM Occurences")

# Right Pie
slices <- c(31,77,35)
lbls <- c("FalseNeg","Match","FalsePos")
pie(slices, labels = lbls, col = c("gray","blue","yellow"), main = "All Metabarcoding Conditions")


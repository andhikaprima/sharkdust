## Modified from 
## Title: Data Analysis code for eDNA Wetmarket
## Author: Shelby E. McIlroy
## Date: 17 Feb 2022

#Load Libraries
library(phyloseq)
library(ggplot2)
library(dplyr)
library(vegan)
library(BiodiversityR)

# working directory
setwd("sharkdust")

##Load data

OTUmat<-read.csv("Data/4_species_curve/OTUmatrix.csv",header=T) #use file "OTUmatrix.csv"
METAmat<-read.csv("Data/4_species_curve/METAmatrix.csv",header=T) #use file "METAmatrix.csv"
TAXAmat<-read.csv("Data/4_species_curve/TAXAmatrix.csv",header=T) #use filt "TAXAmatrix.csv"

##Format data for phyloseq
OTU_rows<-OTUmat$Sample
OTUmat1<-OTUmat[,-1]
rownames(OTUmat1)<-OTU_rows
OTUmat2<-as.matrix(OTUmat1)
OTU<-otu_table(OTUmat2,taxa_are_rows=T)

TAXA_rows<-TAXAmat$X
TAXAmat1<-TAXAmat[,-1]
rownames(TAXAmat1)<-TAXA_rows
TAX<-tax_table(as.matrix(TAXAmat1))

META_rows<-METAmat$X
METAmat1<-METAmat[,-1]
rownames(METAmat1)<-META_rows
METAmat1$ID2<-as.character(METAmat1$ID2)
META<-sample_data(METAmat1)

##Make a phyloseq object
ps<-phyloseq(OTU,META,TAX)
ps<-prune_taxa(taxa_sums(ps)>0,ps)

##Species Accumulation Curves
ps_spec<-t(as(otu_table(ps),"matrix"))
ps_spec<-as.data.frame(ps_spec)
ps_samp<-as(sample_data(ps),"data.frame")

Accum.1 <- accumcomp(ps_spec, y=ps_samp, factor='Type', 
                     method='exact', conditioned=FALSE, plotit=FALSE)
Accum.1

accum.long1 <- accumcomp.long(Accum.1, ci=NA, label.freq=5)
head(accum.long1)

##Plot Species Accumulation Curve
# Colorblind friendly
SAC03 <- ggplot(data=accum.long1, aes(x = Sites, y = Richness, ymax = UPR, ymin = LWR)) + 
  geom_line(aes(colour=Grouping), size=1.25) +
  geom_point(data=subset(accum.long1, labelit==TRUE), 
             aes(colour=Grouping), size=0.5) +
  geom_ribbon(aes(fill=Grouping), alpha=0.2, show.legend=FALSE) + 
  theme_bw()+
  theme(plot.background = element_rect("white"),
        # text = element_text(size = 12),
        panel.background = element_rect("white"),
        panel.grid = element_line("white"),
        axis.line = element_line("gray25"),
        axis.text = element_text(size = 12, color = "gray25",  family="sans"),
        axis.title = element_text(color = "gray25",  family="sans"),
        legend.text = element_text(size = 12,  family="sans"),
        legend.key = element_rect("white")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = c(0.25,0.91),legend.background = element_rect(fill='transparent'))+
  scale_colour_manual(name = element_blank(),
                      labels = c("Tissue","Dust"),
                      values = c("#7C97D2","lightsalmon1"))  +
  scale_fill_manual(values = c("#7C97D2","lightsalmon1"))+
  labs(x = "Samples", y = "Taxa")
SAC03

# PERMANOVA, the differences between dust adn tissue
Sample <- OTUmat$Sample

OTUmat_adonis <- as.data.frame(t(OTUmat[,-1]))
colnames(OTUmat_adonis) <- Sample

METAmat1_adonis <- METAmat1

adonis(OTUmat_adonis ~ Type, data=METAmat1_adonis, permutations=999)


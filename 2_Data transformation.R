########################
# Modified from McInnes JC, Jarman SN, Lea M-A, et al. (2017c) 
# DNA metabarcoding as a marine conservation and management tool: a circumpolar 
# examination of fishery discards in the diet of threatened albatrosses. 
# Frontiers in Marine Science 4, 277.

# working directory
setwd("sharkdust")

## Pre data treatment
# 28 samples pooled to 7 based on location and using excel were squared rooted

BBA_Da=read.csv("Data/2_transformation/0_RawDustPooled_SqRt.csv",header=T,stringsAsFactors=F)
row.names(BBA_Da)=BBA_Da[,1]
BBA_D=as.matrix(BBA_Da[,-1])# get rid of info column
head(BBA_D)
dim(BBA_D)  #[1] 61 7
#Data from a single population - 7 samples and 62 shark taxa 
# amplified with fish-specific PCR

# How many seq per scat?
BBA_RC=apply(BBA_D,2,sum);BBA_RC;mean(BBA_RC)

# number of taxa against read depth
BBA_TaxC=apply(BBA_D,2,function(x) length(unique(x)))
plot(BBA_RC,BBA_TaxC)

plot(BBA_TaxC, BBA_RC)

#Turn into proportions
BBA_D
BBA_Dprop=prop.table(BBA_D,2)
dim(BBA_Dprop)
apply(BBA_Dprop,2,sum)

# Make PA dataset
PA_BBA_D= BBA_Dprop 
PA_BBA_D[PA_BBA_D>0]=1
PA_BBA_D[PA_BBA_D=0]=0

# PA_BBA_D[PA_BBA_D>.0099999999]=1
# PA_BBA_D[PA_BBA_D<.01]=0

#####  FIGURE 1a
BBA_col=rainbow(61)[seq(1,61,by=3)] # colour
##### Look at food in 10 samples
#samps=sample(ncol(BBA_D), 10, replace=F)

#samps= all (7)
samps=c(1:7)

PA_BBA_7=as.matrix(PA_BBA_D[,samps])
# barplot(PA_BBA_7,col=BBA_col,xaxt='n')

BBA_D_7=as.matrix(BBA_D[,samps])

# Plot
barplot(PA_BBA_7,col=BBA_col,xaxt='n', main = "Occurrence", ylab = "Count")
barplot(prop.table(PA_BBA_7,2)*100,col=BBA_col,xaxt='n', main = "Occurrence (weighted)", ylab = "%")
barplot(prop.table(BBA_D_7,2)*100,col=BBA_col,xaxt='n', main = "Relative read abundance", ylab = "%", xlab = "Sample")


# Stacked RRA
# library
library(ggplot2)
require(dplyr)
library(reshape2)
library(RColorBrewer)
library(viridis)

data_rra <- data.frame(prop.table(BBA_D_7,2)*100)
data_rra <- tibble::rownames_to_column(data_rra, "Species")

# rename column/ sample
names(data_rra)[names(data_rra) == "DUST_JKT_PP1_0001"] <- "dJKT1"
names(data_rra)[names(data_rra) == "DUST_IDM_PP2_000P"] <- "dIDM2"
names(data_rra)[names(data_rra) == "DUST_IDM_PP3_000P"] <- "dIDM3"
names(data_rra)[names(data_rra) == "DUST_CLP_PP4_000P"] <- "dCLP4"
names(data_rra)[names(data_rra) == "DUST_SBY_BP1_0027"] <- "dSBY5"
names(data_rra)[names(data_rra) == "DUST_SBY_PP5_000P"] <- "dSBY6"
names(data_rra)[names(data_rra) == "DUST_BYW_PP6_000P"] <- "dBYW7"


# reformat to long dataset
data_rra_long <- melt(data_rra, id.vars=c("Species"))

# rearrange species
data_rra_long$Species <- factor(data_rra_long$Species, levels = data_rra$Species)


# color
RRA_col2=rainbow(100)[seq(1,100,by=1)]
RRA_col=rainbow(100)[runif(100, min=1, max=100)]
# RRA_col=rainbow(100)[rnorm(100, mean=10, sd=2)]
df_RRA_col <- data.frame(RRA_col)
df_RRA_col1 <- as.character(df_RRA_col)
#write.csv(df_RRA_col,"df_RRA_col.csv")

#reload  color
#df_RRA_col <- read.csv("df_RRA_col.csv")
df_RRA_col <- df_RRA_col$RRA_col
df_RRA_col1 <- as.character(df_RRA_col)


# plot
ggplot(data_rra_long, aes(fill=Species, y=value, x=variable)) + 
  geom_bar(position="fill", stat="identity") +
  # scale_fill_manual(values = RRA_col) +
  labs(title="", x="Sample", y = "% of Relative Read Abundance (RRA)")+
  theme(legend.position="bottom")
  # scale_fill_viridis(discrete = TRUE, option = "viridis")
  )

ggplot(data_rra_long, aes(fill=Species, y=value, x=variable)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values = df_RRA_col1) +
  labs(title="", x="Sample", y = "% of Relative Read Abundance (RRA)")+
  theme(legend.position="bottom")
# scale_fill_viridis(discrete = TRUE, option = "viridis")
)

ggplot(data_rra_long, aes(fill=Species, y=value, x=variable)) + 
  geom_bar(position="fill", stat="identity") +
  # scale_fill_manual(values = RRA_col) +
  labs(title="", x="Sample", y = "% of Relative Read Abundance (RRA)")+
  theme(legend.position="bottom") +
  scale_fill_viridis(discrete = TRUE, option = "viridis")
)

# Create csv of the output
write.csv(PA_BBA_7,"Data/2_transformation/1_DustOccurence_sqrt.csv")
write.csv(prop.table(PA_BBA_7,2)*100,"Data/2_transformation/2_DustOccurenceWeighted_sqrt.csv")
write.csv(prop.table(BBA_D_7,2)*100,"Data/2_transformation/3_DustRRA_sqrt.csv")

RRA_test <- prop.table(BBA_D_7,2)*100

#####  FIGURE 1b
# Plot FOO barplot
FOO_All_BBA_=apply(PA_BBA_D,1,sum)/length(PA_BBA_D[1,])
Ord_BBA_=order(-FOO_All_BBA_) #
#dev.new(width=4, height=5)
barplot(FOO_All_BBA_[Ord_BBA_ ][1:15],ylim=c(0,1),col="blue",xaxt='n', xlab="", ylab = "Frequency of occurrence")

P_FOO_BBA_=  apply(PA_BBA_D,1,sum)/sum(apply(PA_BBA_D,1,sum))
P_FOO_BBA_
sum(P_FOO_BBA_)

P_wFOO_BBA_=apply(prop.table(PA_BBA_D,2),1,mean,na.rm=T)
P_wFOO_BBA_
sum(P_wFOO_BBA_)

P_RRA_BBA_=apply(prop.table(as.matrix(BBA_D),2),1,mean)
P_RRA_BBA_
sum(P_RRA_BBA_)

#dev.new(width=7, height=4)
plot(P_FOO_BBA_[Ord_BBA_][1:15],pch=15,col=2,type="b",ylim= c(0,0.3),ylab="% in sampla", xlab="Taxa",xlim=c(0,16))
points(16, sum(P_FOO_BBA_[Ord_BBA_][16:length(P_FOO_BBA_)]),pch=15,col=2,cex=1.2 )

points(P_wFOO_BBA_[Ord_BBA_][1:15],pch=16,col=3,type="b")
points(16, sum(P_wFOO_BBA_[Ord_BBA_][16:length(P_wFOO_BBA_)]),pch=16,col=3,cex=1.2 )
sum(P_wFOO_BBA_[Ord_BBA_][1:15])


points(P_RRA_BBA_[Ord_BBA_][1:15],pch=17,type="b")
points(16, sum(P_RRA_BBA_[Ord_BBA_][16:length(P_RRA_BBA_)]),pch=17,cex=1.2 )

legend('top',c('POO','wPOO','RRA'),pch=c(15,16,17),col=c(2,3,1), box.lty=0)


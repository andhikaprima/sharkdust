library(dplyr)
library(tidyr)
library(data.table)
library(vegan)
library(ggplot2)
library(tidyverse)
library(extrafont)

# working directory
setwd("sharkdust")

dust <- read.csv("Data/2_transformation/3_DustRRA_sqrt.csv",header=T,stringsAsFactors=F)
htb <- read.csv("Data/3_nmds/2_HTB_pooledPA.csv",header=T,stringsAsFactors=F)

dat_all <- merge(dust, htb, by = "Sample", all= TRUE)
dat_all[is.na(dat_all)] <- 0
dat_all


# transpose
# first remember the names
n <- dat_all$Sample

# transpose all but the first column (name)
T_dat_all <- as.data.frame(t(dat_all[,-1]))
colnames(T_dat_all) <- n

str(T_dat_all) # Check the column types

write.csv(T_dat_all,"Data/3_nmds/3_Combined data dust and tissue.csv")


# Metadata
info <- read.csv("1_metadata.csv")

# set up a "metadata" frame - will be useful for plotting later!
site_type <- info %>% 
  select(Sample,	Type,	Site,	Location, ID_ori,	ID2, Richness)

# recount the richness
site_type$Richness <- rowSums(T_dat_all)



nmds_dust_htb <- metaMDS(T_dat_all, distance="jaccard", autotransform = T)
# Square root transformation
# Wisconsin double standardization
# Run 0 stress 0.1739537 
# Run 1 stress 0.1844417 
# Run 2 stress 0.1589202 
# ... New best solution
# ... Procrustes: rmse 0.1533272  max resid 0.2503789 
# Run 3 stress 0.1592174 
# ... Procrustes: rmse 0.03196921  max resid 0.0999626 
# Run 4 stress 0.1751368 
# Run 5 stress 0.1739538 
# Run 6 stress 0.1918296 
# Run 7 stress 0.1751392 
# Run 8 stress 0.1902937 
# Run 9 stress 0.1592174 
# ... Procrustes: rmse 0.03196097  max resid 0.09992847 
# Run 10 stress 0.1872984 
# Run 11 stress 0.1905629 
# Run 12 stress 0.1877347 
# Run 13 stress 0.1844417 
# Run 14 stress 0.1844417 
# Run 15 stress 0.2019675 
# Run 16 stress 0.1751392 
# Run 17 stress 0.1905628 
# Run 18 stress 0.1921971 
# Run 19 stress 0.1905629 
# Run 20 stress 0.1844417 
# *** No convergence -- monoMDS stopping criteria:
#   18: stress ratio > sratmax
# 2: scale factor of the gradient < sfgrmin

nmds_dust_htb
# Call:
#   metaMDS(comm = T_dat_all, distance = "jaccard", autotransform = T) 
# 
# global Multidimensional Scaling using monoMDS
# 
# Data:     wisconsin(sqrt(T_dat_all)) 
# Distance: jaccard 
# 
# Dimensions: 2 
# Stress:     0.1589202 
# Stress type 1, weak ties
# No convergent solutions - best solution after 20 tries
# Scaling: centring, PC rotation, halfchange scaling 
# Species: expanded scores based on ‘wisconsin(sqrt(T_dat_all))’ 

stressplot(nmds_dust_htb)

plot(nmds_dust_htb)

# Original
#pal <- c("lightsalmon1", "palegreen4")

# Colorblind friendly
pal <- c("lightsalmon1", "#7C97D2")

clean_background <- theme(plot.background = element_rect("white"),
                          panel.background = element_rect("white"),
                          panel.grid = element_line("white"),
                          axis.line = element_line("gray25"),
                          axis.text = element_text(size = 14, color = "gray25",  family="sans"),
                          axis.title = element_text(color = "gray25",  family="sans"),
                          legend.text = element_text(size = 14,  family="sans"),
                          legend.key = element_rect("white"))

plot_df <- scores(nmds_dust_htb, display = "sites") %>% 
  as.data.frame() %>% 
  rownames_to_column("Sample") %>% 
  full_join(site_type, by = "Sample")

## Plot
plot_nmds8a <- ggplot(plot_df, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(size = Richness, color = Type)) +
  scale_color_manual(values = pal) +
  annotate("text", x = -2.5, y = 1.5, label = paste0("stress: ", format(nmds_dust_htb$stress, digits = 4)), hjust = 0) +
  geom_text(data = plot_df, size = 4, nudge_x = 0.25, nudge_y = 0, check_overlap = FALSE, aes(x = NMDS1, y = NMDS2, label = ID2)) +
  theme_bw() +
  theme(plot.background = element_rect("white"),
        # text = element_text(size = 12),
        panel.background = element_rect("white"),
        panel.grid = element_line("white"),
        axis.line = element_line("gray25"),
        axis.text = element_text(size = 12, color = "gray25",  family="sans"),
        axis.title = element_text(color = "gray25",  family="sans"),
        legend.text = element_text(size = 12,  family="sans"),
        legend.key = element_rect("white")) +
  labs(title = "NMDS")
plot_nmds8a


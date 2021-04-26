library(tidyverse)
library(vegan)
source("nes8010.R")
data(varespec)
varespec.ca <- ordi_ca(varespec)
ordi_plot(varespec.ca, display="sites")
ordi_plot(varespec.ca, display="species")
summary(varespec.ca)
varespec28 <- varespec[-9,]  # delete the ninth row (sample code 28)
varespec28.ca <- ordi_ca(varespec28)

# Now compare the plots with and without sample 28:
ordi_plot(varespec.ca, display="sites", geom="text") +
  ggtitle("Original")
ordi_plot(varespec28.ca, display="sites", geom="text") +
  ggtitle("Without sample code 28")
data(varechem)   # Load up the environmental data
summary(varechem)   # Summary of varechem
varespec.cca <- ordi_cca(varespec ~ K + P + Al + pH + Baresoil, data=varechem)
ordi_plot(varespec.cca)
ordi_plot(varespec.cca, display=c("sites", "bp"), geom="text")
ordi_plot(varespec.cca, display=c("species", "bp"), geom="text")
varespec_sites_plt <- ordi_plot(varespec.cca, display=c("sites", "bp"), geom="text")
varespec_spp_plt <- ordi_plot(varespec.cca, display=c("species", "bp"), geom="text")
multi_plot(varespec_sites_plt, varespec_spp_plt, cols=2)
varespec_spp_plt <- ordi_plot(varespec.cca, display=c("species", "bp"), geom="point")
varespec_spp_plt
ordi_identify(varespec_spp_plt)
anova(varespec.cca, by="margin")
anova(varespec.cca, by="axis")
varespec.bigcca <- ordi_cca(varespec ~ . , data=varechem)
varespec.mincca <- ordi_step(varespec.bigcca)
data(dune)
data(dune.env)
dune.rda <- ordi_rda(dune ~ A1 + Management, data=dune.env)
ordi_plot(dune.rda, display=c("sites", "cn"))
ordi_plot(dune.rda, display=c("species", "cn"))
ordi_pca(Y)
summary(dune.env)  # Notice the three categories for Use
dune.pca <- ordi_pca(dune)
data(dune)
data(dune.env)
dune.pca <- ordi_pca(dune)
dune_sco <- ordi_scores(dune.pca, display="sites")

dune_sco <- mutate(dune_sco, Use = dune.env$Use)

ggplot(dune_sco, aes(x=PC1, y=PC2, fill=Use)) +
  geom_point() +
  geom_chull(alpha=0.5) + # What do you think alpha does??
  theme_classic()  # Add this if you don't like the grey background!
ggplot(dune_sco, aes(x=PC1, y=PC2, fill=Use)) +
  geom_point() +
  geom_chull(alpha=0.5) + 
  scale_fill_discrete(name = "Field Use",
                      labels = c("Hayfield", "Hay pasture", "Pasture")) +
  theme_classic()  
library(vegan)
data(dune)
dune.dis <- vegdist(dune)
dune.dis
dune.cla <- hclust(dune.dis, method="average")
plot(dune.cla)
dune.cls <- hclust(dune.dis, method="single")
dune.clc <- hclust(dune.dis, method="complete")
cor(dune.dis, cophenetic(dune.cla))
cor(dune.dis, cophenetic(dune.cls))
cor(dune.dis, cophenetic(dune.clc))
plot(dune.clc)
rect.hclust(dune.clc, 3)
dune.grp <- cutree(dune.clc, 3)
data(dune.env)

ggplot(dune.env, aes(x=as.factor(dune.grp), y=A1)) +
  geom_boxplot() +
  xlab("Classification group") +
  ylab("Depth of A1 soil horizon") +
  theme_classic()
dune.pca <- ordi_pca(dune)

dune_sco <- ordi_scores(dune.pca, display="sites")
dune_sco <- mutate(dune_sco, group_code = as.factor(dune.grp))

ggplot(dune_sco, aes(x=PC1, y=PC2, fill=group_code)) +
  geom_point() +
  geom_chull(alpha=0.5) +
  scale_fill_discrete(name = "Classification",
                      labels = c("Group 1", "Group 2", "Group 3"))
dunespp.dis <- vegdist(t(dune), method = "raup")
dunespp.clc <- hclust(dunespp.dis, method="complete")
plot(dunespp.clc)
dune.pca <- ordi_pca(dune)
vegemite(dune, dune.pca)
# Create the dendrogram for the samples using default Bray-Curtis
quadrat_tree <- hclust(vegdist(dune), "average")

# Use Raup-Crick for the species, plus t() to transpose
spp_tree <- hclust(vegdist(t(dune), "raup"), "average")

# Create heatmap
tabasco(dune, quadrat_tree, spp_tree)
library(readxl)
rocky_shore_wide <- read_xlsx("rocky_shore_wide.xlsx", sheet="data")
rocky_shore_df <- data.frame(rocky_shore_wide)
rownames(rocky_shore_df) <- rocky_shore_df[,"Sname"]
rocky_shore_df <- rocky_shore_df[, -1] # Delete the Sname column as now tagged as rownames
# First move the rownames into a column, then convert into a tibble, and finally rename the
# new column 'rowname' into 'Sname'
rocky_spp_tbl <- rocky_spp  %>%
  rownames_to_column() %>%
  as_tibble() %>%
  rename(Sname = rowname)
rocky_spp_tbl

dune.pca <- ordi_pca(dune)
summary(dune.pca)
ordi_plot(dune.pca)
ordi_plot(dune.pca, display="sites")
ordi_plot(dune.pca, display="species")

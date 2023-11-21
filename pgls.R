#install.packages("ape")
#install.packages("geiger")
#install.packages("nlme")
#install.packages("caper")

library("ape")
library("geiger")
library("nlme")
library("phytools")
library("caper")

# --- load phylogenetic tree --- #
tree<-read.newick("./Betoideae_PlastomePhylogeny.newick")
tree

# --- load all datasets --- #
DNAtransposons_Abs_dat<-read.csv("./DNAtransposons_Abs.txt",header=TRUE, sep="\t")
DNAtransposons_Prop_dat<-read.csv("./DNAtransposons_Prop.txt",header=TRUE, sep="\t")
LTRretrotransposons_Abs_dat<-read.csv("./LTRretrotransposons_Abs.txt",header=TRUE, sep="\t")
LTRretrotransposons_Prop_dat<-read.csv("./LTRretrotransposons_Prop.txt",header=TRUE, sep="\t")
overallRep_Abs_dat<-read.csv("./overallRep_Abs.txt",header=TRUE, sep="\t")
overallRep_Prop_dat<-read.csv("./overallRep_Prop.txt",header=TRUE, sep="\t")
SatDNAs_Abs_dat<-read.csv("./SatDNAs_Abs.txt",header=TRUE, sep="\t")
SatDNAs_Prop_dat<-read.csv("./SatDNAs_Prop.txt",header=TRUE, sep="\t")


# --- check species IDs in loaded data sets --- #
#dat$acc


# --- run analysis for DNAtransposons_Abs --- #
DNAtransposons_Abs_data <- comparative.data(tree, DNAtransposons_Abs_dat, acc, vcv=TRUE, vcv.dim=3)
mod <- pgls(GenomeSize ~ DNAtransposonsMbp, DNAtransposons_Abs_data, lambda='ML')
summary(mod)

# --- run analysis for DNAtransposons_Prop --- #
DNAtransposons_Prop_data <- comparative.data(tree, DNAtransposons_Prop_dat, acc, vcv=TRUE, vcv.dim=3)
mod <- pgls(GenomeSize ~ DNAtransposonProportion, DNAtransposons_Prop_data, lambda='ML')
summary(mod)

# --- run analysis for LTRretrotransposons_Abs --- #
LTRretrotransposons_Abs_data <- comparative.data(tree, LTRretrotransposons_Abs_dat, acc, vcv=TRUE, vcv.dim=3)
mod <- pgls(GenomeSize ~ LTRMbp, LTRretrotransposons_Abs_data, lambda='ML')
summary(mod)

# --- run analysis for LTRretrotransposons_Prop --- #
LTRretrotransposons_Prop_data <- comparative.data(tree, LTRretrotransposons_Prop_dat, acc, vcv=TRUE, vcv.dim=3)
mod <- pgls(GenomeSize ~ LTRProportion, LTRretrotransposons_Prop_data, lambda='ML')
summary(mod)

# --- run analysis for overallRep_Abs --- #
overallRep_Abs_data <- comparative.data(tree, overallRep_Abs_dat, acc, vcv=TRUE, vcv.dim=3)
mod <- pgls(GenomeSize ~ RepeatMbp, overallRep_Abs_data, lambda='ML')
summary(mod)

# --- run analysis for overallRep_Prop --- #
overallRep_Prop_data <- comparative.data(tree, overallRep_Prop_dat, acc, vcv=TRUE, vcv.dim=3)
mod <- pgls(GenomeSize ~ RepeatProportion, overallRep_Prop_data, lambda='ML')
summary(mod)

# --- run analysis for SatDNAs_Abs --- #
SatDNAs_Abs_data <- comparative.data(tree, SatDNAs_Abs_dat, acc, vcv=TRUE, vcv.dim=3)
mod <- pgls(GenomeSize ~ SatelliteMbp, SatDNAs_Abs_data, lambda='ML')
summary(mod)

# --- run analysis for SatDNAs_Prop --- #
SatDNAs_Prop_data <- comparative.data(tree, SatDNAs_Prop_dat, acc, vcv=TRUE, vcv.dim=3)
mod <- pgls(GenomeSize ~ SatelliteProportions, SatDNAs_Prop_data, lambda='ML')
summary(mod)


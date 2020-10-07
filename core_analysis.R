
library(phyloseq)
library(ggplot2)
library(plyr)
library(gridExtra)
library(metagMisc)
library(tidyr)
library(ggridges)

#### first make sure you load the functions in functions.R


################ core microbiome paper



merged_7species_unrarefied<-readRDS("DATA/data_7species_unrarefied.rds") #phyloseq object containing data for 7 species
shorebird_unrarefied<-readRDS("DATA/data_shorebird_unrarefied.rds") #phyloseq object containing data for Red-necked stint


### merge for plotting (without tree, because stint were sequenced using different primers)
## this means for some part of the processing stint have to be processed seperately from the other 7 species


table<-otu_table(merged_7species_unrarefied)
map<-sample_data(merged_7species_unrarefied)
taxonomy<-tax_table(merged_7species_unrarefied)

merged_8species <-merge_phyloseq(table, map, taxonomy)


table<-otu_table(shorebird_unrarefied)
map<-sample_data(shorebird_unrarefied)
taxonomy<-tax_table(shorebird_unrarefied)

stint <-merge_phyloseq(table, map, taxonomy)

merged_8species<-merge_phyloseq(merged_8species, stint) #phylo object 8 species with no tree


############################ rarefaction curves ########################
############################ rarefaction curves ########################
############################ rarefaction curves ########################
############################ rarefaction curves ########################



sample_data(merged_8species)$Species<- factor(sample_data(merged_8species)$Species, levels = c("human", "meerkat", "deer", "carollia","spinyrat","mouselemur","flamingo","Red-necked stint"))


p <- ranacapa::ggrare(merged_8species, step = 500,   se = FALSE)+ xlim(c(0,50000))

p + xlim(c(0,50000))

p1<- p + facet_wrap(~Species, scales="free_y", ncol=4)+geom_vline(xintercept=10000)+theme_bw()+theme(legend.position = "none")

rarefaction_fig<-p1


rarefaction_fig1<- rarefaction_fig+geom_line(aes(col=Species))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  # theme(strip.background = element_blank(),strip.text.x = element_blank())+
  theme(text=element_text(size=14))+ylab("Number of ASVs")+xlab("Number of reads")+
  scale_color_brewer(palette = "Dark2")

rarefaction_fig1



############################### rarefy to 10000 ###################
############################### rarefy to 10000 ###################
############################### rarefy to 10000 ###################
############################### rarefy to 10000 ###################

merged_8species_rare<- rarefy_even_depth(merged_8species, sample.size = 10000, rngsee = 100, replace = TRUE, trimOTUs=TRUE,verbose=TRUE)

merged_8species_rare


merged_7species_rare<- rarefy_even_depth(merged_7species_unrarefied, sample.size = 10000, rngsee = 100, replace = TRUE, trimOTUs=TRUE,verbose=TRUE)
merged_7species_rare<-prune_taxa(taxa_sums(merged_7species_rare)>0, merged_7species_rare)

###############


shorebird_rare<- rarefy_even_depth(shorebird_unrarefied, sample.size = 10000, rngsee = 100, replace = TRUE, trimOTUs=TRUE,verbose=TRUE)
shorebird_rare<-prune_taxa(taxa_sums(shorebird_rare)>0, shorebird_rare)


##############################################################
##############################################################
##############################################################
##############################################################




#################### loop to calculate prevalence and abundance of every ASV per host species

Prevlist<-list()

uniq <- unique(sample_data(merged_8species_rare)$Species)

for (i in 1:length(uniq)){
  
  data_1<-subset_samples(merged_8species_rare, Species == uniq[i])
  data_1<-prune_taxa(taxa_sums(data_1)>0, data_1)
  
  occupancy_abundance<-prevalence(data_1)
  occupancy_abundance$host_species <- as.character(sample_data(data_1)$Species[1])
  occupancy_abundance$RelAbundance<- (occupancy_abundance$TotalAbundance/sum(occupancy_abundance$TotalAbundance))
  occupancy_abundance$RelPrev<-(occupancy_abundance$Prevalence / length(sample_data(data_1)$feature.id))
  occupancy_abundance$RelAbundanceInd<-(occupancy_abundance$TotalAbundance/ occupancy_abundance$Prevalence/10000)
  occupancy_abundance$ASV<-taxa_names(data_1) #this is important! Don't use 'row.names' as it sneakily adds a 1 at the end if it occurs twice
  occupancy_abundance$Sample_size<-length(unique(sample_data(data_1)$Sample))
  
  
  Prevlist[[i]]<-occupancy_abundance
  
}

#combine

occupancy_abundance_df<-do.call(rbind, Prevlist)

str(occupancy_abundance_df)

occupancy_abundance_df$host_species<-factor(occupancy_abundance_df$host_species, level = c("human", "meerkat","deer", "carollia","spinyrat","mouselemur", "flamingo", "Red-necked stint"))



P1<-ggplot(occupancy_abundance_df, aes(x = RelPrev))+
  geom_histogram( col = "black", binwidth = 0.1, aes(fill = host_species), alpha=0.5)+
  theme(axis.text = element_text(size = 10))+
  theme(axis.title = element_text(size = 14))+
  xlab("Occupancy frequency")+
  ylab("Count")+
  facet_wrap(~host_species, ncol = 8, scales = "free_y")+
  theme_bw(base_size = 16)+
  theme(strip.background = element_blank(),strip.text.x = element_blank())+
  theme(axis.text.x = element_blank())+
  theme(legend.position = "none")

P2<-ggplot(occupancy_abundance_df, aes(y = RelPrev, x = host_species))+
  geom_boxplot(aes(fill = host_species), alpha=0.5)+
  scale_y_log10()+
  theme_bw(base_size = 16)+
  ylab("Prevalence (log10 scale)")+xlab("Host species")+
  theme(legend.position = "none")+
  theme(axis.text.x = element_blank())+
  annotate("rect", ymin=summary(occupancy_abundance_df$RelPrev)[2], 
    ymax=summary(occupancy_abundance_df$RelPrev)[5], xmin=0, xmax=Inf, alpha = .2)+
  geom_boxplot(aes(fill = host_species), alpha=0.5)
 



grid.arrange(P1,P2, heights = c(1,1.5))

print(mean(occupancy_abundance_df$RelPrev))*100
print(median(occupancy_abundance_df$RelPrev))*100
print(max(occupancy_abundance_df$RelPrev))


########### test is sample size is skewing things
########### test is sample size is skewing things
########### test is sample size is skewing things
########### test is sample size is skewing things
########### test is sample size is skewing things

metadata<-data.frame(sample_data(merged_8species))

head(metadata)
unique(metadata$Species)

metadata1<-subset(metadata, Species == "meerkat" | Species== "Red-necked stint")
metadata2<-subset(metadata, Species != "meerkat" & Species!= "Red-necked stint")

metadata1$Species<-factor(metadata1$Species)
metadata2$Species<-factor(metadata2$Species)
table(metadata2$Species)

reduced<-ddply(metadata2,.(Species),function(x) x[sample(nrow(x),100),])
#reduced<-ddply(metadata,.(Species),function(x) x[sample(nrow(x),44),])

to_keep<-as.character(reduced$feature.id)

to_keep1<- as.character(metadata1$feature.id)

to_keep_final<-c(to_keep, to_keep1)


merged_8species_rare_reduced<-prune_samples(to_keep_final, merged_8species_rare)
#merged_8species_rare_reduced<-prune_samples(to_keep, merged_8species_rare)

################# repeat

Prevlist<-list()

uniq <- unique(sample_data(merged_8species_rare_reduced)$Species)

for (i in 1:length(uniq)){
  
  data_1<-subset_samples(merged_8species_rare_reduced, Species == uniq[i])
  data_1<-prune_taxa(taxa_sums(data_1)>0, data_1)
  
  occupancy_abundance<-prevalence(data_1)
  occupancy_abundance$host_species <- as.character(sample_data(data_1)$Species[1])
  occupancy_abundance$RelAbundance<- (occupancy_abundance$TotalAbundance/sum(occupancy_abundance$TotalAbundance))
  occupancy_abundance$RelPrev<-(occupancy_abundance$Prevalence / length(sample_data(data_1)$feature.id))
  occupancy_abundance$RelAbundanceInd<-(occupancy_abundance$TotalAbundance/ occupancy_abundance$Prevalence/10000)
  occupancy_abundance$ASV<-taxa_names(data_1) #this is important! Don't use 'row.names' as it sneakily adds a 1 at the end if it occurs twice
  occupancy_abundance$Sample_size<-length(unique(sample_data(data_1)$Sample))
  
  
  Prevlist[[i]]<-occupancy_abundance
  
}

#combine

occupancy_abundance_df_reduced<-do.call(rbind, Prevlist)

str(occupancy_abundance_df_reduced)

occupancy_abundance_df_reduced$host_species<-factor(occupancy_abundance_df_reduced$host_species, level = c("human", "meerkat","deer", "carollia","spinyrat","mouselemur", "flamingo", "Red-necked stint"))


P1<-ggplot(occupancy_abundance_df_reduced, aes(x = RelPrev))+
  geom_histogram( col = "black", binwidth = 0.1, aes(fill = host_species), alpha=0.5)+
  theme(axis.text = element_text(size = 10))+
  theme(axis.title = element_text(size = 14))+
  xlab("Occupancy frequency")+
  ylab("Count")+
  facet_wrap(~host_species, ncol = 8, scales = "free_y")+
  theme_bw(base_size = 16)+
  theme(strip.background = element_blank(),strip.text.x = element_blank())+
  theme(axis.text.x = element_blank())+
  theme(legend.position = "none")

P2<-ggplot(occupancy_abundance_df_reduced, aes(y = RelPrev, x = host_species))+
  geom_boxplot(aes(fill = host_species), alpha=0.5)+
  scale_y_log10()+
  theme_bw(base_size = 16)+
  ylab("Prevalence (log10 scale)")+xlab("Host species")+
  theme(legend.position = "none")+
  theme(axis.text.x = element_blank())+
  annotate("rect", ymin=summary(occupancy_abundance_df_reduced$RelPrev)[2], 
           ymax=summary(occupancy_abundance_df_reduced$RelPrev)[5], xmin=0, xmax=Inf, alpha = .2)+
  geom_boxplot(aes(fill = host_species), alpha=0.5)




grid.arrange(P1,P2, heights = c(1,1.5))

########################  sup figure #######
########################  sup figure #######
########################  sup figure #######
########################  sup figure #######

str(occupancy_abundance_df)

NotFancy <- function(l) {
  l <- format(l, scientific = FALSE)
  parse(text=l)
}

#sup fig 1

ggplot(occupancy_abundance_df, aes(x = RelPrev, y =RelAbundanceInd))+
  geom_point(aes(fill = host_species), pch=21, size=2, colour = "black", alpha = 0.7)+
  facet_wrap(~host_species, ncol = 4)+
  theme_bw(base_size = 16)+
  theme(legend.position="none")+
  theme(strip.background = element_blank(),strip.text.x = element_blank())+
  ylab("Relative abundance (occupied individuals)")+
  xlab("Prevalence")+
  scale_y_log10(labels = NotFancy)+
  geom_smooth()

#sup fig 2

ggplot(occupancy_abundance_df, aes(x = RelPrev, y =RelAbundance))+
  geom_point(aes(fill = host_species), pch=21, size=2, colour = "black", alpha = 0.7)+
  facet_wrap(~host_species, ncol = 4)+
  theme_bw(base_size = 16)+
  theme(legend.position="none")+
  theme(strip.background = element_blank(),strip.text.x = element_blank())+
  ylab("Relative abundance (all individuals)")+
  xlab("Prevalence")+
  scale_y_log10(labels = NotFancy)+
  geom_smooth()

#
grid.arrange(P1,P2, ncol = 2)





############ make list of phyloseq objects #############
############ make list of phyloseq objects #############
############ make list of phyloseq objects #############
############ make list of phyloseq objects #############

sample_data(merged_7species_rare)
phylo_list<-phyloseq_sep_variable(merged_7species_rare, variable = "Species", drop_zeroes = T)


phylo_list$stint<-shorebird_rare


##################################################alpha ##################
##################################################alpha ##################
##################################################alpha ##################
##################################################alpha ##################
##################################################alpha ##################

list1<-list()

uniq<-names(phylo_list)


thresholds<-c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)

# for loop 1 to repeat for every species

for (j in 1:length(uniq)){

data<-phylo_list[[j]]  #for loop 1 (uniq)

data<-prune_taxa(taxa_sums(data) > 0, data)


#### for loop 2 to repeat for every prevalence threshold

list2<-list()


## for loop 2 (thresholds)

for (i in 1:length(thresholds)){
 
   tryCatch({ #catch errors

data_subset<-phyloseq_filter_prevalence(data, prev.trh = thresholds[i]) #edit

alpha<-estimate_richness(data_subset, measures=c("Observed", "Shannon"))
alpha$Faiths<-metagMisc::phyloseq_phylo_div(data_subset, measures = c("PD"))$PD
alpha$BWPD<-estimate_bwpd(data_subset)$PSEs

alpha$Sample<-sample_names(data_subset)
alpha$Species<-as.character(unique(sample_data(data_subset)$Species))
alpha$Prevalence<-thresholds[i] #edit

list2[[i]]<-alpha

   }, error=function(e){})
  
}

alpha_df<-do.call(rbind, list2)

list1[[j]]<-alpha_df

}

## combine 

alpha_df_all<-do.call(rbind, list1)

head(alpha_df_all)
tail(alpha_df_all)


############# check which levels didn't work due to lack of high prevalence taxa


uniq
table(subset(alpha_df_all, Species == uniq[1])$Prevalence) #needs 0.8 and 0.9
table(subset(alpha_df_all, Species == uniq[2])$Prevalence) #fine
table(subset(alpha_df_all, Species == uniq[3])$Prevalence) #fine
table(subset(alpha_df_all, Species == uniq[4])$Prevalence)#fine
table(subset(alpha_df_all, Species == uniq[5])$Prevalence)#fine
table(subset(alpha_df_all, Species == uniq[6])$Prevalence)#fine
table(subset(alpha_df_all, Species == uniq[7])$Prevalence)#fine
table(subset(alpha_df_all, Species == "Red-necked stint")$Prevalence) #needs 0.9


############################################### manually add missing data (where zeros were generated)
############################################### manually add missing data (where zeros were generated)
############################################### manually add missing data (where zeros were generated)
############################################### manually add missing data (where zeros were generated)

## annoying but need to add the zeros to the missing data for when prevalence thresholds leave empty dataset
# ie diversity == 0

### add carollia dfs 0.8

data<-phylo_list$carollia

data<-prune_taxa(taxa_sums(data) > 0, data)

    
    data_subset<-phyloseq_filter_prevalence(data, prev.trh = 0.8) #edit
    
    alpha<-estimate_richness(data_subset, measures=c("Observed", "Shannon"))
    alpha$Faiths<-metagMisc::phyloseq_phylo_div(data_subset, measures = c("PD"))$PD
    alpha$Faiths<-0
    alpha$BWPD<-0
    
    alpha$Sample<-sample_names(data_subset)
    alpha$Species<-as.character(unique(sample_data(data_subset)$Species))
    alpha$Prevalence<-0.8 #edit
    
alpha_df_all1<-rbind(alpha_df_all, alpha)

### add carollia dfs 0.9 (all zeros)

data<-phylo_list$carollia

data<-prune_taxa(taxa_sums(data) > 0, data)


#data_subset<-phyloseq_filter_prevalence(data, prev.trh = 0.9) #edit

alpha0.9<-alpha

alpha0.9$Observed<-0
alpha0.9$Shannon<-0
alpha0.9$Faiths<-metagMisc::phyloseq_phylo_div(data_subset, measures = c("PD"))$PD
alpha0.9$Faiths<-0
alpha0.9$BWPD<-0

alpha0.9$Sample<-sample_names(data_subset)
alpha0.9$Species<-as.character(unique(sample_data(data_subset)$Species))
alpha0.9$Prevalence<-0.9 #edit

alpha_df_all2<-rbind(alpha_df_all1, alpha0.9)

### add stint 0.9

data<-phylo_list$stint

data<-prune_taxa(taxa_sums(data) > 0, data)


data_subset<-phyloseq_filter_prevalence(data, prev.trh = 0.9) #edit

alpha<-estimate_richness(data_subset, measures=c("Observed", "Shannon"))
alpha$Faiths<-metagMisc::phyloseq_phylo_div(data_subset, measures = c("PD"))$PD
alpha$Faiths<-0
alpha$BWPD<-0

alpha$Sample<-sample_names(data_subset)
alpha$Species<-as.character(unique(sample_data(data_subset)$Species))
alpha$Prevalence<-0.9 #edit

alpha_df_all3<-rbind(alpha_df_all2, alpha)

table(subset(alpha_df_all3, Species == uniq[1])$Prevalence) #needs 0.8 and 0.9
table(subset(alpha_df_all3, Species == uniq[2])$Prevalence) #fine
table(subset(alpha_df_all3, Species == uniq[3])$Prevalence) #fine
table(subset(alpha_df_all3, Species == uniq[4])$Prevalence)#fine
table(subset(alpha_df_all3, Species == uniq[5])$Prevalence)#fine
table(subset(alpha_df_all3, Species == uniq[6])$Prevalence)#fine
table(subset(alpha_df_all3, Species == uniq[7])$Prevalence)#fine
table(subset(alpha_df_all3, Species == "Red-necked stint")$Prevalence) #needs 0.9

alpha_df_all<-alpha_df_all3

head(alpha_df_all)


### change NAs into zeros

table(is.na(alpha_df_all))

alpha_df_all$Faiths[is.na(alpha_df_all$Faiths)] <- 0
alpha_df_all$BWPD[is.na(alpha_df_all$BWPD)] <- 0


################# generate standarised values per species
################# generate standarised values per species
################# generate standarised values per species

alpha_df_all<-transform(alpha_df_all, Observed_scaled=ave(Observed, Species, FUN=scale))
alpha_df_all<-transform(alpha_df_all, Faiths_scaled=ave(Faiths, Species, FUN=scale))
alpha_df_all<-transform(alpha_df_all, Shannon_scaled=ave(Shannon, Species, FUN=scale))
alpha_df_all<-transform(alpha_df_all, BWPD_scaled=ave(BWPD, Species, FUN=scale))

head(alpha_df_all)

alpha_df_all<-alpha_df_all[,c(5,6,7,1,2,3,4,8,9,10,11)]





##################### BETA DIVERSITY #######################
##################### BETA DIVERSITY #######################
##################### BETA DIVERSITY #######################
##################### BETA DIVERSITY #######################
##################### BETA DIVERSITY #######################
##################### BETA DIVERSITY #######################
##################### BETA DIVERSITY #######################
##################### BETA DIVERSITY #######################


list1<-list()

uniq<-names(phylo_list)


thresholds<-c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)

# for loop 1 to repeat for every species

for (j in 1:length(thresholds)){
  
  data<-phylo_list[[j]]  #for loop 1 (uniq)
 # data<-phylo_list$carollia  #for loop 1 (uniq)
  
  data<-prune_taxa(taxa_sums(data) > 0, data)
  
  
  #### for loop 2 to repeat for every prevalence threshold
  
  list2<-list()

  
  ## for loop 2 (thresholds)
  
  for (i in 1:length(thresholds)){
    
    tryCatch({ #catch errors
      
      data_subset<-phyloseq_filter_prevalence(data, prev.trh = thresholds[i]) #edit
    #  data_subset<-phyloseq_filter_prevalence(data, prev.trh = 0.5) #edit
      
      unifrac<-as.matrix(phyloseq::distance(data_subset, method = "unifrac"))
      wunifrac<-as.matrix(phyloseq::distance(data_subset, method = "wunifrac"))
      morisita<-as.matrix(phyloseq::distance(data_subset, method = "morisita"))
      jaccard<-as.matrix(phyloseq::distance(data_subset, method = "jaccard"))
      
      beta<-data.frame(colMeans(unifrac, na.rm=T))
      names(beta)<-"unifrac"
      beta$wunifrac<-colMeans(wunifrac, na.rm=T)
      beta$morisita<-colMeans(morisita, na.rm=T)
      beta$jaccard<-colMeans(jaccard, na.rm=T)
      
      beta$Sample<-sample_names(data_subset)
      beta$Species<-as.character(unique(sample_data(data_subset)$Species))
      beta$Prevalence<-thresholds[i] #edit
      
      list2[[i]]<-beta
      
    }, error=function(e){})
    
  }
  
  beta_df<-do.call(rbind, list2)
  
  list1[[j]]<-beta_df
  
}



beta_df_all<-do.call(rbind, list1)

head(beta_df_all)



table(subset(beta_df_all, Species == uniq[1])$Prevalence) #needs and 0.8 and 0.9
table(subset(beta_df_all, Species == uniq[2])$Prevalence) 
table(subset(beta_df_all, Species == uniq[3])$Prevalence) #fine
table(subset(beta_df_all, Species == uniq[4])$Prevalence)#fine
table(subset(beta_df_all, Species == uniq[5])$Prevalence)#fine
table(subset(beta_df_all, Species == uniq[6])$Prevalence)#fine
table(subset(beta_df_all, Species == uniq[7])$Prevalence)#fine
table(subset(beta_df_all, Species == "Red-necked stint")$Prevalence) #needs 0.9


############################### missing values ###################
############################### missing values ###################
############################### missing values ###################
############################### missing values ###################
############################### missing values ###################


###add carollia 0.8
###add carollia 0.8
###add carollia 0.8


data<-phylo_list$carollia

data<-prune_taxa(taxa_sums(data) > 0, data)

data_subset<-phyloseq_filter_prevalence(data, prev.trh = 0.8)

unifrac<-as.matrix(phyloseq::distance(data_subset, method = "unifrac"))
wunifrac<-as.matrix(phyloseq::distance(data_subset, method = "wunifrac"))
morisita<-as.matrix(phyloseq::distance(data_subset, method = "morisita"))
jaccard<-as.matrix(phyloseq::distance(data_subset, method = "jaccard"))

beta<-data.frame(colMeans(jaccard, na.rm=T))
names(beta)<-"unifrac"
beta$unifrac<-0
beta$wunifrac<-0
beta$morisita<-0
beta$jaccard<-colMeans(jaccard, na.rm=T)

beta$Sample<-sample_names(data_subset)
beta$Species<-as.character(unique(sample_data(data_subset)$Species))
beta$Prevalence<-0.8

beta_df_all2<-rbind(beta_df_all, beta)


######## add carollia 0.9
######## add carollia 0.9
######## add carollia 0.9


data<-phylo_list$carollia

data<-prune_taxa(taxa_sums(data) > 0, data)

data_subset<-phyloseq_filter_prevalence(data, prev.trh = 0.9)

unifrac<-as.matrix(phyloseq::distance(data_subset, method = "unifrac"))
wunifrac<-as.matrix(phyloseq::distance(data_subset, method = "wunifrac"))
morisita<-as.matrix(phyloseq::distance(data_subset, method = "morisita"))
jaccard<-as.matrix(phyloseq::distance(data_subset, method = "jaccard"))


beta$unifrac<-0
beta$wunifrac<-0
beta$morisita<-0
beta$jaccard<-0

beta$Sample<-sample_names(data_subset)
beta$Species<-as.character(unique(sample_data(data_subset)$Species))
beta$Prevalence<-0.9

beta_df_all3<-rbind(beta_df_all2, beta)


######## add stint 0.9
######## add stint 0.9
######## add stint 0.9


data<-phylo_list$stint

data<-prune_taxa(taxa_sums(data) > 0, data)

data_subset<-phyloseq_filter_prevalence(data, prev.trh = 0.9)

unifrac<-as.matrix(phyloseq::distance(data_subset, method = "unifrac"))
wunifrac<-as.matrix(phyloseq::distance(data_subset, method = "wunifrac"))
morisita<-as.matrix(phyloseq::distance(data_subset, method = "morisita"))
jaccard<-as.matrix(phyloseq::distance(data_subset, method = "jaccard"))


beta<-data.frame(colMeans(jaccard, na.rm=T))
names(beta)<-"unifrac"
beta$unifrac<-0
beta$wunifrac<-0
beta$morisita<-0
beta$jaccard<-colMeans(jaccard, na.rm=T)

beta$Sample<-sample_names(data_subset)
beta$Species<-as.character(unique(sample_data(data_subset)$Species))
beta$Prevalence<-0.9

beta_df_all4<-rbind(beta_df_all3, beta)


############ check all prevalence levels are there

table(subset(beta_df_all4, Species == uniq[1])$Prevalence) #needs and 0.8 and 0.9
table(subset(beta_df_all4, Species == uniq[2])$Prevalence) 
table(subset(beta_df_all4, Species == uniq[3])$Prevalence) #fine
table(subset(beta_df_all4, Species == uniq[4])$Prevalence)#fine
table(subset(beta_df_all4, Species == uniq[5])$Prevalence)#fine
table(subset(beta_df_all4, Species == uniq[6])$Prevalence)#fine
table(subset(beta_df_all4, Species == uniq[7])$Prevalence)#fine
table(subset(beta_df_all4, Species == "Red-necked stint")$Prevalence) #needs 0.9


beta_df_all<-beta_df_all4




################# generate standarised values per species
################# generate standarised values per species
################# generate standarised values per species




########################## FIGURES ###############################
########################## FIGURES ###############################



##################################### beta diversity ######################
##################################### beta diversity ######################
##################################### beta diversity ######################
##################################### beta diversity ######################


head(beta_df_all)

beta_df_all$Species<-factor(beta_df_all$Species, level = c("human", "meerkat","deer", "carollia","spinyrat","mouselemur", "flamingo", "Red-necked stint"))


beta_df_all$Prevalence<-factor(beta_df_all$Prevalence)

beta_short_normal<-beta_df_all[,c(1:7)]


beta_long_normal<-gather(beta_short_normal, Index, Distance, jaccard, unifrac,  morisita,  wunifrac, factor_key = TRUE)


#beta1$Index<-factor(beta1$Index, level = c("Jaccard_dissimilarity","Unifrac_dissimilarity","BrayCurtis_dissimilarity", "Weighted_unifrac_dissimilarity"))
beta_long_normal$Index<-factor(beta_long_normal$Index, level = c("jaccard", "unifrac",  "morisita",  "wunifrac"))


##############################

ggplot(beta_long_normal, aes(x =Prevalence, y = Distance))+
  geom_jitter(aes(fill = Species), pch=21, size=1, colour = "black", alpha = 0.7, width =0.3)+
 # geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  facet_wrap(~Index+Species, ncol = 8)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme_bw(base_size = 16)+
  ylab("")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  # theme( axis.text.y = element_text(size=10), axis.text.x = element_text(size=8), axis.title=element_text(size=14))+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+theme(legend.position="none")+
  theme(legend.position="none")+
  theme(strip.background = element_blank(), strip.text.x = element_blank())





##################################################### alpha and beta correlations ########
##################################################### alpha and beta correlations ########
##################################################### alpha and beta correlations ########
##################################################### alpha and beta correlations ########

unique(beta_df_all$Prevalence)


beta_sub<-subset(beta_short_normal, Prevalence == "0")
alpha_sub<-subset(alpha_short_normal, Prevalence == "0")

head(alpha_sub)

alpha_sub<-transform(alpha_sub, Observed_scaled=ave(Observed, Species, FUN=scale))
alpha_sub<-transform(alpha_sub, Faiths_scaled=ave(Faiths, Species, FUN=scale))
alpha_sub<-transform(alpha_sub, Shannon_scaled=ave(Shannon, Species, FUN=scale))
alpha_sub<-transform(alpha_sub, BWPD_scaled=ave(BWPD, Species, FUN=scale))



merged<-beta_sub



library(expss)
library(ggcorrplot)

merged$Observed<-vlookup(merged$Sample, alpha_sub, lookup_column = "Sample", result_column = "Observed_scaled")
merged$Shannon<-vlookup(merged$Sample, alpha_sub, lookup_column = "Sample", result_column = "Shannon_scaled")
merged$Faiths<-vlookup(merged$Sample, alpha_sub, lookup_column = "Sample", result_column = "Faiths_scaled")
merged$BWPD<-vlookup(merged$Sample, alpha_sub, lookup_column = "Sample", result_column = "BWPD_scaled")


head(merged)

merged_sub<-merged[,c(4:11)]

head(merged_sub)

merged_sub<-na.omit(merged_sub)

merged_sub<-merged_sub[,c(2,3,1,4,8,6,7,5)]

names(merged_sub)[1]<-"W Unifrac"
names(merged_sub)[2]<-"Morisita"
names(merged_sub)[3]<-"UW Unifrac"
names(merged_sub)[4]<-"Jaccard"

corr <- round(cor(merged_sub), 1)
p.mat <- cor_pmat(merged_sub)

ggcorrplot(corr, lab = T, p.mat = p.mat)+ 
  scale_fill_gradient2(limit = c(-1,1), low = "navyblue", high =  "darkred", mid = "white", midpoint = 0)+
  theme(legend.position = "none")


####################

corr_samples<-beta[,c(1,2,3)]

str(corr_samples)

subset(corr_samples, Sample == "1094")

corr_samples_long <- spread(corr_samples, Prevalence, Unifrac_dissimilarity)
data_long

############################################ stats ########################################################
############################################ stats ########################################################
############################################ stats ########################################################
############################################ stats ########################################################
############################################ stats ########################################################



library(glmmTMB)
library(effects)
library(scales)
library(gamm4)

hist(alpha_short_scaled$Observed_scaled)

head(alpha_long_scaled)



alpha_short_scaled$Prevalence<-factor(alpha_short_scaled$Prevalence)
alpha_short_scaled$Prevalence<-as.numeric(alpha_short_scaled$Prevalence)


#########################

observed_model<-glmmTMB(Observed_scaled~Prevalence +Species+ (1|Sample), 
  family = gaussian, 
  data = alpha_short_scaled)


summary(observed_model)
#Prevalence  -2.80702    0.01531  -183.3   <2e-16 ***

#########################



faiths_model<-glmmTMB(Faiths_scaled~Prevalence+Species+ (1|Sample), 
  family = gaussian, 
  data = alpha_short_scaled)


summary(faiths_model)
#Prevalence  -2.995888   0.012824  -233.6   <2e-16 ***


##############################


shannon_model<-glmmTMB(Shannon_scaled~Prevalence+Species+ (1|Sample), 
  family = gaussian, 
  data = alpha_short_scaled)


summary(shannon_model)
#Prevalence  -2.29062    0.01668 -137.32   <2e-16 ***

##################################

BWPD_model<-glmmTMB(BWPD_scaled~Prevalence+Species+ (1|Sample), 
  family = gaussian, 
  data = alpha_short_scaled)


summary(BWPD_model)

#Prevalence  -0.37078    0.02147 -17.267  < 2e-16 ***

#####################################
#####################################
#####################################


######################### beta

#beta_short_normal$Prevalence<-as.numeric(beta_short_normal$Prevalence) #for slopes
beta_short_normal$Prevalence<-factor(beta_short_normal$Prevalence)



jaccard_model<-glmmTMB(jaccard~Prevalence +Species+ (1|Sample), 
  family = gaussian, 
  data = beta_short_normal)



summary(jaccard_model)

#########################



unifrac_model<-glmmTMB(unifrac~Prevalence+Species+ (1|Sample), 
  family = gaussian, 
  data = beta_short_normal)


summary(unifrac_model)


##############################

morisita_model<-glmmTMB(morisita~Prevalence+Species+ (1|Sample), 
  family = gaussian, 
  data = beta_short_normal)


summary(morisita_model)


##################################

wunifrac_model<-glmmTMB(wunifrac~Prevalence+Species+ (1|Sample), 
  family = gaussian, 
  data = beta_short_normal)


summary(wunifrac_model)

######################################plot
######################################plot
######################################plot
######################################plot


P1<-plot(effects::Effect("Prevalence",observed_model), main = "Observed", ylim=c(-2,2), xlab = "", ylab="")
P2<-plot(effects::Effect("Prevalence",faiths_model), main = "Faiths", ylim=c(-2,2), xlab = "", ylab="")
P3<-plot(effects::Effect("Prevalence",shannon_model), main = "Shannon", ylim=c(-2,2), xlab = "", ylab="")
P4<-plot(effects::Effect("Prevalence",BWPD_model), main = "BWPD", ylim=c(-2,2), xlab = "", ylab="")

P5<-plot(effects::Effect("Prevalence",jaccard_model), main = "Jaccard",  ylim=c(0,1), xlab = "", ylab="")
P6<-plot(effects::Effect("Prevalence",unifrac_model), main = "Unifrac",  ylim=c(0,1), xlab = "", ylab="")
P7<-plot(effects::Effect("Prevalence",morisita_model), main = "Morisita",  ylim=c(0,1), xlab = "", ylab="")
P8<-plot(effects::Effect("Prevalence",wunifrac_model), main = "Weighted Unifrac",  ylim=c(0,1), xlab = "", ylab="")

#library(gridExtra)
grid.arrange(P1,P2,P3,P4, P5, P6, P7, P8, ncol=4)



##################################################

head(alpha_long_scaled)

mean_alpha<-ddply(alpha_long_scaled, .(Species,Prevalence,Index), summarize, mean= mean(Measure))

head(mean_alpha)

P1<-ggplot(alpha_long_scaled, aes(x =Prevalence, y = Measure, group = Species))+
 # geom_point(aes(fill = Species), pch=21, size=3, colour = "black", alpha = 0.7)+
 # geom_line(alpha = 0.5)+
  facet_wrap(~Index, ncol = 8)+
  # theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme_bw(base_size = 16)+
  ylab("")+xlab("")+  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  # theme( axis.text.y = element_text(size=10), axis.text.x = element_text(size=8), axis.title=element_text(size=14))+
  # theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+ 
  theme(panel.grid.major = element_blank(),  panel.grid.minor = element_blank())+theme(legend.position="none")+
  theme(legend.position="none")+
  theme(strip.background = element_blank(),
    strip.text.x = element_blank())+
  geom_smooth(aes(col = Species), alpha = 0.2)+
  geom_vline(xintercept = "0.7", linetype = "dashed")+
  theme(axis.text.x = element_text(angle = 45))




###########################



P2<-ggplot(beta_long_normal, aes(x =Prevalence, y = Distance, group = Species))+
 # geom_point(aes(fill = Species), pch=21, size=3, colour = "black", alpha = 0.7)+
  #geom_line(alpha = 0.5)+
  facet_wrap(~Index, ncol = 8)+
  # theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme_bw(base_size = 16)+
  ylab("")+
  xlab("")+  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  # theme( axis.text.y = element_text(size=10), axis.text.x = element_text(size=8), axis.title=element_text(size=14))+
  # theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+ 
  theme(panel.grid.major = element_blank(),  panel.grid.minor = element_blank())+
  theme(legend.position="none")+
  theme(strip.background = element_blank(),strip.text.x = element_blank())+
  ylim(c(0,1))+
  geom_smooth(aes(col = Species), alpha = 0.2)+
  geom_vline(xintercept = "0.7", linetype = "dashed")+
  theme(axis.text.x = element_text(angle = 45))

grid.arrange(P1,P2, ncol = 1)

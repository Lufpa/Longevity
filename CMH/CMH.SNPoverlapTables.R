### Longevity experiment
### Making counts table where SNPs that overlap all replicates (2 or 3) 
### This table will be used to downsample 1 read per individual (CMH.downsampling.R)
### The downsampled tables will be used for CMH

#setwd("C:/Users/lamaya/Documents/ViabilityExperiment/Analysis/CMH/chr2R/")

#### LOAD DATA FOR T0 CAGE A AND MAKE TWO RANDOMIZED DATASETS WITH HALF THE INDIVDIUALS EACH
#### IN THIS WAY WE'LL HAVE TWO INDEPENDENT T0 CAGE A DATASETS TO RUN CMH WITHOUT VIOLATING INDEPENDENCE OF EACH REPLICATE COMPARISON

  # chr2R - CAGE A T0
  {
# Cage A T0, batch 1 and 2 (Amanda split the SNP calling in two batches of ~500 individuals)
#cage_a_T0_alt1 <- read.table("../alternate_counts_2R_cageA.T0.batch1.txt",h=T, check.names = F)
#cage_a_T0_alt2 <- read.table("../alternate_counts_2R_cageA.T0.batch2.txt",h=T, check.names = F)
#cage_a_T0_ref1 <- read.table("../reference_counts_2R_cageA.T0.batch1.txt",h=T, check.names = F)
#cage_a_T0_ref2 <- read.table("../reference_counts_2R_cageA.T0.batch2.txt",h=T, check.names = F)

# Combine the two T0 batches
#cage_a_T0_alt <- merge(cage_a_T0_alt1, cage_a_T0_alt2, by="site", all=T)
#cage_a_T0_ref <- merge(cage_a_T0_ref1, cage_a_T0_ref2, by ="site", all=T)

#  Make two batches again randomly sampling from the full dataset, to get two independent T0 cages
#cage_a_T0_alt1 <- sample(colnames(cage_a_T0_alt[,-1]), round(ncol(cage_a_T0_alt[,-1])/2,0)) #550 indiv
#cage_a_T0_alt2 <- setdiff(colnames(cage_a_T0_alt[,-1]), cage_a_T0_alt1) #551 indiv
#cage_a_T0_ref1 <- cage_a_T0_alt1
#cage_a_T0_ref2 <- cage_a_T0_alt2

#cage_a_T0_alt1 <- cage_a_T0_alt[,colnames(cage_a_T0_alt) %in% cage_a_T0_alt1]
#cage_a_T0_alt1 <- cbind.data.frame("site"=cage_a_T0_alt$site, cage_a_T0_alt1)
#cage_a_T0_alt2 <- cage_a_T0_alt[,colnames(cage_a_T0_alt) %in% cage_a_T0_alt2]
#cage_a_T0_alt2 <- cbind.data.frame("site"=cage_a_T0_alt$site, cage_a_T0_alt2)
#cage_a_T0_ref1 <- cage_a_T0_ref[,colnames(cage_a_T0_ref) %in% cage_a_T0_ref1]
#cage_a_T0_ref1 <- cbind.data.frame("site"=cage_a_T0_alt$site, cage_a_T0_ref1)
#cage_a_T0_ref2 <- cage_a_T0_ref[,colnames(cage_a_T0_ref) %in% cage_a_T0_ref2]
#cage_a_T0_ref2 <- cbind.data.frame("site"=cage_a_T0_alt$site, cage_a_T0_ref2)

#rm(cage_a_T0_alt, cage_a_T0_ref)

# write the data
#write.table(cage_a_T0_alt1, "alternate_counts_2R_cageA.T0_1.randomized.txt", sep="\t", quote=F, row.names = F, col.names = T)
#write.table(cage_a_T0_alt2, "alternate_counts_2R_cageA.T0_2.randomized.txt", sep="\t", quote=F, row.names = F, col.names = T)
#write.table(cage_a_T0_ref1, "reference_counts_2R_cageA.T0_1.randomized.txt", sep="\t", quote=F, row.names = F, col.names = T)
#write.table(cage_a_T0_ref2, "reference_counts_2R_cageA.T0_2.randomized.txt", sep="\t", quote=F, row.names = F, col.names = T)
}

 
#### LOAD DATA FOR T0 CAGE D

  # chr2R - CAGE D T0
  {
#cage_d_T0_alt1 <- read.table("../alternate_counts_2R_cageD.T0.batch1.txt", h=T, check.names = F)
#cage_d_T0_alt2 <- read.table("../alternate_counts_2R_cageD.T0.batch2.txt", h=T, check.names = F)
#cage_d_T0_ref1 <- read.table("../reference_counts_2R_cageD.T0.batch1.txt", h=T, check.names = F)
#cage_d_T0_ref2 <- read.table("../reference_counts_2R_cageD.T0.batch2.txt", h=T, check.names = F)

# Combine the two T0 batches
#cage_d_T0_alt <- merge(cage_d_T0_alt1, cage_d_T0_alt2, by="site", all=T)
#cage_d_T0_ref <- merge(cage_d_T0_ref1, cage_d_T0_ref2, by ="site", all=T)

#rm(cage_d_T0_alt1, cage_d_T0_alt2)
#rm(cage_d_T0_ref1, cage_d_T0_ref2)

# write the data
#write.table(cage_d_T0_alt, "alternate_counts_2R_cageD.T0.merged.txt", sep="\t", quote=F, row.names = F, col.names = T)
#write.table(cage_d_T0_ref, "reference_counts_2R_cageD.T0.merged.txt", sep="\t", quote=F, row.names = F, col.names = T)
}

### LOAD T0 CAGES AFTER SHUFFLING AT0, AND MERGING DT0
cage_d_T0_alt <- read.table("T0newFiles/alternate_counts_2R_cageD.T0.merged.txt", h=T, check.names=F)
cage_d_T0_ref <- read.table("T0newFiles/reference_counts_2R_cageD.T0.merged.txt", h=T, check.names=F)

cage_a_T0_alt1 <- read.table("T0newFiles/alternate_counts_2R_cageA.T0_1.randomized.txt", h=T, check.names=F)
cage_a_T0_ref1 <- read.table("T0newFiles/reference_counts_2R_cageA.T0_1.randomized.txt", h=T, check.names=F)

cage_a_T0_alt2 <- read.table("T0newFiles/alternate_counts_2R_cageA.T0_2.randomized.txt", h=T, check.names=F)
cage_a_T0_ref2 <- read.table("T0newFiles/reference_counts_2R_cageA.T0_2.randomized.txt", h=T, check.names=F)

  

#### LOAD DATA FOR CAGE A CONTROL 1 (Tend)

  # chr 4 - CA1
{
cage_a_c1_alt <- read.table("../alternate_counts_2R_cageA.C1.txt",h=T, check.names = F)
cage_a_c1_ref <- read.table("../reference_counts_2R_cageA.C1.txt",h=T, check.names = F)
}

#### LOAD DATA FOR CAGE A CONTROL 2 (Tend)
{
  cage_a_c2_alt <- read.table("../alternate_counts_2R_cageA.C2.txt",h=T, check.names = F)
  cage_a_c2_ref <- read.table("../reference_counts_2R_cageA.C2.txt",h=T, check.names = F)
}

#### LOAD DATA FOR CAGE A HIGHSUGAR 1 (Tend)
{
  cage_a_hs1_alt <- read.table("../alternate_counts_2R_cageA.HS1.txt",h=T, check.names = F)
  cage_a_hs1_ref <- read.table("../reference_counts_2R_cageA.HS1.txt",h=T, check.names = F)
}

#### LOAD DATA FOR CAGE A HIGHSUGAR 2 (Tend)
{
  cage_a_hs2_alt <- read.table("../alternate_counts_2R_cageA.HS2.txt",h=T, check.names = F)
  cage_a_hs2_ref <- read.table("../reference_counts_2R_cageA.HS2.txt",h=T, check.names = F)
}


#### LOAD DATA FOR CAGE D CONTROL 1 (Tend)
{
  cage_d_c1_alt <- read.table("../alternate_counts_2R_cageD.C1.txt",h=T, check.names = F)
  cage_d_c1_ref <- read.table("../reference_counts_2R_cageD.C1.txt",h=T, check.names = F)
}

#### LOAD DATA FOR CAGE D HIGHSUGAR 1 (Tend)
{
  cage_d_hs1_alt <- read.table("../alternate_counts_2R_cageD.HS1.txt",h=T, check.names = F)
  cage_d_hs1_ref <- read.table("../reference_counts_2R_cageD.HS1.txt",h=T, check.names = F)
}  

# removing flies from ref tables
allfiles<- c("cage_a_c1_ref","cage_a_c2_ref","cage_a_T0_ref1","cage_a_T0_ref2","cage_d_T0_ref","cage_d_c1_ref","cage_a_hs1_ref","cage_a_hs2_ref","cage_d_hs1_ref")
#for (i in allfiles){ t <- get(i); print(i); print(dim(t))} #dimensions of each table, 1st col is SITE
#### FILTER INDIVIDUALS WITH >90% MISSING DATA
#{
missperfly <- as.list(NULL)
for (i in allfiles){
  t <- get(i)
  missperfly[[i]]  <- apply(t, 2, function(x) sum(is.na(x)/nrow(t)))
  }
#par(mfrow=c(3,3))
#for (i in 1:length(missperfly)){ hist(missperfly[[i]], main=names(missperfly)[i])}
fliestokeep<- lapply(missperfly, function(x) x <- x[x<0.9])
n<- lapply(missperfly, length)
n_filtered <- lapply(fliestokeep, length)
print("total and filtered flies <90% missing data")
print (cbind(unlist(n),unlist(n_filtered)))
#plot(unlist(n) ~ unlist(n_filtered))
#abline(0,1)

#for (i in 1:length(allfiles)){
#  t<-allfiles[i]
#  tt<-get(t)
#  tt <- tt[,colnames(tt) %in% names(fliestokeep[[i]])] #only ref files were modified, dont forget!
#  assign(t,tt)
#    }
#}

# removing flies from alt tables
#allfilesALT<- c("cage_a_c1_alt","cage_a_c2_alt","cage_a_T0_alt1","cage_a_T0_alt2","cage_d_T0_alt","cage_d_c1_alt","cage_a_hs1_alt","cage_a_hs2_alt","cage_d_hs1_alt")
#for (i in 1:length(allfilesALT)){
#  t<-allfilesALT[i]
#  tt<-get(t)
#  tt <- tt[,colnames(tt) %in% names(fliestokeep[[i]])] #only alt files were modified, dont forget!
#  assign(t,tt)
#}


### OVERLAP SNPS BETWEEN CAGES TO BE RUN TOGETHER IN CMH (CTRL VS HS)

# CTRL cages

#rep1<-intersect(cage_a_c1_ref$site,cage_a_T0_ref1$site) #1837 SNPs
#rep2<- intersect(cage_a_c2_ref$site, cage_a_T0_ref2$site) # 1580 SNPs
#rep3<-intersect(cage_d_T0_ref$site, cage_d_c1_ref$site) # 904
#ctrl_overlap <- intersect(intersect(rep1,rep2),rep3) #488 SNPs overlap all cages
#ctrl_overlap_2reps_A <- intersect(rep1,rep2) # 1050
#ctrl_overlap_2reps_A1D <- intersect(rep1,rep3) # 604
#ctrl_overlap_2reps_A2D <- intersect(rep2,rep3) # 608


# Make a unique table with all CTRL replicates - use this table to downsample reads
#ref_ctrl1<- merge(cage_a_c1_ref, cage_a_c2_ref, by="site", all=F)
#ref_ctrl2<- merge(cage_a_T0_ref1, cage_a_T0_ref2, by="site", all=F)
#ref_ctrl3<- merge(cage_d_T0_ref, cage_d_c1_ref, by="site", all=F)
#allcages_ref_ctrl_2reps <- merge(ref_ctrl1, ref_ctrl2, by="site", all=F)
#allcages_ref_ctrl_3reps <- merge(allcages_ref_ctrl_2reps, ref_ctrl3, by="site", all=F)
#print (c("dim overalp 3 replicates", dim(allcages_ref_ctrl_3reps))) #488 SNPs, 3280 Indv
#print (c("dim overalp 2 replicates", dim(allcages_ref_ctrl_2reps))) #1050 NSPs, 1858 indv

#alt_ctrl1<- merge(cage_a_c1_alt, cage_a_c2_alt, by="site", all=F)
#alt_ctrl2<- merge(cage_a_T0_alt1, cage_a_T0_alt2, by="site", all=F)
#alt_ctrl3<- merge(cage_d_T0_alt, cage_d_c1_alt, by="site", all=F)
#allcages_alt_ctrl_2reps <- merge(alt_ctrl1, alt_ctrl2, by="site", all=F)
#allcages_alt_ctrl_3reps <- merge(allcages_alt_ctrl_2reps, alt_ctrl3, by="site", all=F)
#dim(allcages_alt_ctrl_3reps) #488 SNPs, 3280 Indv
#dim(allcages_alt_ctrl_2reps) #1050 NSPs, 1858 indv

#write.table(allcages_alt_ctrl_2reps, "ctrl_alt_2reps_chr2R.txt", quote=F, sep="\t", col.names = T, row.names = F)
#write.table(allcages_alt_ctrl_3reps, "ctrl_alt_3reps_chr2R.txt", quote=F, sep="\t", col.names = T, row.names = F)
#write.table(allcages_ref_ctrl_2reps, "ctrl_ref_2reps_chr2R.txt", quote=F, sep="\t", col.names = T, row.names = F)
#write.table(allcages_ref_ctrl_3reps, "ctrl_ref_3reps_chr2R.txt", quote=F, sep="\t", col.names = T, row.names = F)

# get order of individuals in the allcages files
# at the end of HS

# HIGH SUGAR cages
#rep1hs<-intersect(cage_a_hs1_ref$site,cage_a_T0_ref1$site) #1723 SNPs
#rep2hs<- intersect(cage_a_hs2_ref$site,cage_a_T0_ref2$site) # 1945 SNPs
#rep3hs<-intersect(cage_d_T0_ref$site, cage_d_hs1_ref$site) # 847
#hs_overlap <- intersect(intersect(rep1hs,rep2hs),rep3hs) #561 SNPs overlap all cages
#hs_overlap_2reps_A <- intersect(rep1hs,rep2hs) # 1244
#hs_overlap_2reps_A1D <- intersect(rep1hs,rep3hs) # 647
#hs_overlap_2reps_A2D <- intersect(rep2hs,rep3hs) # 648

# Make a unique table with all HS replicates - use this table to downsample reads
#ref_hs1<- merge(cage_a_hs1_ref, cage_a_hs2_ref, by="site", all=F)
#ref_hs2<- merge(cage_a_T0_ref1, cage_a_T0_ref2, by="site", all=F)
#ref_hs3<- merge(cage_d_T0_ref, cage_d_hs1_ref, by="site", all=F)
#allcages_ref_hs_2reps <- merge(ref_hs1, ref_hs2, by="site", all=F)
#allcages_ref_hs_3reps <- merge(allcages_ref_hs_2reps, ref_hs3, by="site", all=F)
#print(c( "dim overalp 3reps HS ref", dim(allcages_ref_hs_3reps))) #561 SNPs, 3679 Indv
#print(c( "dim overalpt 2reps HS ref",dim(allcages_ref_hs_2reps))) #1244 NSPs, 2087 indv

#alt_hs1<- merge(cage_a_hs1_alt, cage_a_hs2_alt, by="site", all=F)
#alt_hs2<- merge(cage_a_T0_alt1, cage_a_T0_alt2, by="site", all=F)
#alt_hs3<- merge(cage_d_T0_alt, cage_d_hs1_alt, by="site", all=F)
#allcages_alt_hs_2reps <- merge(alt_hs1, alt_hs2, by="site", all=F)
#allcages_alt_hs_3reps <- merge(allcages_alt_hs_2reps, alt_hs3, by="site", all=F)
#print(c( "dim overalp 3reps HS alt", dim(allcages_alt_hs_3reps))) #561 SNPs, 3679 Indv
#print(c( "dim overalp 2reps HS alt", dim(allcages_alt_hs_2reps))) #1244 NSPs, 2087 indv

#write.table(allcages_alt_hs_2reps, "hs_alt_2reps_chr2R.txt", quote=F, sep="\t", col.names = T, row.names = F)
#write.table(allcages_alt_hs_3reps, "hs_alt_3reps_chr2R.txt", quote=F, sep="\t", col.names = T, row.names = F)
#write.table(allcages_ref_hs_2reps, "hs_ref_2reps_chr2R.txt", quote=F, sep="\t", col.names = T, row.names = F)
#write.table(allcages_ref_hs_3reps, "hs_ref_3reps_chr2R.txt", quote=F, sep="\t", col.names = T, row.names = F)

#order of individuals in the allcages files HS
flypercage <- unlist(n_filtered)  # number of individuals is x-1 because n_filtered contains hte columns "site"
firstflyC <- c(1, 1+cumsum(flypercage[1:5]-1))
names(firstflyC) <- names(flypercage[1:6])
lastflyC <- cumsum(flypercage[1:6]-1)
indC <-cbind.data.frame(firstflyC, lastflyC)
indC


firstflyHS <- c(1, 1+cumsum(flypercage[c(7,8,3,4,5)]-1)) # fixed this, it was 7,3,8,4,5, so the order was hs1 t01 hs2 t02, instead of hs1 hs2 t01 t02 as 
#in control tables. This doesnt matter much because teh samples per group are pretty well equilibrated. but now it is actually correct. the names in the 
#indpercate_HS tables is fixed in teh "Clean code" folder, the original files are left as they were
names(firstflyHS) <- names(flypercage[c(7,8,3,4,5,9)])
lastflyHS <- cumsum(flypercage[c(7,8,3,4,5,9)]-1)
indHS <-cbind.data.frame(firstflyHS, lastflyHS)
indHS

write.table(indC, "indpercage_CTRL_chr2R.txt", col.names = T,sep="\t", quote=F,row.names = T) 
write.table(indHS, "indpercage_HS_chr2R.txt", col.names = T,sep="\t",quote=F, row.names = T) 


# aHS1 2-525 , aHS2 526-1048, aT01 1049-1570, aT02 1571-2087, dT0 2088-2986, dHS1 2987-3679




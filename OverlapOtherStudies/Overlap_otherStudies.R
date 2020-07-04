##
## Longevity project
## Overlap of age-related SNPs in our data with previous studies
## 

#install.packages("SuperExactTest")
library(SuperExactTest)

setwd("C:/Users/lfpal/Documents/ViabilityExperiment/CleanCode/OverlapOtherStudies/")

#Load genes from previous studies and ours
#most studies extend the lengt of the gene to include potential regulatory regions
#our genes were extended by 1Kb as in several previous studies
  dat <- read.table("AgeRelatedGenes_previouspublications.txt",h=T)
  names <- colnames(dat)

#Make them into a list to be read by the SuperExactTest package
  datlist <- as.list(dat)
  datlist$ours_all_1kb <- c(as.character(dat$ours_HS_C_1kb),as.character( dat$ours_HS_1kb),as.character(dat$ours_C_1kb))
  datlist$ours_all_1kb <- unique(datlist$ours_all_1kb)
  datlist <- lapply(datlist, function(x) subset(x, !is.na(x)))  
#Lenght of each dataset inside the list
  length.gene.sets=sapply(datlist,length)
 
#Background size for significance estimation
  #Previous studies either do not report the background set,
  #or use a somehow arbitrary number (12k) without specifying how many 
  #genes were actually tested in their study. 
  total=15700  #number of genes tagged by SNPs in our study
  
#Testing all possible pairwise comparisons
  res=supertest(datlist, n=total, degree=2)  
  write.table(summary(res)$Table, file="SuperExactTest_output.txt", row.names=FALSE,
              sep="\t", quote=F, col.names = T)
  
#####
  
#Load supertest reslts, only including the pairwise comparisons with our results
  sub <- read.table("SuperExactTest_output_subset.txt",h=T)

#Adjust pvalues for multiple testing
  #sub$p.adj <- p.adjust(sub$P.value, method="BH" )
  #write.table(sub, "SuperExactTest_output_subset.txt", col.names = T, quote=F, sep="\t")
  
  
##### 
#Plots
  
#Overlap with Fecundity genes (from Durham et al)- using all our genes
#The results are the same when dividing genes into diet specific and boht-diets groups, show that inthe supp. tables
  
  fec <- c("durham_fecundityw1","durham_fecundityw3","durham_fecundityw5","durham_fecundityw7")
  b<-barplot(sub$Observed.Overlap[sub$group1=="ours_all_1kb" & sub$group2 %in% fec],ylim=c(0,75),
      names.arg = c("week 1", "week 3", "week 5", "wee 7") , ylab = "Overlap Fecundity genes & Longevity genes", xlab="Fecundity at Week")
  text (b , 2+sub$Observed.Overlap[sub$group1=="ours_all_1kb" & sub$group2 %in% fec],
        c(signif(sub$p.adj[sub$group1=="ours_all_1kb" & sub$group2=="durham_fecundityw1"], digits = 2),
          signif(sub$p.adj[sub$group1=="ours_all_1kb" & sub$group2=="durham_fecundityw3"], digits = 2),
          signif(sub$p.adj[sub$group1=="ours_all_1kb" & sub$group2=="durham_fecundityw5"], digits = 2),
          signif(sub$p.adj[sub$group1=="ours_all_1kb" & sub$group2=="durham_fecundityw7"], digits = 2)))
  text (b , 4+sub$Observed.Overlap[sub$group1=="ours_all_1kb" & sub$group2 %in% fec],
        c("" , "", round(sub$FE[sub$group1=="ours_all_1kb" & sub$group2=="durham_fecundityw5"],digits = 2),
          round(sub$FE[sub$group1=="ours_all_1kb" & sub$group2=="durham_fecundityw7"], digits=2)))
  
  
#Overlap with diet and diet*longevity genes (from Hoedjes et al)
  #There are differences depending on the group, so plot the 3 groups
  
  diet <- "hoedjes_diet"
  dietlife <- "hoedjes_diet_lifespan"
  
  alldiet <- cbind(c(sub$Observed.Overlap[sub$group1=="ours_C_1kb" & sub$group2 %in% diet],
               sub$Observed.Overlap[sub$group1=="ours_HS_1kb" & sub$group2 %in% diet],
               sub$Observed.Overlap[sub$group1=="ours_HS_C_1kb" & sub$group2 %in% diet]),
               c(sub$Observed.Overlap[sub$group1=="ours_C_1kb" & sub$group2 %in% dietlife],
                 sub$Observed.Overlap[sub$group1=="ours_HS_1kb" & sub$group2 %in% dietlife],
                 sub$Observed.Overlap[sub$group1=="ours_HS_C_1kb" & sub$group2 %in% dietlife]))
  alldiet.p <- cbind(c(sub$p.adj[sub$group1=="ours_C_1kb" & sub$group2 %in% diet],
                     sub$p.adj[sub$group1=="ours_HS_1kb" & sub$group2 %in% diet],
                     sub$p.adj[sub$group1=="ours_HS_C_1kb" & sub$group2 %in% diet]),
                   c(sub$p.adj[sub$group1=="ours_C_1kb" & sub$group2 %in% dietlife],
                     sub$p.adj[sub$group1=="ours_HS_1kb" & sub$group2 %in% dietlife],
                     sub$p.adj[sub$group1=="ours_HS_C_1kb" & sub$group2 %in% dietlife]))
  alldiet.e <- cbind(c(sub$FE[sub$group1=="ours_C_1kb" & sub$group2 %in% diet],
                       sub$FE[sub$group1=="ours_HS_1kb" & sub$group2 %in% diet],
                       sub$FE[sub$group1=="ours_HS_C_1kb" & sub$group2 %in% diet]),
                     c(sub$FE[sub$group1=="ours_C_1kb" & sub$group2 %in% dietlife],
                       sub$FE[sub$group1=="ours_HS_1kb" & sub$group2 %in% dietlife],
                       sub$FE[sub$group1=="ours_HS_C_1kb" & sub$group2 %in% dietlife]))
  
  
 
  d<-barplot(alldiet, beside = T , ylim=c(0,50),
             ylab = "Overlap Diet genes & Longevity genes", xlab="", legend.text = c("C_only", "HS_only", "boht_diets"),
             args.legend=list(bty="n", x="topleft"))
  axis(side = 1, at=c(2.5,6.5), labels=c("Diet genes", "Diet-by-lifespan"), cex=0.8)
  text (d , 2+alldiet,
        signif(alldiet.p, digits = 2))
  text (d , 4+alldiet, 
        signif(alldiet.e, digits = 2))
        
  
  
#Overlap with other longevity studies, including selection experimetns and standing variation
  selection <- c("fabian", "carnes","remolina", "hoedjes_lifespan")
  standing <- c("durham_lifespan", "huang_DGRP", "huang_AIP")
  mutant <-  c("genAge_dmel")
  ssm <- c("fabian", "carnes","remolina", "hoedjes_lifespan", "durham_lifspan", "huang_DGRP", "huang_AIP","genAge_dmel")
  ssm2 <- c("Fabian", "Carnes","Remolina", "Hoedjes", "Durham", "Huang_DGRP", "Huang_AIP","GenAge")
  
  alllifespan <- sub$Observed.Overlap[sub$group1=="ours_all_1kb" & sub$group2 %in% ssm]
  alllifespan.p <- sub$p.adj[sub$group1=="ours_all_1kb" & sub$group2 %in% ssm]
  alllifespan.e <- sub$FE[sub$group1=="ours_all_1kb" & sub$group2 %in% ssm]
  
  l<-barplot(alllifespan[c(8,6,7,2,5,4,3,1)], ylim=c(0,355),
              ylab = "Overlap with previous studies (# genes)", xlab="Previous studies on lifespan",
             col=c(rep("cadetblue",4), rep("salmon",3), "orange"),
             legend=c("Selection experiments", "Standing variation", "Mutational experiments"),
             args.legend=list(bty="n", fill=c("cadetblue","salmon","orange"), x="topleft"),
                        names.arg = ssm2)
  #show pvalue info
    # text (l , 5+alllifespan[c(8,6,7,2,5,4,3,1)],
        signif(alllifespan.p[c(8,6,7,2,5,4,3,1)], digits = 2))
  #show enrichment info
        text (l , 15+alllifespan[c(8,6,7,2,5,4,3,1)], 
        signif(alllifespan.e[c(8,6,7,2,5,4,3,1)], digits = 2))
  #show significance of the enrichment
    text (l+0.3 , 15+alllifespan[c(8,6,7,2,5,4,3,1)], 
          c("***","***","ns", "***", "ns", "***", "***", "ns" ))
 
    
  
    
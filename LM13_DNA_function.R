##### Functions for Denef et al, submitted to Frontiers in Microbiology ########

##### Normalization #######

# Better rounding function than R's base round
matround <- function(x){trunc(x+0.5)}

# Scales reads by 
# 1) taking proportions
# 2) multiplying by a given library size of n
# 3) rounding 
# Default for n is the minimum sample size in your library
# Default for round is floor
scale_reads <- function(physeq, n = min(sample_sums(physeq)), round = "floor") {
  
  # transform counts to n
  physeq.scale <- transform_sample_counts(physeq, 
                                          function(x) {(n * x/sum(x))}
  )
  
  # Pick the rounding functions
  if (round == "floor"){
    otu_table(physeq.scale) <- floor(otu_table(physeq.scale))
  } else if (round == "round"){
    otu_table(physeq.scale) <- matround(otu_table(physeq.scale))
  }
  
  # Prune taxa and return new phyloseq object
  physeq.scale <- prune_taxa(taxa_sums(physeq.scale) > 0, physeq.scale)
  return(physeq.scale)
}


######## ADONIS ###########

# Function to run adonis test on a phyloseq object and a variable from metadata
# Make sure OTU data is standardized/normalized before 
phyloseq_to_adonis <- function(physeq, distmat = NULL, dist = "bray", formula) {
  
  if(!is.null(distmat)){
    phydist <- distmat
  } else {
    phydist <- phyloseq::distance(physeq, dist)
  }
  
  metadata <- as(sample_data(physeq), "data.frame")
  
  # Adonis test
  f <- reformulate(formula, response = "phydist")
  adonis.test <- adonis(f, data = metadata)
  print(adonis.test)
  
  # Run homogeneity of dispersion test if there is only 1 variable
  if (length(formula) == 1) {
    
    group <- metadata[,formula]
    beta <- betadisper(phydist, group)
    disper.test = permutest(beta)
    print(disper.test)
    
    l <- list(
      dist = phydist, 
      formula = f, 
      adonis = adonis.test, 
      disper = disper.test
    )
    
  } else {
    
    l <- list(
      dist = phydist, 
      formula = f, 
      adonis = adonis.test
    )
  }
  return (l)
}

###### Merge functions ############

# Merge samples by averaging OTU counts instead of summing
merge_samples_mean <- function(physeq, group, round){
  # Calculate the number of samples in each group
  group_sums <- as.matrix(table(sample_data(physeq)[ ,group]))[,1]
  
  # Merge samples by summing
  merged <- merge_samples(physeq, group)
  
  # Divide summed OTU counts by number of samples in each group to get mean
  # Calculation is done while taxa are columns
  x <- as.matrix(otu_table(merged))
  if(taxa_are_rows(merged)){ x<-t(x) }
  
  # Pick the rounding functions
  if (round == "floor"){
    out <- floor(t(x/group_sums))
  } else if (round == "round"){
    out <- matround(t(x/group_sums))
  }
  
  # Return new phyloseq object with taxa as rows
  out <- otu_table(out, taxa_are_rows = TRUE)
  otu_table(merged) <- out
  return(merged)
}

# Merge samples, just including OTUs that were present in all merged samples
# Call this function before running merge_samples()
merge_OTU_intersect <- function(physeq, group){
  
  # Make sure we're not starting with more taxa than we need 
  physeq <- prune_taxa(taxa_sums(physeq) > 0, physeq)
  
  s <- data.frame(sample_data(physeq))
  l <- levels(s[,group])
  o <- otu_table(physeq)
  
  # Loop through each category
  for (cat in 1:length(l)) {
    
    # Get the index of all samples in that category
    w <- which(s[,group]==l[cat])
    
    # subset to just those columns of OTU table
    cat.sub<-o[,w]
    print(dim(cat.sub))
    
    # Find the indices of 0's in the OTU table
    zeros <- apply(cat.sub, 1, function(r) any(r == 0))
    
    # If an OTU had a 0 in at least one sample, change all samples to 0
    cat.sub[zeros,] <- 0
  }
  
  o[,w] <- cat.sub
  otu_table(physeq) <- o
  
  return(physeq)
  
}

####### DEseq ############
## Marain Schmidt wrote these 2 functions off of the following tutorial from the Phyloseq GitHub page:
#http://joey711.github.io/phyloseq-extensions/DESeq2.html

deSEQ <- function(data, valuetest){
  data_pruned=prune_taxa(taxa_sums(data)>(147*nrow(sample_data(data))),data)
  de_data = phyloseq_to_deseq2(data_pruned, valuetest)
  de_data2 = DESeq(de_data, test="Wald", fitType="parametric")
  res_data = results(de_data2, cooksCutoff = FALSE)
  alpha = 0.01
  sig_data = res_data[which(res_data$padj < alpha), ]
  sigtab_sherm = cbind(as(sig_data, "data.frame"), as(tax_table(data_pruned)[rownames(sig_data), ], "matrix"))
} 

deSEQ_noprune <- function(data, valuetest){
  data=prune_taxa(taxa_sums(data)>0,data)
  de_data = phyloseq_to_deseq2(data, valuetest)
  de_data2 = DESeq(de_data, test="Wald", fitType="parametric")
  res_data = results(de_data2, cooksCutoff = FALSE, contrast=c("DNA","cD","D"))
  plotMA(res_data)
  alpha = 1
  sig_data = res_data[which(res_data$padj < alpha), ]
  sigtab_sherm = cbind(as(sig_data, "data.frame"), as(tax_table(data)[rownames(sig_data), ], "matrix"))
} 

plot_deSEQ <- function(deSEQdata, title){
  y = tapply(deSEQdata$log2FoldChange, deSEQdata$Species, function(x) max(x))
  y = sort(y, TRUE)
  deSEQdata$Species = factor(as.character(deSEQdata$Species), levels=names(y))
  ggplot(deSEQdata, aes(x=Phylum, y=log2FoldChange, color=Phylum, size=plotvalue, shape=sig)) + 
    geom_point(alpha=0.9) + theme_bw() +  ggtitle(title) +  
    scale_color_manual(values = phylum.colors,name="Phylum") +
    #    scale_color_manual(name="p-value",  breaks = c("0", "1"), labels = c("p>1", "p<0.01"), values = c("0" = "grey", "1"="black")) +
    scale_shape_manual(name = "p-value", breaks = c("0", "1"), 
                       labels = c("p>0.01", "p<0.01"),
                       values = c("0" = 21, "1"= 19)) +
    scale_size_continuous("abundance") +
    scale_y_continuous(breaks=seq(-15,15,1),limits=c(-7,7)) + 
    theme(axis.title.x = element_text(face="bold", size=16),
          axis.text.x = element_text(angle=30, colour = "black", vjust=1, hjust = 1, size=12),
          axis.text.y = element_text(colour = "black", size=16),
          axis.title.y = element_text(face="bold", size=16),
          plot.title = element_text(face="bold", size = 20),
          legend.title = element_text(size=12, face="bold"),
          legend.text = element_text(size = 12),
          strip.background = element_rect(colour="black"),
          legend.position="right")
}

lot_deSEQ_combo <- function(deSEQdata, title){
  y = tapply(deSEQdata$log2FoldChange, deSEQdata$Species, function(x) max(x))
  y = sort(y, TRUE)
  deSEQdata$Species = factor(as.character(deSEQdata$Species), levels=names(y))
  ggplot(deSEQdata, aes(x=Phylum_plot, y=log2FoldChange, color=Phylum, size=plotvalue, shape=sig)) + 
    geom_point(alpha=0.9) + theme_bw() +  ggtitle(title) +  
    scale_color_manual(values = phylum.colors,name="Phylum") +
    #    scale_color_manual(name="p-value",  breaks = c("0", "1"), labels = c("p>0.05", "p<0.05"), values = c("0" = "grey", "1"="black")) +
    scale_shape_manual(name = "p-value", breaks = c("0", "1"), 
                       labels = c("p>0.05", "p<0.05"),
                       values = c("0" = 21, "1"= 19)) +
    scale_size_continuous("abundance") +
    scale_y_continuous(breaks=seq(-15,15,1),limits=c(-6,5)) + 
    theme(axis.title.x = element_text(face="bold", size=16),
          axis.text.x = element_text(angle=90, colour = "black", vjust=1, hjust = 1, size=12),
          axis.text.y = element_text(colour = "black", size=16),
          axis.title.y = element_text(face="bold", size=16),
          plot.title = element_text(face="bold", size = 20),
          legend.title = element_text(size=12, face="bold"),
          legend.text = element_text(size = 12),
          strip.background = element_rect(colour="black"),
          legend.position="right")
}





### functions used in core analysis


prevalence <- function(physeq, add_tax = TRUE){
  
  ## Check if taxa are rows
  trows <- taxa_are_rows(physeq)
  
  ## Extract OTU table
  otutab <- as.data.frame(otu_table(physeq))
  
  ## Transpose OTU table (species should be arranged by rows)
  if(trows == FALSE){
    otutab <- t(otutab)
  }
  
  ## Estimate prevalence (number of samples with OTU present)
  prevdf <- apply(X = otutab,
                  #MARGIN = ifelse(trows, yes = 1, no = 2),  # for a non-transposed data
                  MARGIN = 1,
                  FUN = function(x){sum(x > 0)})
  
  ## Add total and average read counts per OTU
  prevdf <- data.frame(Prevalence = prevdf,
                       TotalAbundance = taxa_sums(physeq),
                       MeanAbundance = rowMeans(otutab),
                       MedianAbundance = apply(otutab, 1, median))
  
  ## Add taxonomy table
  if(add_tax == TRUE && !is.null(tax_table(physeq, errorIfNULL = F))){
    prevdf <- cbind(prevdf, tax_table(physeq))
  }
  return(prevdf)
}




###phyloseq filter prevalence function code (from metagMisc)

phyloseq_filter_prevalence <- function(physeq, prev.trh = 0.05, abund.trh = NULL, threshold_condition = "OR", abund.type = "total"){
  
  ## Threshold validation
  if(prev.trh > 1 | prev.trh < 0){ stop("Prevalence threshold should be non-negative value in the range of [0, 1].\n") }
  if(!is.null(abund.trh)){ 
    if(abund.trh <= 0){ stop("Abundance threshold should be non-negative value larger 0.\n") }
  }
  
  ## Check for the low-prevalence species (compute the total and average prevalences of the features in each phylum)
  prevdf_smr <- function(prevdf){
    ddply(prevdf, "Phylum", function(df1){ data.frame(Average = mean(df1$Prevalence), Total = sum(df1$Prevalence))})
  }
  # prevdf_smr( prevalence(physeq) )
  
  ## Check the prevalence threshold
  # phyloseq_prevalence_plot(prevdf, physeq)
  
  ## Define prevalence threshold as % of total samples
  ## This function is located in 'phyloseq_prevalence_plot.R' file
  prevalenceThreshold <- prev.trh * phyloseq::nsamples(physeq)
  
  ## Calculate prevalence (number of samples with OTU) and OTU total abundance
  prevdf <- prevalence(physeq)
  
  ## Get the abundance type
  if(abund.type == "total") { prevdf$AbundFilt <- prevdf$TotalAbundance }
  if(abund.type == "mean")  { prevdf$AbundFilt <- prevdf$MeanAbundance }
  if(abund.type == "median"){ prevdf$AbundFilt <- prevdf$MedianAbundance }
  
  ## Which taxa to preserve
  if(is.null(abund.trh)) { tt <- prevdf$Prevalence >= prevalenceThreshold }
  if(!is.null(abund.trh)){
    ## Keep OTU if it either occurs in many samples OR it has high abundance
    if(threshold_condition == "OR"){
      tt <- (prevdf$Prevalence >= prevalenceThreshold | prevdf$AbundFilt >= abund.trh)
    }
    
    ## Keep OTU if it occurs in many samples AND it has high abundance
    if(threshold_condition == "AND"){
      tt <- (prevdf$Prevalence >= prevalenceThreshold & prevdf$AbundFilt >= abund.trh)
    }
  }
  
  ## Extract names for the taxa we whant to keep
  keepTaxa <- rownames(prevdf)[tt]
  
  ## Execute prevalence filter
  res <- phyloseq::prune_taxa(keepTaxa, physeq)
  return(res)
}




##create function that estimates abundance weighted phylogenetic diversity (modified from estimate_pd() function)


estimate_bwpd <- function(phylo){
  
  # Error if input is not of class phylo
  if(class(phylo) != "phyloseq"){
    stop("Input file is not of class 'phyloseq'.")
  }
  
  # Error if no class phy_tree
  if(!(.hasSlot(phylo, "phy_tree"))){
    stop("Could not find tree slot in phylo object.")
  }
  
  # Transpose if needed
  # Adapted from phyloseq/vegan import
  OTU <- phyloseq::otu_table(phylo)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  
  # Get matrix version of OTU table
  otutable <- as(OTU, "matrix")
  
  # Get phylogenetic tree from pyloseq object
  tree <- phyloseq::phy_tree(phylo)
  
  # Print status message
  message("Calculating Faiths PD-index...")
  
  # If object is greater than 10mb, then print status message
  if(object.size(otutable) > 10000000){
    message("This is a large object, it may take awhile...")
  }
  
  # Calculate Faith's PD-index
  pdtable <- picante::pse(otutable, tree)
  
  # Return data frame of results
  return(pdtable)
}


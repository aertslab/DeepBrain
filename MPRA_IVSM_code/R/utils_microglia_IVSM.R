process.enh.cDNA.plasmid <- function(env, 
                                     line, 
                                     cDNA.counts.file.path, 
                                     plasmid.counts.file.path, 
                                     min.cDNA.bcs = 5,
                                     min.plasmid.bcs = 5,
                                     min.counts = 2,
                                     keep.shuffled = T,
                                     gauss.est = "median", # or "max.dens"
                                     loess.norm = F,
                                     out.dir = NULL) {
  ########### Load packages ###########
  library(dplyr)
  library(edgeR)
  library(preprocessCore)
  library(affy)
  library(ggplot2)
  library(ggpubr)
  library(reshape2)
  library(stringr)
  
  ########### Load, aggregate and CPM normalize data ###########
  # Load data  
  message("Processing ", line, " with filtering set to ", min.cDNA.bcs," barcodes for cDNA and ", min.plasmid.bcs," barcodes for plasmid.\n")
  env[[line]] <- new.env()
  env[[line]]$plasmid.raw <- read.table(file = plasmid.counts.file.path, header = F, sep = "\t", stringsAsFactors = F)
  env[[line]]$cDNA.raw <- read.table(file = cDNA.counts.file.path, header = F, sep = "\t", stringsAsFactors = F)

  # Aggregate the cDNA counts by enhancer: per enhancer add all the counts for all coupled rBCs together, 
  # and only keep enhancers is they have a min number of rBC linked, and store in plasmid.aggregated
  env[[line]]$plasmid.raw %>%
    group_by(V1) %>% 
    dplyr::summarise(
      counts = sum(V3),
      nb_barcodes = length(V1)
    ) %>%
    filter(nb_barcodes >= min.plasmid.bcs & counts >= min.counts) %>%
    magrittr::set_colnames(c("synthetic.region", "counts", "nb_barcodes")) -> env[[line]]$plasmid.aggregated
  env[[line]]$cDNA.raw %>% 
    group_by(V1) %>% 
    dplyr::summarise(
      counts = sum(V3),
      nb_barcodes = length(V1)
    ) %>%
    filter(nb_barcodes >= min.cDNA.bcs & counts >= min.counts) %>%
    magrittr::set_colnames(c("synthetic.region", "counts", "nb_barcodes")) -> env[[line]]$cDNA.aggregated

  # CPM Normalization (normalise using the number counts in this final counts matrix)
  message("CPM normalizing...")
  env[[line]]$plasmid.aggregated$CPM.counts <- edgeR::cpm(y = env[[line]]$plasmid.aggregated[,"counts"])[,1]
  env[[line]]$cDNA.aggregated$CPM.counts <- edgeR::cpm(y = env[[line]]$cDNA.aggregated[,"counts"])[,1]
  message("[OK]\n")

  ########### Merge data then Loess, Input and Basal normalization ###########
  # Merge plasmid and cDNA counts
  message("Merging plasmid and cDNA counts...")
  env[[line]]$merged <- inner_join(env[[line]]$plasmid.aggregated, 
                                   env[[line]]$cDNA.aggregated, 
                                   by = "synthetic.region", 
                                   suffix = c(".plasmid", ".cDNA"))
  
  # Assign enhancer class
  # First assign for all enhancers 'NA'
  # Then grep specific characters and put a global name for the groups of enhancers
  env[[line]]$merged$name <- str_sub(string = env[[line]]$merged$synthetic.region, 1, -14)
  env[[line]]$merged$class <- "NA"
  env[[line]]$merged[grep(pattern = "FIRE-ST",env[[line]]$merged$name),"class"] <- "FIRE_SM"
  env[[line]]$merged[grep(pattern = "Shuffle_",env[[line]]$merged$name),"class"] <- "Shuffled"

  message("[OK]\n")
  message("Number of enhancers remaining after combining cDNA and Plasmid: ", nrow(env[[line]]$merged))
  message("cDNA reads:    Median = ", median(x = env[[line]]$merged$counts.cDNA),"  Mean = ", mean(x = env[[line]]$merged$counts.cDNA))
  message("Plasmid reads: Median = ", median(x = env[[line]]$merged$counts.plasmid),"  Mean = ", mean(x = env[[line]]$merged$counts.plasmid), "\n")
  message("Normalizing activity...")

  # loess Normalization
  loess.matrix <- normalize.loess(as.matrix(env[[line]]$merged[,c("CPM.counts.plasmid","CPM.counts.cDNA")]))
  env[[line]]$merged$loess.CPM.counts.plasmid <- loess.matrix[,1]
  env[[line]]$merged$loess.CPM.counts.cDNA <- loess.matrix[,2]
  
  # Input normalization
  env[[line]]$merged.input.norm <- data.frame("synthetic.region" = env[[line]]$merged$synthetic.region,
                                              "class" = env[[line]]$merged$class,
                                              stringsAsFactors = F)
  if (loess.norm) {
    ## Use loess-CPM normalized data
    env[[line]]$merged.input.norm[["CPM.InputNorm"]] <- unlist(env[[line]]$merged[,"loess.CPM.counts.cDNA"]/env[[line]]$merged[,"loess.CPM.counts.plasmid"])
  } else {
    ## Use CPM normalized data
    env[[line]]$merged.input.norm[["CPM.InputNorm"]] <- unlist(env[[line]]$merged[,"CPM.counts.cDNA"]/env[[line]]$merged[,"CPM.counts.plasmid"])
  }
  ## Add log2 value
  env[[line]]$merged.input.norm$log.CPM.InputNorm <- log2(env[[line]]$merged.input.norm$CPM.InputNorm)
  
  # Basal expression normalization
  shuffled <- env[[line]]$merged.input.norm[grep(pattern = "Shuffle_", x = env[[line]]$merged.input.norm$synthetic.region),]
  basal.ex <- median(shuffled[,"CPM.InputNorm"])
  env[[line]]$merged.input.basal.norm <- data.frame("synthetic.region" = env[[line]]$merged.input.norm$synthetic.region,
                                                    "class" = env[[line]]$merged.input.norm$class,
                                                    stringsAsFactors = F)
  env[[line]]$merged.input.basal.norm$CPM.Input.BasalNorm <- env[[line]]$merged.input.norm[,"CPM.InputNorm"]/basal.ex
  ## Add log2 value
  env[[line]]$merged.input.basal.norm$log.CPM.Input.BasalNorm <- log2(env[[line]]$merged.input.basal.norm$CPM.Input.BasalNorm)
  
  # Remove shuffled tiles
  if (keep.shuffled == F) {
    env[[line]]$merged.input.basal.norm <- env[[line]]$merged.input.basal.norm[!(env[[line]]$merged.input.basal.norm$class == "Shuffled"),]
  }
  message("[OK]\n")
  
  ########### Calculate p-value and adjusted p-value of enhancer activity ###########
  # robust fit of a gaussian by robust estimates of it's two params
  message("Estimation of active enhancers...")
  sample <- env[[line]]$merged.input.basal.norm
  subsettype <- c("Shuffled")
  print(paste0("Number of sequences used for curve fitting: ", nrow(sample[sample$class %in% subsettype,])))
  if (gauss.est == "median"){
    ## Use median
    location <- median(x = log2(sample[sample$class %in% subsettype, "CPM.Input.BasalNorm"]))
    print(paste0("Median: ",location))
  } else if (gauss.est == "max.dens"){
    ## Use max density value
    location <- density(log2(sample[sample$class %in% subsettype, "CPM.Input.BasalNorm"]))$x[which.max(density(log2(sample[sample$class %in% subsettype, "CPM.Input.BasalNorm"]))$y)]
    print(paste0("Max density: ",location))
  } else {
    stop("Incorrect value for 'gauss.est'")
  }
  scale <- mad(x = log2(sample[sample$class %in% subsettype, "CPM.Input.BasalNorm"]))
  
  # Determine p-values and adjussted p-values
  env[[line]]$merged.input.basal.norm$pvalue <- pnorm(q = sample$log.CPM.Input.BasalNorm, mean = location, sd = scale, lower.tail = F)
  env[[line]]$merged.input.basal.norm$padj <- p.adjust(p = env[[line]]$merged.input.basal.norm$pvalue, method = "fdr")
  message("[OK]\n")
  
  # controlling the FDR on this should result in < 5% false discoveries
  print("Repartition of active enhancers:")
  tbl_act <- table(env[[line]]$merged.input.basal.norm[env[[line]]$merged.input.basal.norm$padj < 0.05, "class"])
  print(tbl_act/sum(env[[line]]$merged.input.basal.norm$padj < 0.05)*100)
  print("Proportion of active enhancers per class:")
  tbl_class <- table(env[[line]]$merged.input.basal.norm[, "class"])
  print(tbl_act/tbl_class[names(tbl_class) %in% names(tbl_act)]*100)
  print(paste0("Number of enhancers with padj < 0.05: ",nrow(env[[line]]$merged.input.basal.norm[env[[line]]$merged.input.basal.norm$padj < 0.05,])))
  print(paste0("Number of Negative control sequences with padj < 0.05: ",nrow(env[[line]]$merged.input.basal.norm[env[[line]]$merged.input.basal.norm$padj < 0.05 & env[[line]]$merged.input.basal.norm$class == "Shuffled",])))

  
  ########### Combined replicates ###########
  message("Combining replicates...")
  if (!("Combined" %in% names(env))){
    # Create combined dataset
    env[["Combined"]] <- new.env()
    env[["Combined"]]$Input.normalised <- data.frame("synthetic.region" = env[[line]]$merged$synthetic.region,
                                                     "class" = env[[line]]$merged$class, 
                                                     stringsAsFactors = F)
    env[["Combined"]]$Basal.normalised <- env[["Combined"]]$Input.normalised 
  }
  if (paste0("CPM.InputNorm.",line) %in% names(env[["Combined"]]$Input.normalised)){
    # Remove data if sample was previously processed
    env[["Combined"]]$Input.normalised <- env[["Combined"]]$Input.normalised[, !names(env[["Combined"]]$Input.normalised) %in% paste0("CPM.InputNorm.",line)]
    env[["Combined"]]$Basal.normalised <- env[["Combined"]]$Basal.normalised[, !names(env[["Combined"]]$Basal.normalised) %in% paste0("CPM.InputNorm.",line)]
  }
  # Save sample data in combined data frames
  env[["Combined"]]$Input.normalised <- full_join(env[["Combined"]]$Input.normalised, 
                                                  env[[line]]$merged.input.norm[, c("synthetic.region", "CPM.InputNorm")], 
                                                  by = "synthetic.region")
  names(env[["Combined"]]$Input.normalised)[length(names(env[["Combined"]]$Input.normalised))] <- paste0("CPM.InputNorm.",line)
  env[["Combined"]]$Basal.normalised <- full_join(env[["Combined"]]$Basal.normalised, 
                                                  env[[line]]$merged.input.basal.norm[, c("synthetic.region", "CPM.Input.BasalNorm")], 
                                                  by = "synthetic.region")
  names(env[["Combined"]]$Basal.normalised)[length(names(env[["Combined"]]$Basal.normalised))] <- paste0("CPM.InputNorm.",line)
  
  # Reassign class
  env[["Combined"]]$Basal.normalised$name <- str_sub(string = env[["Combined"]]$Basal.normalised$synthetic.region, 1, -14)
  env[["Combined"]]$Basal.normalised$class <- "NA"
  env[["Combined"]]$Basal.normalised[grep(pattern = "FIRE-ST",env[["Combined"]]$Basal.normalised$name),"class"] <- "FIRE_SM"
  env[["Combined"]]$Basal.normalised[grep(pattern = "Shuffle_",env[["Combined"]]$Basal.normalised$name),"class"] <- "Shuffled"
  message("[OK]\n")

  ########### QC ###########
  # Histogram plasmid barcode distribution
  # hist(x = env[[line]]$merged$nb_barcodes.plasmid, breaks = 200)
  # abline(v = min.plasmid.bcs, col = "red")
  
  # Histogram cDNA barcode distribution
  # hist(x = env[[line]]$merged$nb_barcodes.cDNA, breaks = 200)
  # abline(v = min.cDNA.bcs, col = "red")
  
  # Plot barcode and read distribution
  plot1 <- ggplot(melt(env[[line]]$merged[,!names(env[[line]]$merged) %in% c("counts.plasmid","class","counts.cDNA")]), aes(x = log2(value), fill = variable)) +
    geom_density(alpha = .3) +
    labs(title = paste0("CHEQ-seq ", line," - Plasmid / cDNA barcode distribution"), subtitle = "Log CPM Normalized", y = "Density", x = "Log2 Counts")
  plot2 <- ggplot(env[[line]]$merged, aes(x = counts.plasmid, y = counts.cDNA)) +
    geom_point()
  plot3 <- ggplot(env[[line]]$merged, aes(x = log(counts.plasmid), y = log(counts.cDNA))) +
    geom_point()
  p <- ggarrange(plot1, ggarrange(plot2, plot3,
                                  labels = c("B", "C"),
                                  ncol = 2),
                 labels = c("A"),
                 nrow = 2)
  print(p)
  
  # MA plots
  ## No normalization
  matrix <- as.matrix(env[[line]]$merged[,c("CPM.counts.plasmid","CPM.counts.cDNA")])
  ma.plot(A = rowMeans(log2(matrix)), 
          M = log2(matrix[,2])-log2(matrix[,1]), 
          cex = 1, pch = 16, )
  title(main = paste0("CHEQ-seq Synthetic Enhancer - ",line), sub = "CPM Normalized MA plot")
  ## Loess normalization
  if (loess.norm == F){
    matrix_norm <- normalize.loess(matrix)
  } else {
    matrix_norm <- as.matrix(env[[line]]$merged[,c("loess.CPM.counts.plasmid","loess.CPM.counts.cDNA")])
  }
  ma.plot(A = rowMeans(log2(matrix_norm)),
          M = log2(matrix_norm[,2])-log2(matrix_norm[,1]),
          cex = 1, pch = 16)
  title(main = paste0("CHEQ-seq Synthetic Enhancer - ",line), sub = "loess CPM Normalized MA plot")
  
  # Violin plots - expression per class (input basal norm)
  plot5 <- ggplot(data = env[[line]]$merged.input.basal.norm, aes(x = class, y = log2(CPM.Input.BasalNorm), fill = class)) +
    geom_violin() +
    geom_jitter(data = env[[line]]$merged.input.basal.norm[env[[line]]$merged.input.basal.norm$padj < 0.05,], shape = 16, position = position_jitter(0.2), cex = 0.8) +
    geom_hline(yintercept = log2(1), linetype = "dashed", color = "red", size = 0.5) +
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_text(size=16, face = "bold"),
          legend.position = "none",
          axis.text.x = element_text(angle = 60, hjust = 1),
          axis.text = element_text(size = 8, face = "bold"),
          axis.title = element_text(size = 14, face = "bold")) +
    labs(x = "Class", y = "Log2 FC Input Basal Normalized", title = paste0(environmentName(env), " - ", line))
  print(plot5)
  
  # check fit gaussian
  plot6 <- ggplot(data = env[[line]]$merged.input.basal.norm, mapping = aes(x = log2(CPM.Input.BasalNorm), fill = class)) +
    geom_density(size = 0, alpha = .4) +  
    stat_function(fun = dnorm, args = list(mean = location, sd = scale)) + 
    theme_minimal() +
    labs(fill = "Class", x = "Log2 FC Input Basal Normalized", y = "Density")
  print(plot6)
  
  # Bar plot active enhancers
  plot7 <- ggplot(data = env[[line]]$merged.input.basal.norm, mapping = aes(x = padj < 0.05, fill = class)) + 
    geom_bar(alpha = .8) +
    theme_minimal() +
    labs(fill = "Class")
  print(plot7)
  
  ########### Adding metadata ###########
  env[[line]]$metadata <- list()
  env[[line]]$metadata$Measured.enhancers <- nrow(env[[line]]$merged)
  env[[line]]$metadata$Basal.expression <- basal.ex
  env[[line]]$metadata$BC.cDNA.filtering <- min.cDNA.bcs
  env[[line]]$metadata$BC.plasmid.filtering <- min.plasmid.bcs  

  
  ########### Saving ###########
  # Save (you can specify an outdir and it will write a table).
  # Two outputs: the first one with the log2 of the CPM input-basal-normalised values; second with the log2 of the CPM input-normalised values.
  if (!is.null(out.dir)){
    write.table(x = env[[line]]$merged.input.basal.norm[,c(1,4)],
                file = file.path(out.dir, paste0(line,"_",environmentName(env),"__filter_",min.cDNA.bcs,"-",min.plasmid.bcs,"__cpm_input_basal_normalized.tsv")),
                quote = F, sep = "\t", row.names = F, col.names = T)
    write.table(x = env[[line]]$merged.input.norm[,c(1,4)],
                file = file.path(out.dir, paste0(line,"_",environmentName(env),"__filter_",min.cDNA.bcs,"-",min.plasmid.bcs,"__cpm_input_normalized.tsv")),
                quote = F, sep = "\t", row.names = F, col.names = T)
  }
  return (env)
} 


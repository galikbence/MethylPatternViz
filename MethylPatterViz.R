#!/usr/local/bin/R
#Load packages
options(warn=-1)
library(argparser)
library(ggplot2)
library(seqinr)
library(dplyr)
library(gridExtra)
library(openxlsx)

#Command line args
parser <- arg_parser(description='Process commandline arguments. Input files are Bismark coverage files (*_trimmed_bismark_bt2.bismark.cov).')
parser <- add_argument(parser, arg=c("--chr", "--start", "--end", "--gene", "--all_samples", "--sample", "--paired", "--path"), 
                       help = c("chromosome", "start site", "stop site", "plot/legen name or gene name", 
                                "process all samples in the folder T(rue) or F(alse)", "process one or more samples, separeated by comma",
                                "plot results as pairs, T(rue) or F(alse), only works with even number of samples", "path of the folder contains the chr sequences"),
                       type = c("numeric", "numeric", "numeric", "character", "character", "character", "character", "character"))

args = parse_args(parser)

#Set up input parameters and input files
if(args$all_samples == "T"){cov_files <- list.files(path = ".", pattern = "*.cov")} else {cov_files <- unlist(strsplit(args$sample, ","))}

chr <- args$chr
start <- args$start
end <- args$end
gene <- args$gene
g_path <- args$path 
  
#Make empty CpG pattern table
cpg_table <- matrix(data = c(start:end), byrow = F, ncol = 1)

#Progress message
cat("######################", sep = "\n")
cat(paste("Read in chromosome: ", chr, sep = ""), sep = "\n")
chrom <-  read.fasta(paste(g_paht, "/", chr, ".fa", sep = ""), as.string = T)

cat("######################", sep = "\n")
cat(paste("Get all possible methylation sites in the ROI", sep = ""), sep = "\n")
region <- unlist(strsplit(substring(chrom[[1]], start, end), ""))
all_cpgs <- length(grep("c|g", region))  

p_files <- cbind(gsub("_trimmed_bismark_bt2.bismark.cov", "", cov_files), 1:length(cov_files))
# Depends on the nameing scheme you may order the samples
#p_files <- p_files[order(p_files[,3]),]

#Process the input files
if(args$paired == "T"){
  for(i in seq(1,nrow(p_files),2)){  
    cat("######################", sep = "\n")
    cat(paste("Processing sample: ", cov_files[as.numeric(p_files[i,2])], sep = ""), sep = "\n")
    
    #Read in current file
    bism <- read.csv(cov_files[as.numeric(p_files[i,2])], header = F, stringsAsFactors = F, sep = "\t")
    bism <- bism[grep(chr, bism$V1),]

    cat("######################", sep = "\n")
    cat(paste("Processing sample: ", cov_files[as.numeric(p_files[i+1,2])], sep = ""), sep = "\n")
    
    bism2 <- read.csv(cov_files[as.numeric(p_files[i+1,2])], header = F, stringsAsFactors = F, sep = "\t")
    bism2 <- bism2[grep(chr, bism2$V1),]
    
    #Get data for the region of interes
    cpg <- matrix(data = c(c(start:end), 
                           rep(0, length(c(start:end)))), byrow = F, ncol = 2)
    
    cpg[,2] <- bism[match(cpg[,1], bism[,2]),4]
    cpg[is.na(cpg[,2]),2] <- 0
    cpg_table <-  cbind(cpg_table, cpg[,2])
    cpg <- as.data.frame(cpg)
    colnames(cpg) <- c(paste("Position"), "Methylation_level")
    
    cpg2 <- matrix(data = c(c(start:end), 
                           rep(0, length(c(start:end)))), byrow = F, ncol = 2)
    
    cpg2[,2] <- bism2[match(cpg2[,1], bism2[,2]),4]
    cpg2[is.na(cpg2[,2]),2] <- 0
    cpg_table <-  cbind(cpg_table, cpg2[,2])
    cpg2 <- as.data.frame(cpg2)
    colnames(cpg2) <- c(paste("Position"), "Methylation_level")
    
    #Plot results as pairs
    p <-ggplot(data=cpg, aes(x=cpg[,1], y=cpg[,2])) +
      geom_bar(stat="identity", fill="black", colour="black", width=1) +
      ylim(0,100) + 
      ggtitle(gsub("_trimmed_bismark_bt2.bismark.cov", "", cov_files[as.numeric(p_files[i,2])])) +
      xlab(paste("Position ", "chr", chr, ":", start, "-", end, sep = "")) + 
      ylab("Methylation level (%)")
    
    p2 <-ggplot(data=cpg2, aes(x=cpg2[,1], y=cpg2[,2])) +
      geom_bar(stat="identity", fill="black", colour="black", width=1) +
      ylim(0,100) + 
      ggtitle(gsub("_trimmed_bismark_bt2.bismark.cov", "", cov_files[as.numeric(p_files[i+1,2])])) +
      xlab(paste("Position ", "chr", chr, ":", start, "-", end, sep = "")) + 
      ylab("Methylation level (%)")
    
    #Create folder for results
    if (file.exists(gene)){
    } else {
      dir.create(gene)
    }
    
    #Save plots
    ggsave(
      filename = paste(gene, "/", gsub("_trimmed_bismark_bt2.bismark.cov", "", cov_files[as.numeric(p_files[i,2])]),
                       "_and_", gsub("_trimmed_bismark_bt2.bismark.cov", "", cov_files[as.numeric(p_files[i+1,2])]), "_", gene, ".png", sep = ""),
      plot = grid.arrange(p, p2, ncol=1),
      path = NULL,
      scale = 1.5,
      units = "cm",
      width = 12,
      height = 8,
      dpi = 320,
      limitsize = TRUE,
    )
  }
} else {
  for(i in 1:length(cov_files)){  
    cat("######################", sep = "\n")
    cat(paste("Processing sample: ", cov_files[i], sep = ""), sep = "\n")
    
    #Read in current file
    bism <- read.csv(cov_files[i], header = F, stringsAsFactors = F, sep = "\t")
    bism <- bism[grep(chr, bism$V1),]
    
    #Get data for the region of interes
    cpg <- matrix(data = c(c(start:end), 
                         rep(0, length(c(start:end)))), byrow = F, ncol = 2)
  
    cpg[,2] <- bism[match(cpg[,1], bism[,2]),4]
    cpg[is.na(cpg[,2]),2] <- 0
    cpg_table <-  cbind(cpg_table, cpg[,2])
    cpg <- as.data.frame(cpg)
    colnames(cpg) <- c(paste("Position"), "Methylation_level")
    
    #Plot results
    p <-ggplot(data=cpg, aes(x=cpg[,1], y=cpg[,2])) +
      geom_bar(stat="identity", fill="black", colour="black", width=1) +
      ylim(0,100) + 
      ggtitle(gsub("_trimmed_bismark_bt2.bismark.cov", "", cov_files[i])) +
      xlab(paste("Position ", "chr", chr, ":", start, "-", end, sep = "")) + 
      ylab("Methylation level (%)")
    
    #Create folder for results
    if (file.exists(gene)){
    } else {
      dir.create(gene)
    }
    
    #Save plots
    ggsave(
      filename = paste(gene, "/",gsub("_trimmed_bismark_bt2.bismark.cov", "", cov_files[i]), "_", gene, ".png", sep = ""),
      plot = p,
      path = NULL,
      scale = 1.5,
      units = "cm",
      width = 12,
      height = 4,
      dpi = 320,
      limitsize = TRUE,
      )
  }
}
#Write output table
cat("######################", sep = "\n")
cat("Calculating methylatin levels for all samples", sep = "\n")
cpg_perc <- (length(grep(TRUE, cpg_table[,2] > 0))/all_cpgs)*100
  
if(ncol(cpg_table) > 2){
  for(i in 3:ncol(cpg_table)){
    cpg_perc <- c(cpg_perc, (length(grep(TRUE, cpg_table[,i] > 0))/all_cpgs)*100)
    }
}

cat("######################", sep = "\n")
cat("Writing output XLXS files contains methtylation levels", sep = "\n")

if(args$paired == "T"){
  cpg_levels <- cbind(c(gsub("_trimmed_bismark_bt2.bismark.cov", "", cov_files[as.numeric(p_files[,2])])), cpg_perc)
  colnames(cpg_levels) <- c("Sample", "Methtylation level (%)")
  write.xlsx(x = cpg_levels, file = paste(gene, "/", gene, "_methylation_levels.xlsx", sep = ""), col.names = F, row.names = F)
  
  cpg_levels <- as.data.frame(cpg_levels)
  cpg_levels$Sample <- factor(cpg_levels$Sample, levels = cpg_levels$Sample)
  
  p <-ggplot(data=cpg_levels, aes(x=cpg_levels[,1], y=cpg_levels[,2])) +
    geom_bar(stat="identity", fill="black", colour="black", width=0.5) +
    ggtitle(paste("Methylation levels in ", "chr", chr, ":", start, "-", end, sep = "")) +
    xlab("Samples") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ylab("Methylation level (%)")
  
  ggsave(
    filename = paste(gene, "/", gene, "_methylatoion_levels.png", sep = ""),
    plot = p,
    path = NULL,
    scale = 1.5,
    units = "cm",
    width = 12,
    height = 4,
    dpi = 320,
    limitsize = TRUE,
  )
} else {
  cpg_levels <- cbind(c(gsub("_trimmed_bismark_bt2.bismark.cov", "", cov_files)), cpg_perc)
  colnames(cpg_levels) <- c("Sample", "Methtylation level (%)")
  write.xlsx(x = cpg_levels, file = paste(gene, "/", gene, "_methylation_levels.xlsx", sep = ""), col.names = F, row.names = F)
  
  cpg_levels <- as.data.frame(cpg_levels)
  p <-ggplot(data=cpg_levels, aes(x=cpg_levels[,1], y=cpg_levels[,2])) +
    geom_bar(stat="identity", fill="black", colour="black", width=0.5) +
    ggtitle(paste("Methylation levels in ", "chr", chr, ":", start, "-", end, sep = "")) +
    xlab("Samples") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ylab("Methylation level (%)")
  
  ggsave(
    filename = paste(gene, "/", gene, "_methylatoion_levels.png", sep = ""),
    plot = p,
    path = NULL,
    scale = 1.5,
    units = "cm",
    width = 12,
    height = 4,
    dpi = 320,
    limitsize = TRUE,
  )
}

cat("######################", sep = "\n")
cat("Writing output XLXS files contains methtylation patterns", sep = "\n")
cat("Only those position will be marked where is a methylated site", sep = "\n")

if(args$paired == "T"){
  colnames(cpg_table) <- c(paste("Position on chr", chr, sep = ""), gsub("_trimmed_bismark_bt2.bismark.cov", "", cov_files[as.numeric(p_files[,2])]))
  empty_rows <- sum(cpg_table[1,2:ncol(cpg_table)])
  for(i in 2:nrow(cpg_table)){empty_rows <- c(empty_rows, sum(cpg_table[i,2:ncol(cpg_table)]))}
  cpg_table <- cpg_table[empty_rows > 0,] 
  write.xlsx(x = cpg_table, file = paste(gene, "/", gene, "_methylation_patterns.xlsx", sep = ""), col.names = T, row.names = F)
} else  {
  colnames(cpg_table) <- c(paste("Position on chr", chr, sep = ""), gsub("_trimmed_bismark_bt2.bismark.cov", "", cov_files))
  empty_rows <- sum(cpg_table[1,2:ncol(cpg_table)])
  for(i in 2:nrow(cpg_table)){empty_rows <- c(empty_rows, sum(cpg_table[i,2:ncol(cpg_table)]))}
  cpg_table <- cpg_table[empty_rows > 0,] 
  write.xlsx(x = cpg_table, file = paste(gene, "/", gene, "_methylation_patterns.xlsx", sep = ""), col.names = T, row.names = F)
  
}

cat("######################", sep = "\n")
cat("DONE!", sep = "\n")

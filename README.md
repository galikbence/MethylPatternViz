# MethylPatternViz
Simply script for visualization and calculation of methylation levels in a genomic region of interest (ROI) unsing Bismark methylation coverage files (named as *_trimmed_bismark_bt2.bismark.cov). 

The methylation levels are presented as percentages and it is based on the number methylated sites in the ROI devided by all possibly CpG sites within the ROI.The number of all possible CpG sites are extracted from the genome sequence. Therefore the user should donwload the right genome/chr sequence in FASTA format. The script can process and visualize samples as pairs. In this case the user should check the samples are in right (alphabetical) order.

User can spicify the following parameters:

# Name of the cro
  -c, --chr          chromosome
  -s, --start        start site
  -e, --end          stop site
  -g, --gene         plot/legen name or gene name
  -a, --all_samples  process all samples in the folder T(rue) or
                     F(alse)
  --sample           process one or more samples, separeated by comma
  -p, --paired       plot results as pairs, T(rue) or F(alse), only
                     works with even number of samples


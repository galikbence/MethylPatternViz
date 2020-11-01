# MethylPatternViz
Simply script for visualization and calculation of methylation levels in a genomic region of interest (ROI) unsing Bismark methylation coverage files (named as *_trimmed_bismark_bt2.bismark.cov). 

The methylation levels are presented as percentages and it is based on the number methylated sites in the ROI devided by the number of all possible CpG sites within the ROI.The number of all possible CpG sites are extracted from the genome sequence. Therefore, the user should donwload the right genome/chr sequence in FASTA format. The script can process and visualize samples as pairs too. In this case the user should check that the samples are in right (alphabetical) order e.g. pair 1A sample follows pair 1B sample.

# User can spicify the following parameters:
-c, --chr         : name of the chromosome, wihtout the "chr" part.

-s, --start       : start site of the ROI.

-e, --end         : end of ROI.

-g, --gene        : name of the ROI/GENE/promoter. It will be used as the name of the generated plots and output folder too.

-a, --all_samples : process all samples in the workig directory. In this case the user should not specify individual samples.

--sample          : process one or more samples, separeated by comma.

-p, --paired      : plot results as pairs, T(rue) or F(alse), only works with even number of samples.

--path            : aboslute path of the folder contains the chr sequences of the input genome.

# Ouput files - methylation patterns

The script plots the methylation patterns (the percentage of methylation at certain position) for all samples in PNG format. If the samples are indiviudal ones the script will create one file per sample. If samples are paired the script will create one figure per sample pair. The file's names will contian the name of the input samples. The script will create an xlsx sheet named as "GENE"_methylation_patterns.xlsx which contains the input data for the methylation pattern figures.

# Ouput files - methylation levels

The script plots the calculated methylations levels of the ROI for all samples on one figure (named as a "GENE"_methylation_levels.png) in PNG format. Also, it will create an xlsx sheet named as "GENE"_methylation_levels.xlsx which contains the input data for the methylation level figure.

# Usage example for paired samples

```
Rscript MethylPatterViz.R --chr 12 --start 120703649 --end 120806190 --gene TEST -all_samples T --paired T --path /absolute/paht/of/genome/folder
```

It will create a folder named TEST contains all output files.

## Copy-Number calling of targeted panel DNA-seq data using panelcn.mops

####################################
### Jeff Bruce - jeffpbruce@gmail.com
### Script to take in a set of test bams and reference ('normal') bams along
### with a target bed file and output figures and copy-number calls
####################################

####################################
### Required libraries
####################################
library(optparse)
library(panelcn.mops)
library(CopyNumberPlots)

####################################
### Custom Functions
####################################



####################################
### Get inputs from command line argumets
####################################
#
option_list <- list(
  make_option(c("-t", "--test_bams"), default="NA",
              help="Default = 'No default, this parameter is mandatory.'"),
  make_option(c("-r", "--ref_bams"), default="NA",
              help="Default = 'No default, this parameter is mandatory.'"),
  make_option(c("-b", "--bed_file"), default="NA",
              help="Default = 'No default, this parameter is mandatory.'"),
  make_option(c("-o", "--output_directory"), default="NA",
              help="Default = 'No default, this parameter is mandatory.'"),
  make_option(c("-c", "--amp_cutoff"), type = "double", default=2.1,
              help="Default cutoff determined in one study and should not be 
              taken as universal default [default %default]")
)


####################################
### Set variables with test data for development

ref_bams_files        <- "/cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/PanelCNmops/ReferenceSet/sscs_sc/"
test_bams_files       <- "/cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/PanelCNmops/TumourSet/sscs_sc/Screening_indel_realined_TGL59_0028_Pl_n_PE_578_TS_210507_NB551056_0187_AHVVFYBGXH_1_CGGTTGTT-GTGGTATG_R1.fastq.gz.sorted.sscs.sc.sorted.bam"
bed_file              <- "/cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/bed_files/CHARM-MMR_plus_EVOLVE_hg38_liftover.bed"
output_directory      <- "/cluster/projects/ovgroup/projects/OV_Superset/EVOLVE/ctDNA/PanelCNmops/Output_sscs_sc"

####################################


opt <- parse_args(OptionParser(option_list=option_list))

test_bams_files       <- opt$test_bams
ref_bams_files        <- opt$ref_bams
bed_file              <- opt$bed_file
output_directory      <- opt$output_directory


## TODO - check if chromosome have chr prefix

countWindows <- getWindows(bed_file, chr = TRUE)


# load test bam file paths (all bams in folder if input is a directory or single
# file if input is one bam file)

if(dir.exists(test_bams_files)){
  sampleFilenames <- list.files(path = test_bams_files,
                                pattern = ".bam$",
                                full.names = TRUE)
}else{
  sampleFilenames <- test_bams_files
}


# load reference bam file paths
referenceFilenames <- list.files(path = ref_bams_files,
                                 pattern = ".bam$",
                                 full.names = TRUE)

## TODO - This will only work for EVOLVE Panel - need to modify to generalize
## Perhaps just mandate a properly formatted gene column in the bed
##############################################################################
##############################################################################
##############################################################################

gene_list <- c("BRCA1", "BRCA2", "PALB2", "TP53", "APC", "EPCAM", "PMS2", 
               "MLH1", "MSH2", "MSH6", "ABCB1", "CCNE1")

tmp_gene <- countWindows$gene

for(gene in gene_list){
  tmp_gene[grepl(gene,countWindows$gene)] <- gene
}

countWindows$gene <- tmp_gene

tmp_gene <- rep("",nrow(countWindows))

for(gene in gene_list){
  tmp_gene[grepl(gene,countWindows$gene)] <- gene
}

countWindows$gene_Label <- tmp_gene

selectedGenes = "CCNE1"

### Region of interest (CCNE1)
df <- data.frame(chr="chr19", start=29809994, end=29826317,
                 strand=c("."), score=1)

ROI <- makeGRangesFromDataFrame(df)  # strand value "." is replaced with "*"


CCNE1_bed <- read.table("/cluster/projects/pughlab/references/intervals/gene_intervals/hg38/CCNE1_gencode_v31_exons.bed", header = T, sep="\t")

CCNE1_df <- data.frame(chr=CCNE1_bed$Chromsome, start=CCNE1_bed$Start, end=CCNE1_bed$End,
                 strand=CCNE1_bed$Strand, Exon=CCNE1_bed$Exon_Number)

CCNE1_gr <- makeGRangesFromDataFrame(CCNE1_df)  

## Full length of the gene to create intron line
CCNE1_full_df <- data.frame(chr=CCNE1_bed$Chromsome[1], 
                            start=min(CCNE1_bed$Start), 
                            end=max(CCNE1_bed$End),
                       strand="*")

CCNE1_full_gr <- makeGRangesFromDataFrame(CCNE1_full_df)   

##############################################################################
##############################################################################
##############################################################################

## panelcn.mops tries to extract the exon number using a "." as
## a separator in the "name" column

countWindows$name <- make.unique(countWindows$gene)


## Load control files

control <- countBamListInGRanges(countWindows = countWindows,
                                 bam.files = referenceFilenames, 
                                 read.width = 100)


## Loop through tumour samples individually

for (sampleFile in sampleFilenames){
  
  print(paste("Processing",basename(sampleFile)))
  
  test <- countBamListInGRanges(countWindows = countWindows,
                                bam.files = sampleFile, 
                                read.width = 100)
  
  
  XandCB <- test
  
  elementMetadata(XandCB) <- cbind(elementMetadata(XandCB),
                                   elementMetadata(control))
  
  
  resultlist <- runPanelcnMops(XandCB,
                               testiv = 1:ncol(elementMetadata(test)),
                               countWindows = countWindows,
                               selectedGenes = selectedGenes)
  
  sampleNames <- colnames(elementMetadata(test))
  
  resulttable_goi <- createResultTable(resultlist = resultlist, XandCB = XandCB,
                                   countWindows = countWindows,
                                   selectedGenes = selectedGenes,
                                   sampleNames = sampleNames)
  
  resulttable_all <- createResultTable(resultlist = resultlist, XandCB = XandCB,
                                       countWindows = countWindows,
                                       sampleNames = sampleNames)
  
  ## Output results table
  output_directory_sam <- paste0(output_directory,"/", sampleNames[1],"/")
  
  dir.create(output_directory_sam, showWarnings = FALSE, recursive=TRUE)
  
  write.table(resulttable_goi[[1]], sep = "\t", row.names = F, quote = F,
              file = paste0(output_directory_sam,"/", sampleNames[1], "_panelcn.mops.gene_of_interest_results.tsv"))
  
  write.table(resulttable_all[[1]], sep = "\t", row.names = F, quote = F,
              file = paste0(output_directory_sam,"/", sampleNames[1], "_panelcn.mops_results.tsv"))
  
  pdf(paste0(output_directory_sam,"/", sampleNames[1], "_",
             selectedGenes[1], "_panelcn.mops.boxplot.pdf"))
  ## This doesn't work sometimes
  try(
  plotBoxplot(result = resultlist[[1]], sampleName = sampleNames[1],
              countWindows = countWindows,
              selectedGenes = selectedGenes, showGene = 1),
  silent = T)
  dev.off()
  
  
  ## Custom plots
  
  gr_ob <- resultlist[[1]]@gr
  
  gr_ob$Chromosome <- seqnames(gr_ob)
  gr_ob$baf <- resultlist[[1]]@localAssessments[,1]
  ### floor results to show up on plot
  gr_ob$baf[gr_ob$baf<0] <- 0
  
  ## Calculate normalized logR for plotting
  
  norm_res <- resultlist[[1]]@normalizedData
  
  control_samples <- colnames(elementMetadata(control))
  test_sample <- colnames(norm_res)[!(colnames(norm_res) %in% control_samples)]
  
  control_mean <- as.numeric(apply(norm_res[,control_samples],1,mean))
  
  gr_ob$lrr <- log2(norm_res[,test_sample]/control_mean)

  gr_ob$cn <-  as.numeric(gsub("CN","",resultlist[[1]]@integerCopyNumber[,1]))
  gr_ob$loh <- NA
  gr_ob$segment.value  <-  gr_ob$cn
  gr_ob$CopyNumberInteger <- gr_ob$cn
  
  cnplot_dat <- loadSNPData(gr_ob)
  
  all_chrom_file = paste0(output_directory_sam,"/", sampleNames[1], "_",
                          "_panelcn.mops.All_chrom_plot.pdf")
  
  pdf(all_chrom_file,
      height = 8, width = 12)
  plot.params <- getDefaultPlotParams(plot.type=3)
  
  plot.params$leftmargin <- 0.1
  
  kp <- plotKaryotype("hg38", plot.type = 3, labels.plotter = NULL, 
                      main="", cex=1.1, plot.params = plot.params)
  kpAddChromosomeNames(kp, srt=45, cex=1.2,data.panel=2)
  kpAddChromosomeSeparators(kp, lwd=2, col = "#666666")
  
  plotBAF(kp, cnplot_dat, labels = NA,  
          points.col = "dodgerblue2", line.at.0 = FALSE,
          points.cex = 1.1, add.axis = F, data.panel=1)
  
  kpAxis(kp, ymin=0, ymax=1, r0=0, r1=1, col="gray50", cex=1, numticks = 5)
  
  
  kpAddLabels(kp, labels = "Local Assessment Score",srt = 90, 
              cex=1, pos=3, label.margin = 0.07)
  
  plotCopyNumberCallsAsLines(kp, cn.calls = cnplot_dat, 
                             ymin=0, ymax=4, labels="", style="line", col="red",
                             lwd = 2,
                             data.panel=2,add.axis = TRUE,
                             r0=1, r1=0)
  kpAddLabels(kp, labels = "Copy-Number Estimate",srt = 90, 
              cex=1.1, pos=3, label.margin = 0.07, data.panel=2)
  
  dev.off()
  
  ## Output plots for each chromosome
  chrs <- as.character(unique(gr_ob$Chromosome))
  
  for(chr in chrs){
     gene_labels <- subset(countWindows,chromosome == chr)
     ## select the first position for each gene
     gene_labels <- subset(gene_labels,!duplicated(gene_Label))
    
    one_chrom_file <- paste0(output_directory_sam,"/", sampleNames[1],
                                              "_panelcn.mops_",chr, "_plot.pdf")
    
    pdf(one_chrom_file,
        height = 8, width = 12)
    
    plot.params <- getDefaultPlotParams(plot.type=3)
    
    plot.params$leftmargin <- 0.1
    plot.params$rightmargin <- 0.1
    
    kp <- plotKaryotype("hg38", chromosomes = chr, plot.type = 3, 
                        labels.plotter = NULL, 
                        main=chr, cex=1.1, plot.params = plot.params)
    
    plotLRR(kp, cnplot_dat, ymin=-2, ymax=2, labels = NA,  
            points.col = "lightgrey", line.at.0 = FALSE,
            points.cex = 1.1, add.axis = F, data.panel=1, 
            out.of.range.col = "lightgrey")
    
    kpAxis(kp, ymin=-2, ymax=2, r0=0, r1=1, col="gray50", cex=1, numticks = 5)
    
    
    kpAddLabels(kp, labels = "Log2 Normalized Count Ratio",srt = 90, 
                cex=1, pos=3, label.margin = 0.07)
    
    plotCopyNumberCallsAsLines(kp, cn.calls = cnplot_dat, 
                               ymin=0, ymax=4, labels="", style="line", col="red",
                               lwd = 2,
                               data.panel=1,add.axis = F,
                               r0=0, r1=1)
    
    kpAxis(kp, ymin=0, ymax=4, r0=0, r1=1, col="red", cex=1, numticks = 5,
           side=2 )
    
    kpAddLabels(kp, labels = "Copy-Number Estimate",srt = 90, 
                cex=1.1, pos=3, label.margin = 0.07, data.panel=1,
                side="right")
    
    ## label genes
    kpText(kp, chr=chr, y=1, x=gene_labels$start, 
           data.panel = 2, r0=0.2, r1=0, srt=90,
           col="#444444", label=gene_labels$gene_Label, cex=1, 
           pos=2, offset =-0.1,font=2)
    

    dev.off()
  }
  
  ## Zoom in on region of interest 
  ## TODO - Make generalizable ot other genes
  
  zoom_file = paste0(output_directory_sam,"/", sampleNames[1], "_",
                          "_panelcn.mops.CCNE1_plot.pdf")
  
  pdf(zoom_file,
      height = 8, width = 12)
  
  kp <- plotKaryotype("hg38", chromosomes = chr, plot.type = 3, 
                      labels.plotter = NULL, 
                      main=chr, cex=1.1, plot.params = plot.params,
                      zoom = ROI)
  
  plotLRR(kp, cnplot_dat, ymin=-2, ymax=2, labels = NA,  
          points.col = "lightgrey", line.at.0 = FALSE,
          points.cex = 1.1, add.axis = F, data.panel=1, 
          out.of.range.col = "lightgrey")
  
  kpAxis(kp, ymin=-2, ymax=2, r0=0, r1=1, col="gray50", cex=1, numticks = 5)
  
  
  kpAddLabels(kp, labels = "Log2 Normalized Count Ratio",srt = 90, 
              cex=1, pos=3, label.margin = 0.07)
  
  plotCopyNumberCallsAsLines(kp, cn.calls = cnplot_dat, 
                             ymin=0, ymax=4, labels="", style="line", col="red",
                             lwd = 2,
                             data.panel=1,add.axis = F,
                             r0=0, r1=1)
  
  kpAxis(kp, ymin=0, ymax=4, r0=0, r1=1, col="red", cex=1, numticks = 5,
         side=2 )
  
  kpAddLabels(kp, labels = "Copy-Number Estimate",srt = 90, 
              cex=1.1, pos=3, label.margin = 0.07, data.panel=1,
              side="right")
  
  #plot exons as rectangles
  kpPlotRegions(kp, data = CCNE1_gr, r0=0, r1=0.1, data.panel = 2)
  kpPlotNames(kp, data = CCNE1_gr, y0=1, y1=2, 
              labels = CCNE1_df$Exon, position = "bottom", 
              r0=0, r1=0.1, data.panel = 2)
  
  kpSegments(kp, data = CCNE1_full_gr, r0=0, r1=0.1, y1=0.5, y0= 0.5, data.panel = 2,
             col = "black")
  
  kpPlotNames(kp, data = ROI, y0=2, y1=3, 
              labels = "CCNE1", position = "bottom", 
              r0=0, r1=0.1, data.panel = 2, 
              font=2)
  dev.off()

}

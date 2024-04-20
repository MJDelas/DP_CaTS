RNA_3_diffexpression
================

# RNA analysis

Differential expression between gates, conditions and timepoints

``` r
rm(list=ls())

library(DESeq2)
library(RColorBrewer)
library(tidyverse)
```

    ## Warning: package 'stringr' was built under R version 4.2.3

``` r
library(ComplexHeatmap)
library(tximport)
library(ggrepel)
```

### Set dirs

``` r
workingdir="~/Dropbox (The Francis Crick)/DP_cisReg/"
subworkinput="inputs_CaTSRNA/"
outdir="outputs_CaTSRNA_3_diffexpression/"

ifelse(!dir.exists(file.path(workingdir,outdir)), dir.create(file.path(workingdir,outdir)), "Directory exists")
```

    ## [1] "Directory exists"

``` r
suboutdir1="output_between_gates/"
ifelse(!dir.exists(file.path(workingdir,outdir,suboutdir1)), dir.create(file.path(workingdir,outdir,suboutdir1)), "Directory exists")
```

    ## [1] "Directory exists"

``` r
suboutdir2="output_between_timepoints/"
ifelse(!dir.exists(file.path(workingdir,outdir,suboutdir2)), dir.create(file.path(workingdir,outdir,suboutdir2)), "Directory exists")
```

    ## [1] "Directory exists"

``` r
suboutdir3="output_between_conditions/"
ifelse(!dir.exists(file.path(workingdir,outdir,suboutdir3)), dir.create(file.path(workingdir,outdir,suboutdir3)), "Directory exists")
```

    ## [1] "Directory exists"

``` r
suboutdir4="output_between_gates_ignoringothers/"
ifelse(!dir.exists(file.path(workingdir,outdir,suboutdir4)), dir.create(file.path(workingdir,outdir,suboutdir4)), "Directory exists")
```

    ## [1] "Directory exists"

## Load data

For RNA analysis, we are using the output of star_salmon so the import
to deseq is a bit more complicated.

``` r
#salmon counts from pipeline, import
path_files =  list.files(paste0(workingdir,subworkinput)) 
samples = data.frame(run=path_files, stringsAsFactors = FALSE) %>%
  filter(str_detect(run, "^D"))

files <- file.path(paste0(workingdir,subworkinput), samples$run, "quant.sf")
names(files) <- samples$run
all(file.exists(files))
```

    ## [1] TRUE

``` r
#from pipeline
tx2gene = read_tsv(paste0(paste0(workingdir,subworkinput),"/salmon_tx2gene.tsv"))
```

    ## New names:
    ## Rows: 35118 Columns: 3
    ## ── Column specification
    ## ──────────────────────────────────────────────────────── Delimiter: "\t" chr
    ## (3): NM_001011874, Xkr4...2, Xkr4...3
    ## ℹ Use `spec()` to retrieve the full column specification for this data. ℹ
    ## Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## • `Xkr4` -> `Xkr4...2`
    ## • `Xkr4` -> `Xkr4...3`

``` r
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)
```

    ## reading in files with read_tsv
    ## 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 
    ## transcripts missing from tx2gene: 1
    ## summarizing abundance
    ## summarizing counts
    ## summarizing length

## Colors and shapes

``` r
sorted_gate <- c("pMN","DP","p3","neurons")
sorted_conditions <- c("500","UPSAG","dRA2UPSAG","dRA")


shapes4_manual = c(18,15,16,17) # these are block
shapes5_manual = c(25,21,22,23,24) # these are filled
shapes4_fill_manual = c(23,21,22,24)
 

color_gates <- c("#e60000","#cd00cd","#0073e6","#696969")

# for Days
colors_greys <- c("#f6f6f6","#808080","#333333")

# conditions

colors_conditions <- c("#e67300","#4d9a00","#cdcd00","#0073e6")
```

## Differential analysis between domains for a given timepoint and condition

Targeted diff analysis in subsets of samples: - for each timepoint: diff
expression between domains

- make a extra lapply to go through conditions

This DESeq2 analysis is done by subsetting samples. I subset the
`txi.salmon$counts` table and then use `DESeqDataSetFromMatrix` with
`round(sub_counts)` to perform the differential analysis.

``` r
count_matrix <- txi.salmon$counts %>%
  as.data.frame()


#subset 
treatment=c("_500_","_UPSAG_","_dRA_","_dRA2UPSAG_")

timepoint=c("D5_","D6_","D7_")

#matrix
allgates=matrix(c("_pMN_","_p3_","_pMN_","_DP_","_p3_","_DP_"),
                 nrow=2,
                 ncol=3)

comparisons <- allgates

# count_matrix <- count_table %>%
#   column_to_rownames("Geneid")


# x=3
# y=2
# i=3

## many comparisons do not have 2 reps and 2 conditions

PairWiseDEseq <-  lapply(c(1:length(treatment)),function (i) {
  lapply(c(1:length(timepoint)),function (x) {
    lapply(c(1:ncol(allgates)), function (y) {
        regime <- treatment[i]
        timepoints <- timepoint[x]
        celltypes <- allgates[,y]
        sub_counts <- count_matrix %>%
          dplyr::select(contains(regime)  & contains(celltypes) & contains(timepoints))
        
      ## if_else statement because many comparisons do not have enough samples
        
        if(ncol(sub_counts)>=4){
      
          ## Make metadata file for DESeq
          genecolData_sub <- data.frame(Sample = colnames(sub_counts))
          genecolData_sub <- genecolData_sub %>% 
            separate(Sample,into=c("Day","Condition","Gate","Rep"), sep="_", remove=FALSE) %>%
            mutate(DayGate=factor(paste(Day,Gate,sep="_")),
                  DayCondition=paste(Day,Condition),
                  Experiment=paste(Condition,Rep,sep="_"),
                  Gate=factor(Gate, levels=sorted_gate))
          genecolData_sub <- as.data.frame(unclass(genecolData_sub))
          
          
          dds_sub <- DESeqDataSetFromMatrix(countData =  round(sub_counts),
                                  colData = genecolData_sub,
                                  design = ~ Gate)
              
          dds_sub <- DESeq(dds_sub)
          
          vsd_sub <- varianceStabilizingTransformation(dds_sub,blind = FALSE)
          
          # Export normalized tables for plotting elsewhere
          dds_sub_counts <- counts(dds_sub, normalized = TRUE)
          vsd_sub_data <- assay(vsd_sub)
          
          results_sub <- results(dds_sub)
    
          # color significant
          results_sub_plot1 <- results_sub %>%
            as.data.frame() %>%
            rownames_to_column("geneid") %>%
            mutate(color_sig=case_when(padj < 0.1 & padj > 0.01 ~ "under01",
                                       padj < 0.01 & padj > 0.001 ~ "under001",
                                       padj < 0.001 & padj >0 ~ "under0001",
                                       TRUE ~ "over01"))
          
          #threshold <- 8
          plot_ma <- ggplot(results_sub_plot1 %>% as.data.frame(), aes(x=baseMean, y=log2FoldChange, color=color_sig, label=geneid)) +
              geom_point(size=1) +
              # geom_point(data= results_sub_plot1[results_sub_plot1$log2FoldChange > threshold,],
              #           aes(x=baseMean, y=threshold), shape = 2, colour="#d83a00") +
              # geom_point(data= results_sub_plot1[results_sub_plot1$log2FoldChange < -threshold,],
              #           aes(x=baseMean, y=-threshold), shape = 2, colour="#d83a00") +
              #ylim(-threshold,threshold) +
              scale_x_log10() +
              scale_color_manual(values = c("gray30","#d83a00","#ff9b76","#ffd4c4")) +
              geom_text_repel(data = subset(results_sub_plot1, color_sig=="under0001" & baseMean > 100 & log2FoldChange > 0),
                      nudge_y = 10,
                      #nudge_x=-10,
                      #force_pull   = 10,
                      force         = 65,
                      max.overlaps = Inf,
                      box.padding = 0.5,
                      segment.color = "grey50",
                      direction     = "both") +
              geom_text_repel(data = subset(results_sub_plot1, color_sig=="under0001" & baseMean > 100 & log2FoldChange < 0),
                      nudge_y = -20,
                      #nudge_x=-10,
                      #force_pull   = 10,
                      force         = 65,
                      max.overlaps = Inf,
                      box.padding = 0.5,
                      segment.color = "grey50",
                      direction     = "both") +
              ylab(paste0("log2 Fold Change ", resultsNames(dds_sub)[2] ," RNA")) +
              ggtitle(paste0("Results ",regime,timepoints,resultsNames(dds_sub)[2])) +
              theme_bw() 
          
          # print plot
          ggsave(paste0(workingdir,outdir,suboutdir1,"MAplot_",regime,timepoints,resultsNames(dds_sub)[2],".pdf"), plot=plot_ma,
           width=5, height=3, units="in", useDingbats=FALSE)
          
    
          ## Export files
          
          write.table(dds_sub_counts,
          file = paste0(workingdir,outdir,suboutdir1,"CountsNormalized_",regime,"_",timepoints,resultsNames(dds_sub)[2],".txt"),
              quote = FALSE, row.names = TRUE)
          write.csv(vsd_sub_data,
              paste0(workingdir,outdir,suboutdir1,"VSData_",regime,"_",timepoints,resultsNames(dds_sub)[2],".csv"),
              quote = FALSE)
          write.table(results_sub,
              file = paste0(workingdir,outdir,suboutdir1,"Results_DESeq_",regime,"_",timepoints,resultsNames(dds_sub)[2],".txt"),
              quote = FALSE, row.names = TRUE)
    
          results_return <- results_sub %>% as.data.frame() %>% rownames_to_column("Geneid")
          #results_return$Comparison <- paste0("Comp_",timepoints,celltypes[1],celltypes[2])
          results_return$Comparison <- paste0("Comp_",regime,timepoints,resultsNames(dds_sub)[2])
          results_return
        } 
  
    })
  })
})
```

## Differential analysis between timepoints for each domain and condition

For each gate: - differential analysis between timepoints for the same
condition: To do, extra lapply for conditions

``` r
#subset 
treatment=c("_500_","_UPSAG_","_dRA_","_dRA2UPSAG_")

allgates=c("_pMN_","_p3_","_DP_")
#matrix
timepoints=matrix(c("D5_","D6_","D5_","D7_","D6_","D7_"),
                 nrow=2,
                 ncol=3)

comparisons <- timepoints

# i=1
# x=2
# y=3

PairWiseDEseq <- lapply(c(1:length(treatment)),function (i) {
  lapply(c(1:length(allgates)),function (x) {
    lapply(c(1:ncol(timepoints)), function (y) {
        regime <- treatment[i]
        celltypes <- allgates[x]
        days <- timepoints[,y]
        sub_counts <- count_matrix %>%
          dplyr::select(contains(regime)  & contains(celltypes) & contains(days))
        
         ## if_else statement because many comparisons do not have enough samples
        
        if(ncol(sub_counts)>=4){
        
        ## Make metadata file for DESeq
        genecolData_sub <- data.frame(Sample = colnames(sub_counts))
        genecolData_sub <- genecolData_sub %>% 
          separate(Sample,into=c("Day","Condition","Gate","Rep"), sep="_", remove=FALSE) %>%
          mutate(DayGate=factor(paste(Day,Gate,sep="_")),
                DayCondition=paste(Day,Condition),
                Experiment=paste(Condition,Rep,sep="_"),
                Gate=factor(Gate, levels=sorted_gate))
        genecolData_sub <- as.data.frame(unclass(genecolData_sub))
           
        dds_sub <- DESeqDataSetFromMatrix(countData = round(sub_counts),
                                      colData = genecolData_sub,
                                      design = ~ Day)
        
        dds_sub <- DESeq(dds_sub)
        
        vsd_sub <- varianceStabilizingTransformation(dds_sub,blind = FALSE)
        
        # Export normalized tables for plotting elsewhere
        dds_sub_counts <- counts(dds_sub, normalized = TRUE)
        vsd_sub_data <- assay(vsd_sub)
        
        results_sub <- results(dds_sub)
        
        ## plot MA custom plot 
        
        #plotMA(results_sub,ylim=c(-8,8))
        
        # color significant
        results_sub_plot1 <- results_sub %>%
          as.data.frame() %>%
          rownames_to_column("geneid") %>%
          mutate(color_sig=case_when(padj < 0.1 & padj > 0.01 ~ "under01",
                                     padj < 0.01 & padj > 0.001 ~ "under001",
                                     padj < 0.001 & padj >0 ~ "under0001",
                                     TRUE ~ "over01"))
        
        #threshold <- 8
         plot_ma <- ggplot(results_sub_plot1 %>% as.data.frame(), aes(x=baseMean, y=log2FoldChange, color=color_sig, label=geneid)) +
            geom_point(size=1) +
            # geom_point(data= results_sub_plot1[results_sub_plot1$log2FoldChange > threshold,],
            #           aes(x=baseMean, y=threshold), shape = 2, colour="#d83a00") +
            # geom_point(data= results_sub_plot1[results_sub_plot1$log2FoldChange < -threshold,],
            #           aes(x=baseMean, y=-threshold), shape = 2, colour="#d83a00") +
            #ylim(-threshold,threshold) +
            scale_x_log10() +
            scale_color_manual(values = c("gray30","#d83a00","#ff9b76","#ffd4c4")) +
            geom_text_repel(data = subset(results_sub_plot1, color_sig=="under0001" & baseMean > 100 & log2FoldChange > 0),
                    nudge_y = 10,
                    #nudge_x=-10,
                    #force_pull   = 10,
                    force         = 65,
                    max.overlaps = Inf,
                    box.padding = 0.5,
                    segment.color = "grey50",
                    direction     = "both") +
            geom_text_repel(data = subset(results_sub_plot1, color_sig=="under0001" & baseMean > 100 & log2FoldChange < 0),
                    nudge_y = -20,
                    #nudge_x=-10,
                    #force_pull   = 10,
                    force         = 65,
                    max.overlaps = Inf,
                    box.padding = 0.5,
                    segment.color = "grey50",
                    direction     = "both") +
            ylab(paste0("log2 Fold Change ",resultsNames(dds_sub)[2]," RNA")) +
            ggtitle(paste0("Results ",regime,celltypes,resultsNames(dds_sub)[2] )) +
            theme_bw() 
        
         # print plot
          ggsave(paste0(workingdir,outdir,suboutdir2,"MAplot_",regime,celltypes,resultsNames(dds_sub)[2],".pdf"), plot=plot_ma,
           width=5, height=3, units="in", useDingbats=FALSE)
        
        ## Export files
        
        write.table(dds_sub_counts,
        file = paste0(workingdir,outdir,suboutdir2,"CountsNormalized_",regime,celltypes,resultsNames(dds_sub)[2],".txt"),
            quote = FALSE, row.names = TRUE)
        write.csv(vsd_sub_data,
            paste0(workingdir,outdir,suboutdir2,"VSData_",regime,celltypes,resultsNames(dds_sub)[2],".csv"),
            quote = FALSE)
        write.table(results_sub,
            file = paste0(workingdir,outdir,suboutdir2,"Results_DESeq_",regime,celltypes,resultsNames(dds_sub)[2],".txt"),
            quote = FALSE, row.names = TRUE)

        results_return <- results_sub %>% as.data.frame() %>% rownames_to_column("Geneid")
        results_return$Comparison <- paste0("Comp_",regime,celltypes,resultsNames(dds_sub)[2])
        results_return
        }
    })
  }) 
}) 
```

## Differential analysis between conditions for each domain and timepoint

For a given timepoint - differential analysis between conditions each
cell type

``` r
#subset
timepoint=c("D5_","D6_","D7_")
allgates=c("_pMN_","_p3_","_DP_","_neurons_")
#matrix
comparisons=matrix(c("_500_","_UPSAG_","_UPSAG_","_dRA2UPSAG_","_dRA2UPSAG_","_dRA_",
                    "_500_","_dRA2UPSAG_","_500_","_dRA_","_UPSAG_","_dRA_"),
                 nrow=2,
                 ncol=6)


# i=3
# x=2
# y=5

PairWiseDEseq <-  lapply(c(1:length(timepoint)),function (i) {
  lapply(c(1:length(allgates)),function (x) {
    lapply(c(1:ncol(comparisons)), function (y) {
        timepoints <- timepoint[i]
        celltypes <- allgates[x]
        conditions <- comparisons[,y]
        sub_counts <- count_matrix %>%
          dplyr::select(contains(timepoints)  & contains(celltypes) & contains(conditions))
        
         ## if_else statement because many comparisons do not have enough samples
        if(ncol(sub_counts)>=4){
          
        ## Make metadata file for DESeq
        genecolData_sub <- data.frame(Sample = colnames(sub_counts))
        genecolData_sub <- genecolData_sub %>% 
          separate(Sample,into=c("Day","Condition","Gate","Rep"), sep="_", remove=FALSE) %>%
          mutate(DayGate=factor(paste(Day,Gate,sep="_")),
                DayCondition=paste(Day,Condition),
                Experiment=paste(Condition,Rep,sep="_"),
                Gate=factor(Gate, levels=sorted_gate))
        genecolData_sub <- as.data.frame(unclass(genecolData_sub))
           
        dds_sub <- DESeqDataSetFromMatrix(countData = round(sub_counts),
                                      colData = genecolData_sub,
                                      design = ~ Condition)
        
        dds_sub <- DESeq(dds_sub)
        
        vsd_sub <- varianceStabilizingTransformation(dds_sub,blind = FALSE)
        
        # Export normalized tables for plotting elsewhere
        dds_sub_counts <- counts(dds_sub, normalized = TRUE)
        vsd_sub_data <- assay(vsd_sub)
        
        results_sub <- results(dds_sub)
        
        ## plot MA custom plot 
        
        #plotMA(results_sub,ylim=c(-8,8))
        
        # color significant
        results_sub_plot1 <- results_sub %>%
          as.data.frame() %>%
          rownames_to_column("geneid") %>%
          mutate(color_sig=case_when(padj < 0.1 & padj > 0.01 ~ "under01",
                                     padj < 0.01 & padj > 0.001 ~ "under001",
                                     padj < 0.001 & padj >0 ~ "under0001",
                                     TRUE ~ "over01"))
        
        #threshold <- 8
        plot_ma <-ggplot(results_sub_plot1 %>% as.data.frame(), aes(x=baseMean, y=log2FoldChange, color=color_sig, label=geneid)) +
            geom_point(size=1) +
            # geom_point(data= results_sub_plot1[results_sub_plot1$log2FoldChange > threshold,],
            #           aes(x=baseMean, y=threshold), shape = 2, colour="#d83a00") +
            # geom_point(data= results_sub_plot1[results_sub_plot1$log2FoldChange < -threshold,],
            #           aes(x=baseMean, y=-threshold), shape = 2, colour="#d83a00") +
            #ylim(-threshold,threshold) +
            scale_x_log10() +
            scale_color_manual(values = c("gray30","#d83a00","#ff9b76","#ffd4c4")) +
            geom_text_repel(data = subset(results_sub_plot1, color_sig=="under0001" & baseMean > 100 & log2FoldChange > 0),
                    nudge_y = 20,
                    #nudge_x=-10,
                    #force_pull   = 10,
                    force         = 65,
                    max.overlaps = Inf,
                    box.padding = 0.5,
                    segment.color = "grey50",
                    direction     = "both") +
            geom_text_repel(data = subset(results_sub_plot1, color_sig=="under0001" & baseMean > 100 & log2FoldChange < 0),
                    nudge_y = -20,
                    #nudge_x=-10,
                    #force_pull   = 10,
                    force         = 65,
                    max.overlaps = Inf,
                    box.padding = 0.5,
                    segment.color = "grey50",
                    direction     = "both") +
            ylab(paste0("log2 Fold Change ",resultsNames(dds_sub)[2]," RNA")) +
            ggtitle(paste0("Results ",timepoints,celltypes,resultsNames(dds_sub)[2])) +
            theme_bw() 
        
        # print plot
          ggsave(paste0(workingdir,outdir,suboutdir3,"MAplot_",timepoints,celltypes,resultsNames(dds_sub)[2],".pdf"), plot=plot_ma,
           width=5, height=3, units="in", useDingbats=FALSE)
        
        
        ## Export files
        
        write.table(dds_sub_counts,
        file = paste0(workingdir,outdir,suboutdir3,"CountsNormalized_",timepoints,celltypes,resultsNames(dds_sub)[2],".txt"),
            quote = FALSE, row.names = TRUE)
        write.csv(vsd_sub_data,
            paste0(workingdir,outdir,suboutdir3,"VSData_",timepoints,celltypes,resultsNames(dds_sub)[2],".csv"),
            quote = FALSE)
        write.table(results_sub,
            file = paste0(workingdir,outdir,suboutdir3,"Results_DESeq_",timepoints,celltypes,resultsNames(dds_sub)[2],".txt"),
            quote = FALSE, row.names = TRUE)

        results_return <- results_sub %>% as.data.frame() %>% rownames_to_column("Geneid")
        results_return$Comparison <- paste0("Comp_",timepoints,celltypes,resultsNames(dds_sub)[2])
        results_return
        }
    })
  })
})
```

## Between cell types: ignoring signaling regime and timepoint

- to get the main genes, but I will also intersect all the different
  comparisons in the next script
- exclude D7 dRA p3 because they are Floor plate

``` r
count_matrix <- txi.salmon$counts %>%
  as.data.frame()


#subset 
treatment=c("_500_","_UPSAG_","_dRA2UPSAG_")
#matrix
allgates=matrix(c("_pMN_","_p3_","_pMN_","_DP_","_p3_","_DP_"),
                 nrow=2,
                 ncol=3)

comparisons <- allgates

y=3


## many comparisons do not have 2 reps and 2 conditions

PairWiseDEseq <-
    lapply(c(1:ncol(allgates)), function (y) {
        celltypes <- allgates[,y]
        sub_counts <- count_matrix %>%
          dplyr::select(contains(celltypes) & contains(treatment))
  
      
          ## Make metadata file for DESeq
          genecolData_sub <- data.frame(Sample = colnames(sub_counts))
          genecolData_sub <- genecolData_sub %>% 
            separate(Sample,into=c("Day","Condition","Gate","Rep"), sep="_", remove=FALSE) %>%
            mutate(DayGate=factor(paste(Day,Gate,sep="_")),
                  DayCondition=paste(Day,Condition),
                  Experiment=paste(Condition,Rep,sep="_"),
                  Gate=factor(Gate, levels=sorted_gate))
          genecolData_sub <- as.data.frame(unclass(genecolData_sub))
          
          
          dds_sub <- DESeqDataSetFromMatrix(countData =  round(sub_counts),
                                  colData = genecolData_sub,
                                  design = ~ Gate)
              
          dds_sub <- DESeq(dds_sub)
          
          vsd_sub <- varianceStabilizingTransformation(dds_sub,blind = FALSE)
          
          # Export normalized tables for plotting elsewhere
          dds_sub_counts <- counts(dds_sub, normalized = TRUE)
          vsd_sub_data <- assay(vsd_sub)
          
          results_sub <- results(dds_sub)
    
          # color significant
          results_sub_plot1 <- results_sub %>%
            as.data.frame() %>%
            rownames_to_column("geneid") %>%
            mutate(color_sig=case_when(padj < 0.1 & padj > 0.01 ~ "under01",
                                       padj < 0.01 & padj > 0.001 ~ "under001",
                                       padj < 0.001 & padj >0 ~ "under0001",
                                       TRUE ~ "over01"))
          
          #threshold <- 8
          plot_ma <- ggplot(results_sub_plot1 %>% as.data.frame(), aes(x=baseMean, y=log2FoldChange, color=color_sig, label=geneid)) +
              geom_point(size=1) +
              # geom_point(data= results_sub_plot1[results_sub_plot1$log2FoldChange > threshold,],
              #           aes(x=baseMean, y=threshold), shape = 2, colour="#d83a00") +
              # geom_point(data= results_sub_plot1[results_sub_plot1$log2FoldChange < -threshold,],
              #           aes(x=baseMean, y=-threshold), shape = 2, colour="#d83a00") +
              #ylim(-threshold,threshold) +
              scale_x_log10() +
              scale_color_manual(values = c("gray30","#d83a00","#ff9b76","#ffd4c4")) +
              geom_text_repel(data = subset(results_sub_plot1, color_sig=="under0001" & baseMean > 100 & log2FoldChange > 2),
                      nudge_y = 30,
                      #nudge_x=-10,
                      #force_pull   = 10,
                      force         = 65,
                      max.overlaps = Inf,
                      box.padding = 0.5,
                      segment.color = "grey50",
                      direction     = "both") +
              geom_text_repel(data = subset(results_sub_plot1, color_sig=="under0001" & baseMean > 100 & log2FoldChange < -2),
                      nudge_y = -30,
                      #nudge_x=-10,
                      #force_pull   = 10,
                      force         = 65,
                      max.overlaps = Inf,
                      box.padding = 0.5,
                      segment.color = "grey50",
                      direction     = "both") +
              ylab(paste0("log2 Fold Change ", resultsNames(dds_sub)[2] ," RNA")) +
              ggtitle(paste0("Results ",resultsNames(dds_sub)[2])) +
              theme_bw() 
          
          # print plot
          ggsave(paste0(workingdir,outdir,suboutdir4,"MAplot_",resultsNames(dds_sub)[2],".pdf"), plot=plot_ma,
           width=5, height=3, units="in", useDingbats=FALSE)
          
    
          ## Export files
          
          write.table(dds_sub_counts,
          file = paste0(workingdir,outdir,suboutdir4,"CountsNormalized_",resultsNames(dds_sub)[2],".txt"),
              quote = FALSE, row.names = TRUE)
          write.csv(vsd_sub_data,
              paste0(workingdir,outdir,suboutdir4,"VSData_",resultsNames(dds_sub)[2],".csv"),
              quote = FALSE)
          write.table(results_sub,
              file = paste0(workingdir,outdir,suboutdir4,"Results_DESeq_",resultsNames(dds_sub)[2],".txt"),
              quote = FALSE, row.names = TRUE)
    
          results_return <- results_sub %>% as.data.frame() %>% rownames_to_column("Geneid")
          #results_return$Comparison <- paste0("Comp_",timepoints,celltypes[1],celltypes[2])
          results_return$Comparison <- paste0("Comp_",resultsNames(dds_sub)[2])
          results_return
})
```

### Targetted p3 analysis

Compare upsag d7 and 500 d6. Is anything different from 500 d7 vs d6?

``` r
samples_wanted <- c("D7_UPSAG_p3","D5_500_p3")
        # timepoints <- timepoint[i]
        # celltypes <- allgates[x]
        # conditions <- comparisons[,y]
        sub_counts <- count_matrix %>%
          dplyr::select(contains(samples_wanted))
        

        ## Make metadata file for DESeq
        genecolData_sub <- data.frame(Sample = colnames(sub_counts))
        genecolData_sub <- genecolData_sub %>% 
          separate(Sample,into=c("Day","Condition","Gate","Rep"), sep="_", remove=FALSE) %>%
          mutate(DayGate=factor(paste(Day,Gate,sep="_")),
                DayCondition=paste(Day,Condition),
                Experiment=paste(Condition,Rep,sep="_"),
                Gate=factor(Gate, levels=sorted_gate))
        genecolData_sub <- as.data.frame(unclass(genecolData_sub))
           
        dds_sub <- DESeqDataSetFromMatrix(countData = round(sub_counts),
                                      colData = genecolData_sub,
                                      design = ~ Day)
        
        dds_sub <- DESeq(dds_sub)
        
        vsd_sub <- varianceStabilizingTransformation(dds_sub,blind = FALSE)
        
        # Export normalized tables for plotting elsewhere
        dds_sub_counts <- counts(dds_sub, normalized = TRUE)
        vsd_sub_data <- assay(vsd_sub)
        
        results_sub <- results(dds_sub)
        
        ## plot MA custom plot 
        
        #plotMA(results_sub,ylim=c(-8,8))
        
        # color significant
        results_sub_plot1 <- results_sub %>%
          as.data.frame() %>%
          rownames_to_column("geneid") %>%
          mutate(color_sig=case_when(padj < 0.1 & padj > 0.01 ~ "under01",
                                     padj < 0.01 & padj > 0.001 ~ "under001",
                                     padj < 0.001 & padj >0 ~ "under0001",
                                     TRUE ~ "over01"))
        
        #threshold <- 8
        plot_ma <-ggplot(results_sub_plot1 %>% as.data.frame(), aes(x=baseMean, y=log2FoldChange, color=color_sig, label=geneid)) +
            geom_point(size=1) +
            # geom_point(data= results_sub_plot1[results_sub_plot1$log2FoldChange > threshold,],
            #           aes(x=baseMean, y=threshold), shape = 2, colour="#d83a00") +
            # geom_point(data= results_sub_plot1[results_sub_plot1$log2FoldChange < -threshold,],
            #           aes(x=baseMean, y=-threshold), shape = 2, colour="#d83a00") +
            #ylim(-threshold,threshold) +
            scale_x_log10() +
            scale_color_manual(values = c("gray30","#d83a00","#ff9b76","#ffd4c4")) +
            geom_text_repel(data = subset(results_sub_plot1, color_sig=="under0001" & baseMean > 100 & log2FoldChange > 2),
                    nudge_y = 20,
                    #nudge_x=-10,
                    #force_pull   = 10,
                    force         = 65,
                    max.overlaps = Inf,
                    box.padding = 0.5,
                    segment.color = "grey50",
                    direction     = "both") +
            geom_text_repel(data = subset(results_sub_plot1, color_sig=="under0001" & baseMean > 100 & log2FoldChange < -2),
                    nudge_y = -20,
                    #nudge_x=-10,
                    #force_pull   = 10,
                    force         = 65,
                    max.overlaps = Inf,
                    box.padding = 0.5,
                    segment.color = "grey50",
                    direction     = "both") +
            ylab(paste0("log2 Fold Change ",resultsNames(dds_sub)[2]," RNA")) +
            ggtitle(paste0("Comparison ",samples_wanted[1]," and ",samples_wanted[2],", order ",resultsNames(dds_sub)[2])) +
            theme_bw() 
        
        # print plot
          ggsave(paste0(workingdir,outdir,"MAplot_",samples_wanted[1]," and ",samples_wanted[2],", order ", resultsNames(dds_sub)[2],".pdf"), plot=plot_ma,
           width=5, height=3, units="in", useDingbats=FALSE)
        
        
        ## Export files
        
        write.table(dds_sub_counts,
        file = paste0(workingdir,outdir,"CountsNormalized_",samples_wanted[1]," and ",samples_wanted[2],", order ", resultsNames(dds_sub)[2],".txt"),
            quote = FALSE, row.names = TRUE)
        write.csv(vsd_sub_data,
            paste0(workingdir,outdir,"VSData_",samples_wanted[1]," and ",samples_wanted[2],", order ",resultsNames(dds_sub)[2],".csv"),
            quote = FALSE)
        write.table(results_sub,
            file = paste0(workingdir,outdir,"Results_DESeq_",samples_wanted[1]," and ",samples_wanted[2],", order ", resultsNames(dds_sub)[2],".txt"),
            quote = FALSE, row.names = TRUE)

        # results_return <- results_sub %>% as.data.frame() %>% rownames_to_column("Geneid")
        # results_return$Comparison <- paste0("Comp_",resultsNames(dds_sub)[2])
        # results_return
```

``` r
sessionInfo()
```

    ## R version 4.2.2 (2022-10-31)
    ## Platform: aarch64-apple-darwin20 (64-bit)
    ## Running under: macOS 14.4.1
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] grid      stats4    stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] ggrepel_0.9.4               tximport_1.26.1            
    ##  [3] ComplexHeatmap_2.15.4       lubridate_1.9.3            
    ##  [5] forcats_1.0.0               stringr_1.5.1              
    ##  [7] dplyr_1.1.3                 purrr_1.0.2                
    ##  [9] readr_2.1.4                 tidyr_1.3.0                
    ## [11] tibble_3.2.1                ggplot2_3.4.4              
    ## [13] tidyverse_2.0.0             RColorBrewer_1.1-3         
    ## [15] DESeq2_1.38.3               SummarizedExperiment_1.28.0
    ## [17] Biobase_2.58.0              MatrixGenerics_1.10.0      
    ## [19] matrixStats_1.1.0           GenomicRanges_1.50.2       
    ## [21] GenomeInfoDb_1.34.9         IRanges_2.32.0             
    ## [23] S4Vectors_0.36.2            BiocGenerics_0.44.0        
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] bitops_1.0-7           bit64_4.0.5            doParallel_1.0.17     
    ##  [4] httr_1.4.7             tools_4.2.2            utf8_1.2.4            
    ##  [7] R6_2.5.1               DBI_1.1.3              colorspace_2.1-0      
    ## [10] GetoptLong_1.0.5       withr_2.5.2            tidyselect_1.2.0      
    ## [13] bit_4.0.5              compiler_4.2.2         textshaping_0.3.7     
    ## [16] cli_3.6.1              DelayedArray_0.24.0    labeling_0.4.3        
    ## [19] scales_1.2.1           systemfonts_1.0.5      digest_0.6.33         
    ## [22] rmarkdown_2.25         XVector_0.38.0         pkgconfig_2.0.3       
    ## [25] htmltools_0.5.7        fastmap_1.1.1          rlang_1.1.2           
    ## [28] GlobalOptions_0.1.2    rstudioapi_0.15.0      RSQLite_2.3.3         
    ## [31] farver_2.1.1           shape_1.4.6            generics_0.1.3        
    ## [34] jsonlite_1.8.7         vroom_1.6.4            BiocParallel_1.32.6   
    ## [37] RCurl_1.98-1.13        magrittr_2.0.3         GenomeInfoDbData_1.2.9
    ## [40] Matrix_1.6-3           Rcpp_1.0.11            munsell_0.5.0         
    ## [43] fansi_1.0.5            lifecycle_1.0.4        stringi_1.8.1         
    ## [46] yaml_2.3.7             zlibbioc_1.44.0        blob_1.2.4            
    ## [49] parallel_4.2.2         crayon_1.5.2           lattice_0.22-5        
    ## [52] Biostrings_2.66.0      annotate_1.76.0        circlize_0.4.15       
    ## [55] hms_1.1.3              KEGGREST_1.38.0        locfit_1.5-9.8        
    ## [58] knitr_1.45             pillar_1.9.0           rjson_0.2.21          
    ## [61] geneplotter_1.76.0     codetools_0.2-19       XML_3.99-0.15         
    ## [64] glue_1.6.2             evaluate_0.23          png_0.1-8             
    ## [67] vctrs_0.6.4            tzdb_0.4.0             foreach_1.5.2         
    ## [70] gtable_0.3.4           clue_0.3-65            cachem_1.0.8          
    ## [73] xfun_0.43              xtable_1.8-4           ragg_1.2.6            
    ## [76] iterators_1.0.14       AnnotationDbi_1.60.2   memoise_2.0.1         
    ## [79] cluster_2.1.4          timechange_0.2.0

R_quick_plot genes
================

# RNA analysis

Gene plotting with output from previous script.

``` r
rm(list=ls())

library(RColorBrewer)
library(tidyverse)
library(ComplexHeatmap)
```

### Set dirs

``` r
workingdir="~/Dropbox (The Francis Crick)/DP_cisReg/"
subworkinput="outputs_CaTSRNA_1/"
```

## Load RNA data

Normalised RNA counts

``` r
dds_counts <- read.table(paste0(workingdir,subworkinput,"featurecounts.normCounts.txt"),stringsAsFactors =FALSE)
```

## Colors and shapes

Before more plotting, let’s get some metadata organised

``` r
sorted_gate <- c("pMN","DP","p3","neurons")
sorted_conditions <- c("500","UPSAG","dRA2UPSAG","dRA")


shapes4_manual = c(18,15,16,17) # these are block
shapes5_manual = c(25,21,22,23,24) # these are filled
shapes4_fill_manual = c(23,21,22,24)
# 
# # for Days
# red_colors <- c("#fadede","#f3aaaa","#e96666","#cf1e1e")


#color_gates <- c("#b30000","#800080","#009640","#696969")
#color_gates <- c("#b30000","#9a009a","#005ab3","#696969")
color_gates <- c("#e60000","#cd00cd","#0073e6","#696969")

# for Days
colors_greys <- c("#f6f6f6","#808080","#333333")

# conditions
#colors_conditions <- c("#b81bb8","#b81b1b","#1bb81b","#1bb8b8")
colors_conditions <- c("#e67300","#4d9a00","#cdcd00","#0073e6")
```

## Genes to look at

``` r
dds_counts_plot <- dds_counts %>% 
  as.data.frame() %>%
  rownames_to_column("geneid") %>%
  gather(sampleid, counts_norm, starts_with("D")) %>%
  separate(sampleid,into=c("Day","Condition","Gate","Rep"), sep="_", remove=FALSE) %>%
  mutate(DayGate=factor(paste(Day,Gate,sep="_")),
         DayCondition=paste(Day,Condition),
         Gate=factor(Gate, levels=sorted_gate),
         Condition=factor(Condition, levels=sorted_conditions))
```

### Basic QC

``` r
geneOI <- c("Sox2",
            "Dbx1","Irx3","Pax6","Nkx6-1","Olig2","Nkx2-2","Foxa2","Shh","Arx",
            "Tubb3","Sim1","Mnx1")


ggplot(dds_counts_plot %>% filter(geneid %in% geneOI) %>% mutate(geneid=factor(geneid, levels=geneOI)), 
       aes(x=Gate,y=counts_norm)) +
  stat_summary(aes(fill=Condition),
    fun = mean, geom="bar", alpha=0.9, width=0.7,position=position_dodge(0.7)) +
  geom_point(aes(fill=Condition), alpha=0.6, position = position_dodge(width = 0.7),color="black") +
  #geom_col(position="dodge",aes(fill=DayGate)) +
  scale_fill_manual(values=colors_conditions) +
  scale_color_manual(values=colors_conditions) +
  scale_shape_manual(values=shapes4_fill_manual) +
  facet_grid(geneid ~ Day, scales = "free_y") +
  theme_bw()
```

![](DPCaTSRNA_2_quick_plotgenes_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
ggplot(dds_counts_plot %>% filter(geneid %in% geneOI) %>% mutate(geneid=factor(geneid, levels=geneOI)), 
       aes(x=Day,y=counts_norm)) +
  stat_summary(aes(fill=Gate),
    fun = mean, geom="bar", alpha=0.9, width=0.7,position=position_dodge(0.7)) +
  geom_point(aes(fill=Gate), alpha=0.6, position = position_dodge(width = 0.7),color="black") +
  #geom_col(position="dodge",aes(fill=DayGate)) +
  scale_fill_manual(values=color_gates) +
  scale_color_manual(values=color_gates) +
  scale_shape_manual(values=shapes4_fill_manual) +
  facet_grid(geneid ~ Condition, scales = "free_y") +
  theme_bw()
```

![](DPCaTSRNA_2_quick_plotgenes_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

### Temporal genes from pilot analysis

p3 d6 to d7

``` r
geneOI <- c("Lin28a","Nfib","Rbp1","H19","Fabp7","Hoxb6","Nkx6-2","Phox2b")


ggplot(dds_counts_plot %>% filter(geneid %in% geneOI) %>% mutate(geneid=factor(geneid, levels=geneOI)), 
       aes(x=Condition,y=counts_norm)) +
  stat_summary(aes(fill=Day),
    fun = mean, geom="bar", alpha=0.9, width=0.7,position=position_dodge(0.7)) +
  geom_point(aes(fill=Day), alpha=0.6, position = position_dodge(width = 0.7),color="black") +
  #geom_col(position="dodge",aes(fill=DayGate)) +
  #scale_fill_manual(values=color_gates) +
  #scale_color_manual(values=color_gates) +
  scale_shape_manual(values=shapes4_fill_manual) +
  facet_grid(geneid ~ Gate, scales = "free_y") +
  theme_bw()
```

![](DPCaTSRNA_2_quick_plotgenes_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

### What are those neurons?

Just D7:

- dRA: Fev and Gata3 point to more anterior (5HT generation?)
- dRA also high in Pou3f1 (Alcam), expressed in spinal accessory column
  (Cervical level)

``` r
geneOI <- c("Tubb3","Mnx1","Isl2","Isl1","Neurog2","Lhx3","Lhx1",
            "Slc18a3","Irx3",
            "Nkx2-2","Sim1","Neurog3","Phox2b","Pou3f1","Alcam",
            "Abca1","Fgf10",
            "Gata3","Fev","Slc17a8","Cyp26b1",
            "Pou4f2","Lpar1")


ggplot(dds_counts_plot %>% filter(geneid %in% geneOI & Day=="D7") %>% mutate(geneid=factor(geneid, levels=geneOI)), 
       aes(x=Gate,y=counts_norm)) +
  stat_summary(aes(fill=Condition),
    fun = mean, geom="bar", alpha=0.9, width=0.7,position=position_dodge(0.7)) +
  geom_point(aes(fill=Condition), alpha=0.6, position = position_dodge(width = 0.7),color="black") +
  #geom_col(position="dodge",aes(fill=DayGate)) +
  scale_fill_manual(values=colors_conditions) +
  scale_color_manual(values=colors_conditions) +
  scale_shape_manual(values=shapes4_fill_manual) +
  facet_wrap(~ geneid, scales = "free_y", ncol=4) +
  theme_bw()
```

![](DPCaTSRNA_2_quick_plotgenes_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

### Hoxes

``` r
#geneOI <- grep("Hoxc", dds_counts_plot$geneid, value = TRUE)

geneOI <- c("Hoxa1","Hoxa3","Hoxa5","Hoxb2","Hoxb5","Hoxb7","Hoxb8","Hoxc4","Hoxc5","Hoxc6")

ggplot(dds_counts_plot %>% filter(geneid %in% geneOI & Day=="D7"), 
       aes(x=Gate,y=counts_norm)) +
  stat_summary(aes(fill=Condition),
    fun = mean, geom="bar", alpha=0.9, width=0.7,position=position_dodge(0.7)) +
  geom_point(aes(fill=Condition), alpha=0.6, position = position_dodge(width = 0.7),color="black") +
  #geom_col(position="dodge",aes(fill=DayGate)) +
  scale_fill_manual(values=colors_conditions) +
  scale_color_manual(values=colors_conditions) +
  scale_shape_manual(values=shapes4_fill_manual) +
  facet_wrap(~ geneid, scales = "free_y", ncol=3) +
  theme_bw()
```

![](DPCaTSRNA_2_quick_plotgenes_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

### What are DPs

``` r
geneOI <- c("Dll1","Dll3","Neurod4","Nkx2-9",
            "Hey1","Hey2","Neurog1","Neurog2","Hes1","Hes4","Notch1")



#geneOI <- c("Nkx2-9","Nkx6-2","Neurod4","Dll1","Dll3", "Hey1","Neurog1","Neurog2")
ggplot(dds_counts_plot %>% filter(geneid %in% geneOI & Condition %in% c("500","UPSAG","dRA2UPSAG")) %>% 
        mutate(geneid=factor(geneid, levels=geneOI)), 
       aes(x=Gate,y=counts_norm)) +
  stat_summary(aes(fill=Condition),
    fun = mean, geom="bar", alpha=0.9, width=0.7,position=position_dodge(0.7)) +
  geom_point(aes(fill=Condition), alpha=0.6, position = position_dodge(width = 0.7),color="black") +
  #geom_col(position="dodge",aes(fill=DayGate)) +
  scale_fill_manual(values=colors_conditions) +
  scale_color_manual(values=colors_conditions) +
  scale_shape_manual(values=shapes4_fill_manual) +
  facet_grid(geneid ~ Day, scales = "free_y") +
  theme_bw()
```

![](DPCaTSRNA_2_quick_plotgenes_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
ggplot(dds_counts_plot %>% filter(geneid %in% geneOI & Condition %in% c("500","UPSAG","dRA2UPSAG")) %>%
         mutate(geneid=factor(geneid, levels=geneOI)), 
       aes(x=Day,y=counts_norm)) +
  stat_summary(aes(fill=Gate),
    fun = mean, geom="bar", alpha=0.9, width=0.7,position=position_dodge(0.7)) +
  geom_point(aes(fill=Gate), alpha=0.6, position = position_dodge(width = 0.7),color="black") +
  #geom_col(position="dodge",aes(fill=DayGate)) +
  scale_fill_manual(values=color_gates) +
  scale_color_manual(values=color_gates) +
  scale_shape_manual(values=shapes4_fill_manual) +
  facet_grid(geneid ~ Condition, scales = "free_y") +
  theme_bw()
```

![](DPCaTSRNA_2_quick_plotgenes_files/figure-gfm/unnamed-chunk-10-2.png)<!-- -->

``` r
ggplot(dds_counts_plot %>% filter(geneid %in% geneOI & Condition %in% c("500","UPSAG","dRA2UPSAG")) %>%
         filter(Gate !="neur") %>%
         mutate(geneid=factor(geneid, levels=geneOI)), 
       aes(x=Gate,y=counts_norm)) +
  stat_summary(aes(fill=Day),
    fun = mean, geom="bar", alpha=0.9, width=0.7,position=position_dodge(0.7)) +
  geom_point(aes(fill=Day), alpha=0.6, position = position_dodge(width = 0.7),color="black") +
  #geom_col(position="dodge",aes(fill=DayGate)) +
  #scale_fill_manual(values=color_gates) +
  #scale_color_manual(values=color_gates) +
  scale_shape_manual(values=shapes4_fill_manual) +
  facet_grid(geneid ~ Condition, scales = "free_y") +
  theme_bw()
```

![](DPCaTSRNA_2_quick_plotgenes_files/figure-gfm/unnamed-chunk-10-3.png)<!-- -->

### plots for Poster

``` r
geneOI <- c("Nkx2-9","Dll1","Hey1","Neurog2","Hes1","Notch1")


ggplot(dds_counts_plot %>% filter(geneid %in% geneOI & Condition %in% c("500") & Gate!="neurons" ) %>%
         mutate(geneid=factor(geneid, levels=geneOI)), 
       aes(x=Day,y=counts_norm)) +
  stat_summary(aes(fill=Gate),
    fun = mean, geom="bar", alpha=0.9, width=0.7,position=position_dodge(0.7)) +
  geom_point(aes(fill=Gate), alpha=0.6, position = position_dodge(width = 0.7),color="black") +
  #geom_col(position="dodge",aes(fill=DayGate)) +
  scale_fill_manual(values=color_gates) +
  scale_color_manual(values=color_gates) +
  scale_shape_manual(values=shapes4_fill_manual) +
  facet_wrap(~ geneid , scales = "free_y") +
  theme_bw()
```

![](DPCaTSRNA_2_quick_plotgenes_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

Temporal program

``` r
geneOI <- c("Nr6a1","Lef1","Sox9","Rfx4","Nfib","Atf3","Fos")

ggplot(dds_counts_plot %>% filter(geneid %in% geneOI & Condition %in% c("500","UPSAG","dRA2UPSAG")) %>%
         mutate(geneid=factor(geneid, levels=geneOI)), 
       aes(x=Day,y=counts_norm)) +
  stat_summary(aes(fill=Gate),
    fun = mean, geom="bar", alpha=0.9, width=0.7,position=position_dodge(0.7)) +
  geom_point(aes(fill=Gate), alpha=0.6, position = position_dodge(width = 0.7),color="black") +
  #geom_col(position="dodge",aes(fill=DayGate)) +
  scale_fill_manual(values=color_gates) +
  scale_color_manual(values=color_gates) +
  scale_shape_manual(values=shapes4_fill_manual) +
  facet_grid(geneid ~ Condition, scales = "free_y") +
  theme_bw()
```

![](DPCaTSRNA_2_quick_plotgenes_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

``` r
ggplot(dds_counts_plot %>% filter(geneid %in% geneOI & Condition %in% c("500","UPSAG","dRA2UPSAG")) %>%
         filter(Gate !="neur") %>%
         mutate(geneid=factor(geneid, levels=geneOI)), 
       aes(x=Gate,y=counts_norm)) +
  stat_summary(aes(fill=Day),
    fun = mean, geom="bar", alpha=0.9, width=0.7,position=position_dodge(0.7)) +
  geom_point(aes(fill=Day), alpha=0.6, position = position_dodge(width = 0.7),color="black") +
  #geom_col(position="dodge",aes(fill=DayGate)) +
  #scale_fill_manual(values=color_gates) +
  #scale_color_manual(values=color_gates) +
  scale_shape_manual(values=shapes4_fill_manual) +
  facet_grid(geneid ~ Condition, scales = "free_y") +
  theme_bw()
```

![](DPCaTSRNA_2_quick_plotgenes_files/figure-gfm/unnamed-chunk-12-2.png)<!-- -->

Genes from Rayon et al

``` r
geneOI <- c("Olig2","Olig1","Nkx2-2","Nkx2-9","Nkx6-1","Nkx6-2","Fabp7",
            "Sox9","Neurod4","Neurog1","Neurog2","Neurog3")


ggplot(dds_counts_plot %>% filter(geneid %in% geneOI & Condition %in% c("500","UPSAG","dRA2UPSAG")) %>%
         mutate(geneid=factor(geneid, levels=geneOI)), 
       aes(x=Day,y=counts_norm)) +
  stat_summary(aes(fill=Gate),
    fun = mean, geom="bar", alpha=0.9, width=0.7,position=position_dodge(0.7)) +
  geom_point(aes(fill=Gate), alpha=0.6, position = position_dodge(width = 0.7),color="black") +
  #geom_col(position="dodge",aes(fill=DayGate)) +
  scale_fill_manual(values=color_gates) +
  scale_color_manual(values=color_gates) +
  scale_shape_manual(values=shapes4_fill_manual) +
  facet_grid(geneid ~ Condition, scales = "free_y") +
  theme_bw()
```

![](DPCaTSRNA_2_quick_plotgenes_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
ggplot(dds_counts_plot %>% filter(geneid %in% geneOI & Condition %in% c("500","UPSAG","dRA2UPSAG")) %>%
         filter(Gate !="neur") %>%
         mutate(geneid=factor(geneid, levels=geneOI)), 
       aes(x=Gate,y=counts_norm)) +
  stat_summary(aes(fill=Day),
    fun = mean, geom="bar", alpha=0.9, width=0.7,position=position_dodge(0.7)) +
  geom_point(aes(fill=Day), alpha=0.6, position = position_dodge(width = 0.7),color="black") +
  #geom_col(position="dodge",aes(fill=DayGate)) +
  #scale_fill_manual(values=color_gates) +
  #scale_color_manual(values=color_gates) +
  scale_shape_manual(values=shapes4_fill_manual) +
  facet_grid(geneid ~ Condition, scales = "free_y") +
  theme_bw()
```

![](DPCaTSRNA_2_quick_plotgenes_files/figure-gfm/unnamed-chunk-13-2.png)<!-- -->

Rayon’s “gliogenic score”

FABP7, SOX9, SOX10, PDGFRA, CSPG4, FGFR3, FGFBP3, DBI, SLC1A3, HOPX,
ALDH1L1

``` r
geneOI <- c("Fabp7","Sox9","Sox10","Pdgfra","Cspg4","Fgfr3","Fgfbp3","Dbi","Slc1a3","Hopx","Aldh1l1")


ggplot(dds_counts_plot %>% filter(geneid %in% geneOI & Condition %in% c("500","UPSAG","dRA2UPSAG")) %>%
         mutate(geneid=factor(geneid, levels=geneOI)), 
       aes(x=Day,y=counts_norm)) +
  stat_summary(aes(fill=Gate),
    fun = mean, geom="bar", alpha=0.9, width=0.7,position=position_dodge(0.7)) +
  geom_point(aes(fill=Gate), alpha=0.6, position = position_dodge(width = 0.7),color="black") +
  #geom_col(position="dodge",aes(fill=DayGate)) +
  scale_fill_manual(values=color_gates) +
  scale_color_manual(values=color_gates) +
  scale_shape_manual(values=shapes4_fill_manual) +
  facet_grid(geneid ~ Condition, scales = "free_y") +
  theme_bw()
```

![](DPCaTSRNA_2_quick_plotgenes_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

``` r
ggplot(dds_counts_plot %>% filter(geneid %in% geneOI & Condition %in% c("500","UPSAG","dRA2UPSAG")) %>%
         filter(Gate !="neur") %>%
         mutate(geneid=factor(geneid, levels=geneOI)), 
       aes(x=Gate,y=counts_norm)) +
  stat_summary(aes(fill=Day),
    fun = mean, geom="bar", alpha=0.9, width=0.7,position=position_dodge(0.7)) +
  geom_point(aes(fill=Day), alpha=0.6, position = position_dodge(width = 0.7),color="black") +
  #geom_col(position="dodge",aes(fill=DayGate)) +
  #scale_fill_manual(values=color_gates) +
  #scale_color_manual(values=color_gates) +
  scale_shape_manual(values=shapes4_fill_manual) +
  facet_grid(geneid ~ Condition, scales = "free_y") +
  theme_bw()
```

![](DPCaTSRNA_2_quick_plotgenes_files/figure-gfm/unnamed-chunk-14-2.png)<!-- -->

### Wnt genes

following Agalliu et al

``` r
geneOI <- c("Wnt4","Wnt5a","Wnt5b","Wnt1","Wnt7a","Wnt7b")


ggplot(dds_counts_plot %>% filter(geneid %in% geneOI & Condition %in% c("500","UPSAG","dRA2UPSAG")) %>%
         mutate(geneid=factor(geneid, levels=geneOI)), 
       aes(x=Day,y=counts_norm)) +
  stat_summary(aes(fill=Gate),
    fun = mean, geom="bar", alpha=0.9, width=0.7,position=position_dodge(0.7)) +
  geom_point(aes(fill=Gate), alpha=0.6, position = position_dodge(width = 0.7),color="black") +
  #geom_col(position="dodge",aes(fill=DayGate)) +
  scale_fill_manual(values=color_gates) +
  scale_color_manual(values=color_gates) +
  scale_shape_manual(values=shapes4_fill_manual) +
  facet_grid(geneid ~ Condition, scales = "free_y") +
  theme_bw()
```

![](DPCaTSRNA_2_quick_plotgenes_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
ggplot(dds_counts_plot %>% filter(geneid %in% geneOI & Condition %in% c("500","UPSAG","dRA2UPSAG")) %>%
         filter(Gate !="neur") %>%
         mutate(geneid=factor(geneid, levels=geneOI)), 
       aes(x=Gate,y=counts_norm)) +
  stat_summary(aes(fill=Day),
    fun = mean, geom="bar", alpha=0.9, width=0.7,position=position_dodge(0.7)) +
  geom_point(aes(fill=Day), alpha=0.6, position = position_dodge(width = 0.7),color="black") +
  #geom_col(position="dodge",aes(fill=DayGate)) +
  #scale_fill_manual(values=color_gates) +
  #scale_color_manual(values=color_gates) +
  scale_shape_manual(values=shapes4_fill_manual) +
  facet_grid(geneid ~ Condition, scales = "free_y") +
  theme_bw()
```

![](DPCaTSRNA_2_quick_plotgenes_files/figure-gfm/unnamed-chunk-15-2.png)<!-- -->

``` r
geneOI <- c("Nkx2-9","Sp8","Foxa2","Phox2b")


ggplot(dds_counts_plot %>% filter(geneid %in% geneOI & Condition %in% c("500","UPSAG","dRA2UPSAG")) %>%
         mutate(geneid=factor(geneid, levels=geneOI)), 
       aes(x=Day,y=counts_norm)) +
  stat_summary(aes(fill=Gate),
    fun = mean, geom="bar", alpha=0.9, width=0.7,position=position_dodge(0.7)) +
  geom_point(aes(fill=Gate), alpha=0.6, position = position_dodge(width = 0.7),color="black") +
  #geom_col(position="dodge",aes(fill=DayGate)) +
  scale_fill_manual(values=color_gates) +
  scale_color_manual(values=color_gates) +
  scale_shape_manual(values=shapes4_fill_manual) +
  facet_grid(geneid ~ Condition, scales = "free_y") +
  theme_bw()
```

![](DPCaTSRNA_2_quick_plotgenes_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

``` r
ggplot(dds_counts_plot %>% filter(geneid %in% geneOI & Condition %in% c("500","UPSAG","dRA2UPSAG")) %>%
         filter(Gate !="neur") %>%
         mutate(geneid=factor(geneid, levels=geneOI)), 
       aes(x=Gate,y=counts_norm)) +
  stat_summary(aes(fill=Day),
    fun = mean, geom="bar", alpha=0.9, width=0.7,position=position_dodge(0.7)) +
  geom_point(aes(fill=Day), alpha=0.6, position = position_dodge(width = 0.7),color="black") +
  #geom_col(position="dodge",aes(fill=DayGate)) +
  #scale_fill_manual(values=color_gates) +
  #scale_color_manual(values=color_gates) +
  scale_shape_manual(values=shapes4_fill_manual) +
  facet_grid(geneid ~ Condition, scales = "free_y") +
  theme_bw()
```

![](DPCaTSRNA_2_quick_plotgenes_files/figure-gfm/unnamed-chunk-16-2.png)<!-- -->

### Genes for Marine

``` r
geneOI <- c("Olig2","Nkx2-2","Nkx2-9","Pappa","Sox9","Sox6","Nav2","Gpm6b")


ggplot(dds_counts_plot %>% filter(geneid %in% geneOI & Condition %in% c("500","UPSAG")) %>%
         mutate(geneid=factor(geneid, levels=geneOI)), 
       aes(x=Day,y=counts_norm, fill=Gate)) +
  stat_summary(aes(fill=Gate),
    fun = mean, geom="bar", alpha=0.9, width=0.7,position=position_dodge2(0.7,preserve = "single")) +
  geom_point(aes(fill=Gate), alpha=0.6, position = position_dodge2(width = 0.7,preserve = "single"),color="black", shape=21) +
  #geom_col(position="dodge",aes(fill=DayGate)) +
  scale_fill_manual(values=color_gates) +
  scale_color_manual(values=color_gates) +
  scale_shape_manual(values=shapes4_fill_manual) +
  facet_grid(geneid ~ Condition, scales = "free") +
  theme_bw()
```

![](DPCaTSRNA_2_quick_plotgenes_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

``` r
ggplot(dds_counts_plot %>% filter(geneid %in% geneOI & Condition %in% c("500","UPSAG")) %>%
         filter(Gate !="neur") %>%
         mutate(geneid=factor(geneid, levels=geneOI)), 
       aes(x=Gate,y=counts_norm)) +
  stat_summary(aes(fill=Day),
    fun = mean, geom="bar", alpha=0.9, width=0.7,position=position_dodge2(0.7,preserve = "single")) +
  geom_point(aes(fill=Day), alpha=0.6, position = position_dodge2(0.7,preserve = "single"),color="black", shape=21) +
  #geom_col(position="dodge",aes(fill=DayGate)) +
  #scale_fill_manual(values=color_gates) +
  #scale_color_manual(values=color_gates) +
  scale_shape_manual(values=shapes4_fill_manual) +
  facet_grid(geneid ~ Condition, scales = "free_y") +
  theme_bw()
```

![](DPCaTSRNA_2_quick_plotgenes_files/figure-gfm/unnamed-chunk-17-2.png)<!-- -->

## genes for James

``` r
geneOI <- c("Tle1","Tle2","Tle3","Tle4","Tle6","Aes")


ggplot(dds_counts_plot %>% filter(geneid %in% geneOI) %>%
         mutate(geneid=factor(geneid, levels=geneOI)), 
       aes(x=Day,y=counts_norm, fill=Gate)) +
  stat_summary(aes(fill=Gate),
    fun = mean, geom="bar", alpha=0.9, width=0.7,position=position_dodge2(0.7,preserve = "single")) +
  geom_point(aes(fill=Gate), alpha=0.6, position = position_dodge2(width = 0.7,preserve = "single"),color="black", shape=21) +
  #geom_col(position="dodge",aes(fill=DayGate)) +
  scale_fill_manual(values=color_gates) +
  scale_color_manual(values=color_gates) +
  scale_shape_manual(values=shapes4_fill_manual) +
  facet_grid(geneid ~ Condition, scales = "free") +
  theme_bw()
```

![](DPCaTSRNA_2_quick_plotgenes_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

``` r
ggplot(dds_counts_plot %>% filter(geneid %in% geneOI & Condition %in% c("500","UPSAG")) %>%
         filter(Gate !="neur") %>%
         mutate(geneid=factor(geneid, levels=geneOI)), 
       aes(x=Gate,y=counts_norm)) +
  stat_summary(aes(fill=Day),
    fun = mean, geom="bar", alpha=0.9, width=0.7,position=position_dodge2(0.7,preserve = "single")) +
  geom_point(aes(fill=Day), alpha=0.6, position = position_dodge2(0.7,preserve = "single"),color="black", shape=21) +
  #geom_col(position="dodge",aes(fill=DayGate)) +
  #scale_fill_manual(values=color_gates) +
  #scale_color_manual(values=color_gates) +
  scale_shape_manual(values=shapes4_fill_manual) +
  facet_grid(geneid ~ Condition, scales = "free_y") +
  theme_bw()
```

![](DPCaTSRNA_2_quick_plotgenes_files/figure-gfm/unnamed-chunk-18-2.png)<!-- -->

``` r
sessionInfo()
```

    ## R version 4.4.0 (2024-04-24)
    ## Platform: aarch64-apple-darwin20
    ## Running under: macOS Sonoma 14.6.1
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRblas.0.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: Europe/London
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] grid      stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] ComplexHeatmap_2.20.0 lubridate_1.9.3       forcats_1.0.0        
    ##  [4] stringr_1.5.1         dplyr_1.1.4           purrr_1.0.2          
    ##  [7] readr_2.1.5           tidyr_1.3.1           tibble_3.2.1         
    ## [10] ggplot2_3.5.1         tidyverse_2.0.0       RColorBrewer_1.1-3   
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] utf8_1.2.4          generics_0.1.3      shape_1.4.6.1      
    ##  [4] stringi_1.8.4       hms_1.1.3           digest_0.6.35      
    ##  [7] magrittr_2.0.3      evaluate_0.23       timechange_0.3.0   
    ## [10] iterators_1.0.14    circlize_0.4.16     fastmap_1.2.0      
    ## [13] foreach_1.5.2       doParallel_1.0.17   GlobalOptions_0.1.2
    ## [16] fansi_1.0.6         scales_1.3.0        codetools_0.2-20   
    ## [19] cli_3.6.2           rlang_1.1.4         crayon_1.5.2       
    ## [22] munsell_0.5.1       withr_3.0.0         yaml_2.3.8         
    ## [25] tools_4.4.0         parallel_4.4.0      tzdb_0.4.0         
    ## [28] colorspace_2.1-0    BiocGenerics_0.50.0 GetoptLong_1.0.5   
    ## [31] vctrs_0.6.5         R6_2.5.1            png_0.1-8          
    ## [34] stats4_4.4.0        matrixStats_1.3.0   lifecycle_1.0.4    
    ## [37] S4Vectors_0.42.0    IRanges_2.38.0      clue_0.3-65        
    ## [40] cluster_2.1.6       pkgconfig_2.0.3     pillar_1.9.0       
    ## [43] gtable_0.3.5        glue_1.7.0          highr_0.11         
    ## [46] xfun_0.44           tidyselect_1.2.1    rstudioapi_0.16.0  
    ## [49] knitr_1.47          farver_2.1.2        rjson_0.2.21       
    ## [52] htmltools_0.5.8.1   labeling_0.4.3      rmarkdown_2.27     
    ## [55] compiler_4.4.0

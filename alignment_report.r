#-------------------------------- library ------------------------------------#
library(dplyr)
library(stringr)
library(ggplot2)
library(viridis)
library(GenomicRanges)
library(chromVAR) ## For FRiP analysis and differential analysis
library(ggpubr) ## For customizing figures
library(corrplot) ## For correlation plot

#---------------------------------- Prepare ----------------------------------#
args = commandArgs(trailingOnly = TRUE)

all_file = c(args)
sampleid = c()
sample_group = c()
align_hg38 = c()
align_spikeIn = c()
rmdup = c()
assess_frag = c()
assess_useBED = c()
seqDepth = c()
for (i in all_file) {
  sampleid = unique(append(sampleid, strsplit(i, split = "_")[[1]][2]))
  
  sample_group = unique(append(
    sample_group,
    paste0(strsplit(i, split = "_")[[1]][2], "_group", strsplit(i, split = "_")[[1]][1])
    ))
  
  if(str_detect(i, "hg38")) {
    align_hg38 = append(align_hg38, i)
  } else if (str_detect(i, "spikeIn.txt")) {
    align_spikeIn = append(align_spikeIn, i)
  } else if (str_detect(i, "rmDup")) {
    rmdup = append(rmdup, i)
  } else if (str_detect(i, "fragmentLen")) {
    assess_frag = append(assess_frag, i)
  } else if (str_detect(i, "fragmentsCount")) {
    assess_useBED = append(assess_useBED, i)
  } else if (str_detect(i, "seqDepth")) {
    seqDepth = append(seqDepth, i)
  }
}
sample_group = sort(sample_group)
sampleid = sort(sampleid)
#----------------------------- Sequencing depth -------------------------------#
alignResult = c()

for (summaryfile in align_hg38) {
  alignRes = read.table(summaryfile, header = FALSE, fill = TRUE)

  alignRate = substr(alignRes$V1[6], 1, nchar(as.character(alignRes$V1[6]))-1)
  histInfo= strsplit(summaryfile, split = "_")[[1]]

  alignResult = data.frame(
    SampleId = histInfo[2],
    Group = histInfo[1],
    SequencingDepth = alignRes$V1[1] %>% as.character %>% as.numeric,
    MappedFragNum_hg38 = alignRes$V1[4] %>% as.character %>% as.numeric + alignRes$V1[5] %>% as.character %>% as.numeric,
    AlignmentRate_hg38 = alignRate %>% as.numeric
  ) %>% rbind(alignResult, .)
}

alignResult$SampleId = factor(alignResult$SampleId, levels = sampleid)
alignResult %>% mutate(AlignmentRate_hg38 = paste0(AlignmentRate_hg38, "%"))

#---------------------------- Spike-in alignment -----------------------------#
spikeAlign = c()

for (summaryfile in align_spikeIn) {
  spikeRes = read.table(summaryfile, header = FALSE, fill = TRUE)

  alignRate = substr(spikeRes$V1[6], 1, nchar(as.character(spikeRes$V1[6]))-1)
  histInfo= strsplit(summaryfile, split = "_")[[1]]

  spikeAlign = data.frame(
    SampleId = histInfo[2],
    Group = histInfo[1],
    SequencingDepth = spikeRes$V1[1] %>% as.character %>% as.numeric,
    MappedFragNum_spikeIn = spikeRes$V1[4] %>% as.character %>% as.numeric + spikeRes$V1[5] %>% as.character %>% as.numeric,
    AlignmentRate_spikeIn = alignRate %>% as.numeric
  ) %>% rbind(spikeAlign, .)
}

spikeAlign$SampleId = factor(spikeAlign$SampleId, levels = sampleid)
spikeAlign %>% mutate(AlignmentRate_spikeIn = paste0(AlignmentRate_spikeIn, "%"))

#------------------ Summarize the alignment to hg38 and E.coli --------------#
alignSummary = left_join(alignResult, spikeAlign, by = c("SampleId", "Group", "SequencingDepth")) %>% mutate (
  AlignmentRate_hg38 = paste0(AlignmentRate_hg38, "%"),
  AlignmentRate_spikeIn = paste0(AlignmentRate_spikeIn, "%")
)

#----------- Visualizing the sequencing depth and alignment result -----------#
fig1A = alignResult %>% ggplot(aes(x = SampleId, y = SequencingDepth/1000000, fill = SampleId)) +
  geom_boxplot() +
  geom_jitter(aes(color = Group), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("Sequencing Depth per Million") +
  xlab("") +
  ggtitle("A. Sequencing Depth")

fig1B = alignResult %>% ggplot(aes(x = SampleId, y = MappedFragNum_hg38/1000000, fill = SampleId)) +
  geom_boxplot() +
  geom_jitter(aes(color = Group), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("Mapped Fragments per Million") +
  xlab("") +
  ggtitle("B. Alignable Fragment (hg38)")

fig1C = alignResult %>% ggplot(aes(x = SampleId, y = AlignmentRate_hg38, fill = SampleId)) +
  geom_boxplot() +
  geom_jitter(aes(color = Group), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("% of Mapped Fragments") +
  xlab("") +
  ggtitle("C. Alignment Rate (hg38)")

fig1D = spikeAlign %>% ggplot(aes(x = SampleId, y = AlignmentRate_spikeIn, fill = SampleId)) +
  geom_boxplot() +
  geom_jitter(aes(color = Group), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("Spike-in Alignment Rate") +
  xlab("") +
  ggtitle("D. Alignment Rate (E.coli)")

png(file = "Sequencing_Depth_Summary.png", width = 1200, height = 900)
ggarrange(fig1A, fig1B, fig1C, fig1D, ncol = 2, nrow=2, common.legend = TRUE, legend="bottom")
dev.off()

#------------------------- remove duplicates summary -------------------------#
dupResult = c()

for (summaryfile in rmdup) {
  dupRes = read.table(summaryfile, header = TRUE, fill = TRUE)

  histInfo= strsplit(summaryfile, split = "_")[[1]]

  dupResult = data.frame(
    SampleId = histInfo[2],
    Group = histInfo[1],
    MappedFragNum_hg38 = dupRes$READ_PAIRS_EXAMINED[1] %>% as.character %>% as.numeric,
    DuplicationRate = dupRes$PERCENT_DUPLICATION[1] %>% as.character %>% as.numeric * 100, 
    EstimatedLibrarySize = dupRes$ESTIMATED_LIBRARY_SIZE[1] %>% as.character %>% as.numeric) %>% mutate(UniqueFragNum = MappedFragNum_hg38 * (1-DuplicationRate/100))  %>% rbind(dupResult, .)
}

dupResult$SampleId = factor(dupResult$SampleId, levels = sampleid)
alignDupSummary = left_join(
  alignSummary,
  dupResult,
  by = c("SampleId", "Group", "MappedFragNum_hg38")
  ) %>% mutate(DuplicationRate = paste0(DuplicationRate, "%"))

#--------------------- Visualizing the duplication rate ----------------------#
fig2A = dupResult %>% ggplot(aes(x = SampleId, y = DuplicationRate, fill = SampleId)) +
  geom_boxplot() +
  geom_jitter(aes(color = Group), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("Duplication Rate (*100%)") +
  xlab("")

fig2B = dupResult %>% ggplot(aes(x = SampleId, y = EstimatedLibrarySize, fill = SampleId)) +
  geom_boxplot() +
  geom_jitter(aes(color = Group), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("Estimated Library Size") +
  xlab("")

fig2C = dupResult %>% ggplot(aes(x = SampleId, y = UniqueFragNum, fill = SampleId)) +
  geom_boxplot() +
  geom_jitter(aes(color = Group), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("# of Unique Fragments") +
  xlab("")

png(file = "Duplicate_Summary.png", width = 1200, height = 900)
ggarrange(fig2A, fig2B, fig2C, ncol = 3, common.legend = TRUE, legend="bottom")
dev.off()

#--------------- Assess mapped fragment size distribution --------------------#
fragLen = c()

for (summaryfile in assess_frag) {
  histInfo = strsplit(summaryfile, "_")[[1]]

  fragLen = read.table(summaryfile, header = FALSE) %>% mutate(
    fragLen = V1 %>% as.numeric,
    fragCount = V2 %>% as.numeric,
    Weight = as.numeric(V2) / sum(as.numeric(V2)),
    SampleId = histInfo[2],
    Group = histInfo[1],
    SampleInfo = paste0(SampleId, "_group", Group)
    ) %>% rbind(fragLen, .)
}

fragLen$SampleInfo = factor(fragLen$SampleInfo, levels = sample_group)
fragLen$SampleId = factor(fragLen$SampleId, levels = sampleid)

#--------------- Visualizing the fragment size density -----------------------#
fig3A = fragLen %>% ggplot(aes(x = SampleInfo, y = fragLen, weight = Weight, fill = SampleId)) +
  geom_violin(bw = 5) +
  scale_y_continuous(breaks = seq(0, 800, 50)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 20) +
  ggpubr::rotate_x_text(angle = 20) +
  ylab("Fragment Length") +
  xlab("")

fig3B = fragLen %>% ggplot(aes(x = fragLen, y = fragCount, color = SampleId, group = SampleInfo, linetype = Group)) +
  geom_line(size = 1) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma") +
  theme_bw(base_size = 20) +
  xlab("Fragment Length") +
  ylab("Count") +
  coord_cartesian(xlim = c(0, 500))

png(file = "fragment_size_summary.png", width = 1200, height = 600)
  ggarrange(fig3A, fig3B, ncol = 2)
dev.off()

#-------------------- Assess replicate reproducibility -----------------------#
reprod = c()

fragCount = NULL

for(summaryfile in assess_useBED) {
  histInfo = strsplit(summaryfile, "_")[[1]]
  SampleInfo = paste0(histInfo[2], "_group", histInfo[1])
  
  if(is.null(fragCount)){
    fragCount = read.table(summaryfile, header = FALSE)
    colnames(fragCount) = c("chrom", "bin", SampleInfo)
  }else{
    fragCountTmp = read.table(summaryfile, header = FALSE)
    colnames(fragCountTmp) = c("chrom", "bin", SampleInfo)
    fragCount = full_join(fragCount, fragCountTmp, by = c("chrom", "bin"))
  }
}

M = cor(fragCount %>% select(-c("chrom", "bin")) %>% log2(), use = "complete.obs")

png(file = "Replicate_Reproducibility.png", width = 1200, height = 900)
corrplot(
  M, method = "color", outline = T, addgrid.col = "darkgray", order="hclust",
  addrect = 3, rect.col = "black", rect.lwd = 3,cl.pos = "b",
  tl.col = "indianred4", tl.cex = 1, cl.cex = 1, addCoef.col = "black",
  number.digits = 2, number.cex = 1,
  col = colorRampPalette(c("midnightblue","white","darkred"))(100))
dev.off()

#------------------------------- Scaling factor ------------------------------#
scaleFactor = c()

multiplier = 10000

for(summaryfile in seqDepth) {
  spikeDepth = read.table(summaryfile, header = FALSE, fill = TRUE)$V1[1]
  histInfo = strsplit(summaryfile, "_")[[1]]

  scaleFactor = data.frame(
    scaleFactor = multiplier/spikeDepth,
    SampleId = histInfo[2], Group = histInfo[1]
  )  %>% rbind(scaleFactor, .)
}

scaleFactor$SampleId = factor(scaleFactor$SampleId, levels = sampleid)
alignSFsummary = left_join(alignDupSummary, scaleFactor, by = c("SampleId", "Group"))

write.csv(
  alignSFsummary[order(alignSFsummary$SampleId, alignSFsummary$Group), ],
  file = "Alignment_summary.csv", row.names = FALSE
)

#--------------------- Visualizing the scaling factor ------------------------#
fig5A = scaleFactor %>% ggplot(aes(x = SampleId, y = scaleFactor, fill = SampleId)) +
  geom_boxplot() +
  geom_jitter(aes(color = Group), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 20) +
  ylab("Spike-in Scalling Factor") +
  xlab("")

normDepth = inner_join(scaleFactor, alignResult, by = c("SampleId", "Group")) %>% mutate(normDepth = MappedFragNum_hg38 * scaleFactor)

fig5B = normDepth %>% ggplot(aes(x = SampleId, y = normDepth, fill = SampleId)) +
  geom_boxplot() +
  geom_jitter(aes(color = Group), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 20) +
  ylab("Normalization Fragment Count") +
  xlab("") +
  coord_cartesian(ylim = c(1000000, 130000000))

png(file = "Sequencing_depth.png", width = 1200, height = 600)
ggarrange(fig5A, fig5B, ncol = 2, common.legend = TRUE, legend="bottom")
dev.off()



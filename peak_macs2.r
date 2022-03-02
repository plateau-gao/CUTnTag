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
peakfile = c()
bamfile = c()
peaktype = c()
sampleid = c()
group = c()
alignResult = c()

for (file in all_file) {
    if (str_detect(file, "narrowPeak")) {
        peakfile = append(peakfile, file)
        peaktype = unique(append(peaktype, strsplit(file, split="_")[[1]][4]))
        sampleid = unique(append(sampleid, strsplit(file, split = "_")[[1]][2]))
        group = unique(append(group, strsplit(file, split = "_")[[1]][1]))
    } else if (str_detect(file, "mapped.bam")) {
        bamfile = append(bamfile, file)
    } else if (str_detect(file, "summary.csv")) {
        alignResult = read.csv(file, header = TRUE)
    }
}

group = sort(group)
sampleid = sort(sampleid)
peaktype = sort(peaktype)

alignResult = select(alignResult, SampleId, Group, SequencingDepth, MappedFragNum_hg38, AlignmentRate_hg38)
alignResult$Group = as.integer(alignResult$Group)
alignResult$SampleId = factor(alignResult$SampleId, levels = sampleid)

#-------------------------- Number of peaks called ---------------------------#
peakN = c()
peakWidth = c()
for (file in peakfile) {
    histInfo = strsplit(file, "_")[[1]]
    peakInfo = read.table(file, header = FALSE, fill = TRUE) %>% mutate(width = abs(V3-V2))    
    peakN = data.frame(peakN = nrow(peakInfo), peakType = histInfo[4], SampleId = histInfo[2], Group = histInfo[1]) %>% rbind(peakN, .)
    peakWidth = data.frame(width = peakInfo$width, peakType = histInfo[4], SampleId = histInfo[2], Group = histInfo[1]) %>% rbind(peakWidth, .)
}

peakN %>% select(SampleId, Group, peakType, peakN)
peakN$SampleId = factor(peakN$SampleId, levels = sampleid)
peakN$Group = factor(peakN$Group, levels = group)

#---------- Reproducibility of the peak across biological replicates ---------#
peakOverlap = c()
for (type in peaktype) {
    for (id in sampleid) {
        overlap.gr = GRanges()
        for (rep in group) {
            peakInfo = read.table(str_subset(peakfile, paste0(rep, "_", id, "_macs2_", type)), header = FALSE, fill = TRUE)
            peakInfo.gr = GRanges(peakInfo$V1, IRanges(start = peakInfo$V2, end = peakInfo$V3), strand = "*")
            if (length(overlap.gr) > 0) {
                overlap.gr = overlap.gr[findOverlaps(overlap.gr, peakInfo.gr)@from]
            } else {
                overlap.gr = peakInfo.gr
            }
        }
        peakOverlap = data.frame(peakReprod = length(overlap.gr), SampleId = id, peakType = type) %>% rbind(peakOverlap, .) 
    }
}

peakReprod = left_join(peakN, peakOverlap, by = c("SampleId", "peakType")) %>% mutate(peakReprodRate = peakReprod / peakN * 100)
peak_summary <- peakReprod %>% select(SampleId, Group, peakType, peakN, peakReprodNum = peakReprod, peakReprodRate)

write.csv(peak_summary[order(peakN$SampleId, peakN$Group, peakN$peakType), ], file = paste0("macs2_peak_summary.csv"), row.names = FALSE)

#--------------- FRagment proportion in Peaks regions ----------------------#
inPeakData = c()

for (id in sampleid) {
    for (rep in group) {
        peakRes = read.table(str_subset(peakfile, paste0(rep, "_", id, "_macs2_baseq")), header = FALSE, fill = TRUE)
        peak.gr = GRanges(seqnames = peakRes$V1, IRanges(start = peakRes$V2, end = peakRes$V3), strand = "*")
        bamfile = paste0(rep, "_", id, "_mapped.bam")
        fragment_counts = getCounts(bamfile, peak.gr, paired = TRUE, by_rg = FALSE, format = "bam")
        inPeakN = counts(fragment_counts)[,1] %>% sum
        inPeakData = rbind(inPeakData, data.frame(inPeakN = inPeakN, SampleId = id, Group = rep))
    }
}

inPeakData$Group = as.integer(inPeakData$Group)

frip = left_join(inPeakData, alignResult, by = c("SampleId", "Group")) %>% mutate(frip = inPeakN/MappedFragNum_hg38 *100)
frip$SampleId = factor(frip$SampleId, levels = sampleid)
frip$Group = factor(frip$Group, levels = group)
frip_summary <- frip %>% select(SampleId, Group, SequencingDepth, MappedFragNum_hg38, AlignmentRate_hg38, FragInPeakNum = inPeakN, FRiPs = frip)

write.csv(frip_summary[order(frip_summary$SampleId, frip_summary$Group), ], file = paste0("macs2_frip_summary.csv"), row.names = FALSE) 

#------------------------- Virualization ----------------------------------#
fig6A = peakN %>% ggplot(aes(x = as.factor(SampleId), y = peakN, fill = SampleId)) + 
    geom_boxplot() +
    geom_jitter(aes(color = Group), position = position_jitter(0.15)) +
    facet_grid(~peakType) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 18) +
    ylab("Number of Peaks") +
    xlab("")

fig6B = peakWidth %>% ggplot(aes(x = as.factor(SampleId), y = width, fill = SampleId)) +
    geom_violin() +
    facet_grid(Group~peakType) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    scale_y_continuous(trans = "log", breaks = c(400, 3000, 22000)) +
    theme_bw(base_size = 18) +
    ylab("Width of Peaks") +
    xlab("")

fig6C = peakReprod %>% ggplot(aes(x = as.factor(SampleId), y = peakReprodRate, fill = SampleId, label = round(peakReprodRate, 2))) +
    geom_bar(stat = "identity") +
    geom_text(vjust = 0.1) +
    facet_grid(Group~peakType) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 18) +
    ylab("% of Peaks Reproduced") +
    xlab("")

fig6D = frip %>% ggplot(aes(x = as.factor(SampleId), y = frip, fill = SampleId, label = round(frip, 2))) +
    geom_boxplot() +
    geom_jitter(aes(color = Group), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 18) +
    ylab("% of Fragments in Peaks") +
    xlab("")
png(file = paste0("macs2_frip_summary.png"), width = 1200, height = 900)
ggarrange(fig6A, fig6B, fig6C, fig6D, ncol = 2, nrow=2, common.legend = TRUE, legend="bottom")
dev.off()



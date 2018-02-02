################################################################################
#
#                                   META
#
################################################################################

########################################
# Arguments
########################################

# Mapping quality
MQ <- 30
# Path
pathCG <- "../genome/CG.RDS"
pathBED <- "./"

########################################
# Sources
########################################

library(data.table)
library(foreach)
library(ggplot2)
library(dplyr)
library(fst)


################################################################################
#
#                                 FUNCTIONCS
#
################################################################################


########################################
# Write BED 
########################################

writeBed <- function(data, name) {
    library(data.table)
    setDT(data)
    write.table(setkey(data, chr, start, end), 
                paste0(name, ".bed"), quote = FALSE, sep = "\t", 
                col.names = FALSE, row.names = FALSE)
}





################################################################################
#
#                                 GET DATA
#
################################################################################


########################################
# CG
########################################

dCG <- setDT(readRDS(pathCG))
if (!file.exists("CG.bed")) {
    writeBed(dCG, "CG")
}


########################################
# Load BED files
########################################

files <- list.files(pattern = ".*Ion.*bed")

dCoverage <- foreach(i = files[1:2], .combine = cbind) %do% {
    ID <- paste0("Sample_", strsplit(i, "_")[[1]][5])
    print(ID)
    dOrig <- fread(i) %>%
        setnames(c("chr", "start", "end", "ID", "mq", "strand"))
    saveRDS(dOrig[, .N, .(chr, strand, mq)], 
            paste0("stats_MQ_", ID, ".RDS"), compress = FALSE)
    dOrig %>%
        .[mq > MQ] %>%
        .[, chr := sub("Chr_", "chr", chr)] %>%
        .[grep("chr", chr)] %>%
        # PADARYTI A
        # .[, .(chr, start = ifelse(strand == "+", start, end), strand, ID)] %>%
        .[, .(chr, start = ifelse(strand == "+", start, end), strand)] %>%
        .[, end := start + 1] %>%
        .[, .(chr, start, end, strand, ID)] %>%
        .[, ID := ID] %>%
        writeBed("starts")
    dDistance <- fread("bedtools closest -a starts.bed -b CG.bed -d")
    saveRDS(dDistance[, .N, .(strand = V4, distance = V9)], 
        paste0("stats_Distance_", ID, ".RDS"), compress = FALSE)
    foo <- dDistance[V9 <= 0, .N, .(chr = V6, start = V7)]
    merge(dCG, foo, c("chr", "start"), all.x = TRUE)[is.na(N), N := 0]$N
}





################################################################################
#
#                              MAPPING STATS
#
################################################################################


########################################
# MQ
########################################

files <- list.files(pattern = ".*MQ.*RDS")

dMQ <- foreach(i = files, .combine = rbind) %do% {
    readRDS(i)
}

# ggplot(dMQ, aes(quality, color = ID)) +
#     geom_density() +
#     labs(x = "Mapping Quality",
#          title = "Mapping Quality (sample ID)") +
#     theme_classic() +
#     scale_fill_brewer(palette = "Dark2")


# ggplot(dMQ, aes(quality, color = ID)) +
#     labs(x = "Mapping Quality",
#          title = "Mapping Quality (sample ID)") + 
#     theme_classic() +
#     geom_histogram()

# ggplot(dMQ, aes(quality, fill = ID)) +
#     labs(x = "Mapping Quality",
#          title = "Mapping Quality (sample ID)") + 
#     theme_classic() +
#     geom_bar()

# ggplot(dMQ, aes(ID, quality, fill = ID)) +
#     geom_boxplot() + 
#     labs(x = "Sample ID",
#          y = "Mapping Quality",
#          title = "Mapping Quality (sample ID)") +
#     theme_classic() +
#     scale_fill_brewer(palette = "Dark2")


# p <- ggplot(dMQ, aes(strand, quality, fill = ID)) +
#     geom_boxplot() +
#     labs(y = "Mapping Quality",
#          x = "Strand", 
#          title = "Mapping Quality (sample ID, strand)") +
#     theme_classic() +
#     scale_fill_brewer(palette = "Dark2")
# ggsave("~/tmp.pdf", p)

# ggplot(dMQ, aes(chr, quality, fill = ID)) +
#     geom_boxplot() +
#     labs(y = "Mapping Quality",
#          x = "Chromosome", 
#          title = "Mapping Quality (sample ID, chromosome)") +
#     theme_classic() +
#     scale_fill_brewer(palette = "Dark2")



########################################
# Number of reads
########################################

pd <- dStarts[, .N, ID]
ggplot(pd, aes(ID, N)) +
    geom_bar(stat = "identity", position = "dodge") +
labs(x = "Sample ID",
     y = "Number of reads",
     title = "Number of reads (sample ID)") +
    theme_classic() +
    scale_fill_brewer(palette = "Dark2")

pchr <- dStarts[, .N, chr]
ggplot(pchr, aes(chr, N)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "Chromosome",
         y = "Number of reads",
         title = "Number of reads (chromosome)") +
    theme_classic() +
    scale_fill_brewer(palette = "Dark2")

pstrand <- dStarts[, .N, strand]
ggplot(pstrand, aes (strand, N)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "Strand",
         y = "Number of reads",
         title = "Number of reads (strand)") +
    theme_classic() +
    scale_fill_brewer(palette = "Dark2")


################################################################################
#
#                               MAP TO CG
#
################################################################################

dCG <- fread(pathCG)
# Distance from read (for each read) start to CG

write.table(dStarts, "starts.bed",
          quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(dCG, "CG.bed",
          quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

system("bedtools sort -i CG.bed > CG2.bed")

dDistance <- fread("bedtools closest -a starts.bed -b CG2.bed -d") %>%
    .[V6 != -1] %>%
    setnames(c(colnames(dStarts), colnames(dCG), "distance"))

########################################
# Distance to CG
########################################
# per strand

#Signal (coverage) statistics 2.1 Correlation 2.1.1 Values 2.1.2 Heatmap 2.1.3 Scatter plot 2.2 Profile per genome

ggplot(dDistance, aes(x = distance)) +
    labs(x = "Distance",
         y = "Number of reads",
         title = "Distance to CG") +
    theme_classic() +
    geom_histogram()

ggplot(dDistance, aes(distance, fill = strand)) +
    labs(x = "Distance",
         y = "Number of reads",
         title = "Distance to CG") +
    geom_bar()

ggplot(dDistance, aes(distance, fill = chr)) +
    labs(x = "Distance",
         y = "Number of reads",
         title = "Distance to CG") +
    geom_bar()

ggplot(dDistance, aes(distance, color = strand)) +
    labs(x = "Distance",
         y = "Density",
         title = "Distance to CG") +
    geom_density() +
    theme_classic() +
    scale_fill_brewer(palette = "Dark2")

ggplot(dDistance, aes(log(distance), color = strand)) +
    geom_density() +
    theme_classic() +
    scale_fill_brewer(palette = "Dark2")

################################################################################
#
#                               COVERAGE PER CG
#
################################################################################


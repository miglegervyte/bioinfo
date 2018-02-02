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
pathCG <- "~/Desktop/bioinfo/Coprinopsis/naujas/bed_html/CG.bed"
pathBED <- "./"

########################################
# Sources
########################################

library(data.table)
library(ggplot2)
library(dplyr)
library(fst)



################################################################################
#
#                                 GET DATA
#
################################################################################

########################################
# Load BED files
########################################

setwd("~/Desktop/bioinfo/Coprinopsis/naujas/bed_html")

files <- list.files(pattern = "bed$")

# for(i in files) {
library(foreach)
dStarts <- foreach(i = files, .combine = rbind) %do% {
    ID <- sub("_.*", "", sub(".*_0", "S", i))
    d <- i %>%
        fread(nrows = 100) %>%
            .[, "V4" := NULL] %>%
        setnames(c("chr", "start", "end", "quality", "strand")) %>%
        .[, chr := sub("Chr_", "chr", chr)] %>%
        .[grep("^chr", chr)]
    saveRDS(d[, .(quality, strand, chr, ID)], paste0("MQ_", ID, ".RDS"))

    d[quality >= MQ] %>%
        .[, coord := ifelse(strand == "+", start, end)] %>%
        .[, .(chr, start = coord, end = coord + 1, ID, strand)]
}
setkey(dStarts, chr, start, end, strand)


################################################################################
#
#                              MAPPING STATS
#
################################################################################


# Mapping statistics; 
# 1.1 MQ 
# 1.2 n of reads (per chr, strand) ...

########################################
# MQ
########################################

dMQ <- foreach(i = files, .combine = rbind) %do% {
    ID <- sub("_.*", "", sub(".*_0", "S", i))
    readRDS(paste0("MQ_", ID, ".RDS"))[, .(ID, quality)]
}

ggplot(dMQ, aes(quality, color = ID)) +
    geom_density() +
    labs(x = "Mapping Quality",
         title = "Mapping Quality (sample ID)") +
    theme_classic() +
    scale_fill_brewer(palette = "Dark2")

ggplot(dMQ, aes(quality, color = ID)) +
    labs(x = "Mapping Quality",
         title = "Mapping Quality (sample ID)") + 
    theme_classic() +
    geom_histogram()

ggplot(dMQ, aes(quality, fill = ID)) +
    labs(x = "Mapping Quality",
         title = "Mapping Quality (sample ID)") + 
    theme_classic() +
    geom_bar()

ggplot(dMQ, aes(ID, quality, fill = ID)) +
    geom_boxplot() + 
    labs(x = "Sample ID",
         y = "Mapping Quality",
         title = "Mapping Quality (sample ID)") +
    theme_classic() +
    scale_fill_brewer(palette = "Dark2")



dMQ <- foreach(i = files, .combine = rbind) %do% {
    ID <- sub("_.*", "", sub(".*_0", "S", i))
    readRDS(paste0("MQ_", ID, ".RDS"))[, .(ID, strand, quality, chr)]
}
ggplot(dMQ, aes(strand, quality, fill = ID)) +
    geom_boxplot() +
    labs(y = "Mapping Quality",
         x = "Strand", 
         title = "Mapping Quality (sample ID, strand)") +
    theme_classic() +
    scale_fill_brewer(palette = "Dark2")

ggplot(dMQ, aes(chr, quality, fill = ID)) +
    geom_boxplot() +
    labs(y = "Mapping Quality",
         x = "Chromosome", 
         title = "Mapping Quality (sample ID, chromosome)") +
    theme_classic() +
    scale_fill_brewer(palette = "Dark2")



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


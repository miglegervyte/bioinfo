################################################################################
#
#                                   META
#
################################################################################

########################################
# Arguments
########################################
oCG <- 0
# Mapping quality
MQ <- 30
# Path
pathCG <- "./CG.RDS"


########################################
# Sources
########################################

library(data.table)
library(dplyr)
library(ggplot2)
library(foreach)
library(fst)
library(ggpubr)


########################################
# Other
########################################

options(scipen = 999)


# Add normal names
dSamplenames <- data.frame(ID = c("Sample_001", "Sample_006", "Sample_007", "Sample_031", 
                                  "Sample_032", "Sample_082"),
                           new = c("K_0-vitC", "hmC_0-vitC_R1", "hmC_0vitC_R2", "hmC_5vitC_R1", 
                                   "hmC_5vitC_R2", "K_5-vitC"))



################################################################################
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


library("ggpubr")

files <- list.files(pattern = ".*mapped.*bed")
ends <- fread("ends.bed", col.names = c("ID", "length"))

    
    dCoverage <- foreach(i = files, .combine = cbind) %do% {
        ID <- sub(".fastq.bed", "", paste0("Sample_", strsplit(i, "_")[[1]][4])) 
        print(ID)
        dOrig <- fread(i) %>%
            setnames(c("chr", "start", "end", "ID", "mq", "strand"))
        foo <- dOrig %>%
            .[mq > 0] %>% 
            .[, chr := sub("Chr_", "chr", chr)] %>%
            .[grep("chr", chr)] 
        bar <- foo[, .N, ID][N == 1, ID]
        foo <- foo[ID %in% bar] %>%
            .[, .(chr, start = ifelse(strand == "+", start, end - 1), strand, ID)] %>%
            .[, end := start + 1] %>%
            .[, .(chr, start, end, strand, ID)] 
        starts <- foo[, c("chr", "start", "end", "strand", "ID")] %>%
            merge(ends, by = "ID")
        starts <- starts[!duplicated(starts[, c("start", "length")]),] %>%
            .[, c("chr", "start", "end", "strand", "ID")] %>%
            setDT() 
        writeBed(starts, "starts")
        dDistance <- fread("bedtools closest -a starts.bed -b CG.bed -d")
        dDistance[, .N, V5][, .N, N]    
        foo <- dDistance[, .N, V5][N == 1, V5]
        saveRDS(dDistance[V5 %in% foo, .N, .(strand = V4, distance = V9)], 
                paste0("stats_Distance_", ID, ".RDS"), compress = FALSE)
        foo <- dDistance[V9 <= 0, .N, .(chr = V6, start = V7)]
        setDT(merge(dCG, foo, c("chr", "start"), all.x = TRUE))[is.na(N), N := 0]$N
    }
    colnames(dCoverage) <- c("K_0-vitC", "hmC_0-vitC_R1", "hmC_0vitC_R2", "hmC_5vitC_R1", "hmC_5vitC_R2", "K_5-vitC")
    dCoverage <- as.data.frame(dCoverage)
    saveRDS(dCoverage, paste0("dCoverage.RDS"), compress = FALSE)
}

d <- readRDS("stats_Distance_Sample_001.RDS")
d[, distance2 := distance]
d[distance > 10, distance2 := "> 10"]

distance1 <- ggplot(d, aes(distance2, N, fill = strand)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "Distances K_0-vitC")


d <- readRDS("stats_Distance_Sample_006.RDS")
d[, distance2 := distance]
d[distance > 10, distance2 := "> 10"]

distance2 <- ggplot(d, aes(distance2, N, fill = strand)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "Distances hmC_0-vitC_R1")


d <- readRDS("stats_Distance_Sample_007.RDS")
d[, distance2 := distance]
d[distance > 10, distance2 := "> 10"]

distance3 <- ggplot(d, aes(distance2, N, fill = strand)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "Distances hmC_0-vitC_R2")


d <- readRDS("stats_Distance_Sample_031.RDS")
d[, distance2 := distance]
d[distance > 10, distance2 := "> 10"]

distance4 <- ggplot(d, aes(distance2, N, fill = strand)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "Distances hmC_5vitC_R1")


d <- readRDS("stats_Distance_Sample_032.RDS")
d[, distance2 := distance]
d[distance > 10, distance2 := "> 10"]

distance5 <- ggplot(d, aes(distance2, N, fill = strand)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "Distances hmC_5vitC_R2")


d <- readRDS("stats_Distance_Sample_082.RDS")
d[, distance2 := distance]
d[distance > 10, distance2 := "> 10"]

distance6 <- ggplot(d, aes(distance2, N, fill = strand)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "Distances K_5-vitC")


distances <- ggarrange(distance1, distance2, distance3, distance4, distance5, distance6, ncol = 2, nrow = 3)
ggsave("distances.pdf")

dCoverage <- readRDS("dCoverage.RDS")
dcor <- cor(dCoverage)
dmelt <- melt(dcor)
setDT(dmelt)
dmelt[Var1 == Var2, value := NA]
ggplot(dmelt, aes(x=Var1, y=Var2, fill=value)) +
    geom_tile()
ggsave("heatmap.pdf")

################################################################################
#
#                                   META
#
################################################################################

########################################
# Arguments
########################################

distanceToCG <- 0
# Mapping quality
MQ <- 30
# Path
<<<<<<< HEAD
pathCG <- "~/Desktop/bioinfo/Coprinopsis/naujas/bed_html/CG.RDS"
pathBED <- "~/Desktop/bioinfo/Coprinopsis/naujas/bed_html"
=======
pathCG <- "../genome/CG.RDS"
pathBED <- "./"
>>>>>>> 2c306c35163cc7a5129a6fe3f4b0875d73cfd562

########################################
# Sources
########################################

library(data.table)
<<<<<<< HEAD
=======
library(foreach)
library(ggplot2)
>>>>>>> 2c306c35163cc7a5129a6fe3f4b0875d73cfd562
library(dplyr)
library(ggplot2)
library(foreach)
library(fst)


################################################################################
#
#                                 FUNCTIONS
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

<<<<<<< HEAD
dCG <- readRDS(pathCG)
=======
dCG <- setDT(readRDS(pathCG))
>>>>>>> 2c306c35163cc7a5129a6fe3f4b0875d73cfd562
if (!file.exists("CG.bed")) {
    writeBed(dCG, "CG")
}


########################################
# Load BED files
########################################

<<<<<<< HEAD
setwd("~/Desktop/bioinfo/Coprinopsis/naujas/bed_html")


files <- list.files(pattern = ".*mapped.*bed")

dCoverage <- foreach(i = files[1:2], .combine = cbind) %do% {
    ID <- sub(".fastq.bed", "", paste0("Sample_", strsplit(i, "_")[[1]][4])) 
    print(ID)
    
    dOrig <- fread(i, nrows = 100) %>%
        setnames(c("chr", "start", "end", "ID", "mq", "strand")) %>%
        
        .[grep("Chr", chr)] 
    
    saveRDS(dOrig[, .N, .(chr, strand, mq)], 
            paste0("stats_MQ_", ID, ".RDS"), compress = FALSE)
    
    dOrig %>%
        .[mq > MQ] %>%
        .[, chr := sub("Chr_", "chr", chr)] %>%
=======
files <- list.files(pattern = ".*Ion.*bed")

dCoverage <- foreach(i = files[1:2], .combine = cbind) %do% {
    ID <- paste0("Sample_", strsplit(i, "_")[[1]][5])
    print(ID)
    dOrig <- fread(i) %>%
        setnames(c("chr", "start", "end", "ID", "mq", "strand"))
    # Perkelto chr filtravima cia!!!
    saveRDS(dOrig[, .N, .(chr, strand, mq)], 
            paste0("stats_MQ_", ID, ".RDS"), compress = FALSE)
    dOrig %>%
        .[mq > MQ] %>%
        .[, chr := sub("Chr_", "chr", chr)] %>%
        .[grep("chr", chr)] %>%
>>>>>>> 2c306c35163cc7a5129a6fe3f4b0875d73cfd562
        # PADARYTI A
        # .[, .(chr, start = ifelse(strand == "+", start, end), strand, ID)] %>%
        .[, .(chr, start = ifelse(strand == "+", start, end), strand)] %>%
        .[, end := start + 1] %>%
        .[, .(chr, start, end, strand, ID)] %>%
<<<<<<< HEAD
        .[, ID := ID] %>% 
        # replace(, V5=="Sample_006", "K_0-vitC") %>%
        writeBed("starts")
    
    system("bedtools sort -i CG.bed > CG2.bed")
    dDistance <- fread("bedtools closest -a starts.bed -b CG2.bed -d")
    
    saveRDS(dDistance[, .N, .(chr = V1, strand = V4, distance = V9)], 
            paste0("stats_Distance_", ID, ".RDS"), compress = FALSE)
    
=======
        .[, ID := ID] %>%
        writeBed("starts")
    dDistance <- fread("bedtools closest -a starts.bed -b CG.bed -d")
    saveRDS(dDistance[, .N, .(chr = V1, strand = V4, distance = V9)], 
        paste0("stats_Distance_", ID, ".RDS"), compress = FALSE)
>>>>>>> 2c306c35163cc7a5129a6fe3f4b0875d73cfd562
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

files <- list.files(pattern = ".*stats.*MQ.*RDS")




################################################################################
#
#                              MAPPING STATS
#
################################################################################


########################################
# MQ
########################################

files <- list.files(pattern = ".*stats.*MQ.*RDS")

dMQ <- foreach(i = files, .combine = rbind) %do% {
    ID <- sub(".RDS", "", paste0("Sample_", strsplit(i, "_")[[1]][4]))
<<<<<<< HEAD
    print(ID)
    readRDS(i)[, ID := ID]
}


mqSample <- ggplot(dMQ, aes(mq, N, color = ID)) +
    geom_smooth(se = FALSE) + 
    labs(x = "Mapping quality",
         y = "Number of reads",
         title = "Mapping Quality"
         subtitle = "Per sample")

mqStrand <- ggplot(dMQ, aes(mq, N, color = strand)) +
    theme_classic() +
    geom_histogram() +
    labs(x = "Mapping quality",
         y = "Number of reads",
         title = "Mapping Quality"
         subtitle = "Per strand")

mqChr <- ggplot(dMQ, aes(mq, N, color = chr)) +
    geom_histogram() +
    theme_classic() +
    labs(x = "Mapping quality",
         y = "Number of reads",
         title = "Mapping Quality"
         subtitle = "Per chromosome")

# labs(x = "Mapping Quality",
#      title = "Mapping Quality (sample ID)") +
# theme_classic() +
# scale_fill_brewer(palette = "Dark2")
=======
    readRDS(i)[, ID := ID]
}

ggplot(dMQ, aes(mq, N, color = ID)) +
    geom_smooth(se = FALSE)
    # labs(x = "Mapping Quality",
    #      title = "Mapping Quality (sample ID)") +
    # theme_classic() +
    # scale_fill_brewer(palette = "Dark2")
>>>>>>> 2c306c35163cc7a5129a6fe3f4b0875d73cfd562

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

files <- list.files(pattern = ".*stats.*Distance.*RDS")
dMQ <- foreach(i = files, .combine = rbind) %do% {
    ID <- sub(".RDS", "", paste0("Sample_", strsplit(i, "_")[[1]][4]))
    readRDS(i)[, ID := ID]
}


<<<<<<< HEAD
pd <- dMQ[, sum(N), ID, strand]
psample <- ggplot(pd, aes(ID, chr)) +
=======
pd <- dMQ[, sum(N), ID]
p <- ggplot(pd, aes(ID, V1)) +
>>>>>>> 2c306c35163cc7a5129a6fe3f4b0875d73cfd562
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "Sample ID",
         y = "Number of reads",
         title = "Number of reads",
         subtitle = "Per sample") +
    theme_classic() +
    scale_fill_brewer(palette = "Dark2")

<<<<<<< HEAD
#pr <- dMQ[, sum(N), strand]
#ggplot(dMQ, aes(ID, chr, fill = strand)) +
#    geom_bar(stat = "identity", position = "dodge") +
#    labs(x = "Sample ID",
#         y = "Number of reads",
#         title = "Number of reads",
#         subtitle = "Per sample") +
#    theme_classic() +
#    scale_fill_brewer(palette = "Dark2")


dMQ[, .N, .(ID, chr)] %>%
    pchr <- ggplot(aes(ID, N, group = chr)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "Chromosome",
         y = "Number of reads",
         title = "Number of reads"
         subtitle = "Per chromosome") +
=======
ggplot(pd, aes(ID, V1, fill = strand)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "Sample ID",
         y = "Number of reads",
         title = "Number of reads",
         subtitle = "Per sample") +
>>>>>>> 2c306c35163cc7a5129a6fe3f4b0875d73cfd562
    theme_classic() +
    scale_fill_brewer(palette = "Dark2")


<<<<<<< HEAD
=======
dMQ[, .N, .(ID, chr)] %>%
    ggplot(aes(ID, N, group = chr)) +
        geom_bar(stat = "identity", position = "dodge") +
        labs(x = "Chromosome",
             y = "Number of reads",
             title = "Number of reads (chromosome)") +
        theme_classic() +
        scale_fill_brewer(palette = "Dark2")


>>>>>>> 2c306c35163cc7a5129a6fe3f4b0875d73cfd562

########################################
# Distance to CG
########################################
<<<<<<< HEAD
dsample <-  ggplot(dMQ, aes(distance, color = ID)) +
=======
p <-  ggplot(dMQ, aes(distance, color = ID)) +
>>>>>>> 2c306c35163cc7a5129a6fe3f4b0875d73cfd562
    labs(x = "Distance",
         y = "Number of reads",
         title = "Distance to CG"
         subtitle = "Per sample") +
    theme_classic() +
    facet_wrap(~ strand) +
    geom_density()
<<<<<<< HEAD

dchr <-  ggplot(dMQ, aes(distance, color = chr)) +
    labs(x = "Distance",
         y = "Number of reads",
         title = "Distance to CG"
         subtitle = "Per chromosome") +
    theme_classic() +
    facet_wrap(~ strand) +
    geom_density()


dMQ[, distanceGroup := "0"]
dMQ[distance > 0 & distance <= 5, distanceGroup := "1-5"]
dMQ[distance > 5 , distanceGroup := "6-Inf"]

dgr <- dMQ[, sum(N), .(distanceGroup, ID)] %>%
    ggplot(aes(distanceGroup, V1, fill = ID)) +
    geom_bar(stat = "identity", position = "dodge",
             color = "black", width = 0.8) +
    labs(title = "Distance to CG",
         x = "Distance",
         y = "Number of reads",
         fill = "Sample") +
    scale_x_discrete(limits = c("0", "1-5", "6-Inf")) +
    scale_fill_brewer(palette = "Dark2") +
    theme_classic()
=======


dMQ[, distanceGroup := "0"]
dMQ[distance > 0 & distance <= 5, distanceGroup := "1-5"]
dMQ[distance > 5 , distanceGroup := "6-Inf"]

p1 <- dMQ[, sum(N), .(distanceGroup, ID)] %>%
    ggplot(aes(distanceGroup, V1, fill = ID)) +
        geom_bar(stat = "identity", position = "dodge",
                 color = "black", width = 0.8) +
        labs(title = "Distance to CG",
             x = "Distance",
             y = "Number of reads",
             fill = "Sample") +
        scale_x_discrete(limits = c("0", "1-5", "6-Inf")) +
        scale_fill_brewer(palette = "Dark2") +
        theme_classic()

foo <- dMQ[, sum(N), .(distanceGroup, ID)]
bar <- dMQ[, sum(N), ID]
p2 <- merge(foo, bar, c("ID"))[, per := V1.x * 100/ V1.y] %>%
    ggplot(aes(distanceGroup, per, fill = ID)) +
        geom_bar(stat = "identity", position = "dodge",
                 color = "black", width = 0.8) +
        labs(title = "Distance to CG",
             x = "Distance",
             y = "Number of reads",
             fill = "Sample") +
        scale_x_discrete(limits = c("0", "1-5", "6-Inf")) +
        scale_fill_brewer(palette = "Dark2") +
        scale_y_continuous(limits = c(0, 100)) +
        theme_classic()

ggsave("~/tmp.pdf", p2)





>>>>>>> 2c306c35163cc7a5129a6fe3f4b0875d73cfd562




foo <- dMQ[, sum(N), .(distanceGroup, ID)]
bar <- dMQ[, sum(N), ID]
dgr100 <- merge(foo, bar, c("ID"))[, per := V1.x * 100/ V1.y] %>%
    ggplot(aes(distanceGroup, per, fill = ID)) +
    geom_bar(stat = "identity", position = "dodge",
             color = "black", width = 0.8) +
    labs(title = "Distance to CG",
         x = "Distance",
         y = "Number of reads",
         fill = "Sample") +
    scale_x_discrete(limits = c("0", "1-5", "6-Inf")) +
    scale_fill_brewer(palette = "Dark2") +
    scale_y_continuous(limits = c(0, 100)) +
    theme_classic()

ggarrange(mqSample, mqStrand, mqChr, pd, pchr, dsample, dchr, dgr, dgr100, ncol = 3, nrow = 3)

cormat <- cor(dDistance)
melted_cormat <- melt(cormat)
ggplot(data = melted_cormat, aes (x = Var1, y = Var2, fill = value)) +
    geom_tile()

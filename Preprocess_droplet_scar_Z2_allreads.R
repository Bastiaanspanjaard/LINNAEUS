# Description ####
# Pre-process scar output files: merge with succesful whole-transcriptome
# barcodes, remove missequenced scars per cell, validate outcome by hand.

# Parameters ####
log2.cutoff <- 3

# Dependencies ####
source("./Scripts/scar_helper_functions.R")
# Getting scar indices
f <- function (i, j, n) {
  ifelse((i > j) & (j <= n), (j - 1) * (2 * n - 2 - j) / 2 + (i - 1), NA_real_)
}
get.dist.index <- function(k, n){
  # Get indices for the dth element of a distance object of distances between
  # n elements.
  # Starting position for each column
  ptr_all_cols <- f(2:n, 1:(n - 1), n)
  # Maximum valid `k`
  k_max <- n * (n - 1) / 2
  
  if (k > k_max) return(c(i = NA_real_, j = NA_real_))
  j <- sum(ptr_all_cols <= k)  ## get column index j
  i <- k - ptr_all_cols[j] + j + 1  ## get row index i
  return(c(i = i, j = j))
}


# Load data ####
scars.in <- read.csv("./Data/2017_10X_2/scar_2_filtered_scars.csv",
                     stringsAsFactors = F, sep = "\t")
wt.barcodes <- read.csv("./Data/2017_10X_2/Z2_1_bc_preprocessed_noshadow.csv",
                        stringsAsFactors = F, sep = ",")
scars.unfiltered <- merge(scars.in, wt.barcodes)
scars.unfiltered$Scar.id <- 1:nrow(scars.unfiltered)
scars.unfiltered$Keep <- T
scars.unfiltered$Pair <- "With"
cells <- unique(scars.unfiltered$Cell)

# Filter scars per cell ####
# Go through all cells one by one. Within each cell, calculate the Hamming
# distances between the sequences; if these are one or two and the read 
# difference is high (log2(read1/read2) > log2.cutoff), flag the scar for 
# removal.
for(c in 1:length(cells)){
  cell <- cells[c]
  cell.scars <- scars.unfiltered[scars.unfiltered$Cell == cell, ]
  if(nrow(cell.scars) < 2){next}
  cell.scars <- cell.scars[order(-cell.scars$Reads), ]
  n.scars <- nrow(cell.scars)
  cell.seqdist <- stringdistmatrix(cell.scars$Sequence, method = "hamming")
  
  # Identify all sequences that have to be removed.
  for(k in 1:length(cell.seqdist)){
    if(cell.seqdist[k] >2){next}
    else{
      # Get scar indices. Calculate log2 ratio between two reads. If too high,
      # flag lowest scar (Keep = F).
      scar.indices <- get.dist.index(k, n.scars)
      high.scar <- min(scar.indices)
      low.scar <- max(scar.indices)
      if(log2(cell.scars$Reads[high.scar]/cell.scars$Reads[low.scar]) > log2.cutoff){
        scars.unfiltered$Keep[scars.unfiltered$Scar.id == cell.scars$Scar.id[low.scar]] <-
          F
      }
    }
  }
}
# From assessments of overlapping scar, the following scar is very likely
# to be a sequencing error.
# scars.unfiltered$Keep[scars.unfiltered$Scar.id == 903] <- F

# Remove scars that have HD 1 to scars with much higher reads.
scars.filter.1 <- scars.unfiltered[scars.unfiltered$Keep, ]

# Go throught the barcodes and calculate the distances between the scars again. 
# If the HD is one and the read difference is low, mark the scars as suspect. 
# If the HD is two and the read difference is high, mark the scars as suspect.
for(c in 1:length(cells)){
  cell <- cells[c]
  cell.scars <- scars.filter.1[scars.filter.1$Cell == cell, ]
  if(nrow(cell.scars) < 2){next}
  cell.scars <- cell.scars[order(-cell.scars$Reads), ]
  n.scars <- nrow(cell.scars)
  cell.seqdist <- stringdistmatrix(cell.scars$Sequence, method = "hamming")

  for(k in 1:length(cell.seqdist)){
    if(cell.seqdist[k] >2){next}
    else if(cell.seqdist[k] == 1){
      # Get scar indices. Flag both (paste(Pair, other
      # scar)).
      scar.indices <- get.dist.index(k, n.scars)
      high.scar <- min(scar.indices)
      low.scar <- max(scar.indices)
      # If the scars are legitimate scars that are known to be close in HD,
      # skip this pair.
      if((cell.scars$Sequence[high.scar] %in% allowed.close.scars) &
         (cell.scars$Sequence[low.scar] %in% allowed.close.scars)){next}
      scars.filter.1$Pair[scars.filter.1$Scar.id == cell.scars$Scar.id[low.scar]] <-
        paste(scars.filter.1$Pair[scars.filter.1$Scar.id == cell.scars$Scar.id[low.scar]],
              scars.filter.1$Scar.id[scars.filter.1$Scar.id == cell.scars$Scar.id[high.scar]],
              sep = ".")
      scars.filter.1$Pair[scars.filter.1$Scar.id == cell.scars$Scar.id[high.scar]] <-
        paste(scars.filter.1$Pair[scars.filter.1$Scar.id == cell.scars$Scar.id[high.scar]],
              scars.filter.1$Scar.id[scars.filter.1$Scar.id == cell.scars$Scar.id[low.scar]],
              sep = ".")
      }
    else if(cell.seqdist[k] == 2){
      # Get scar indices. Calculate log2 ratio between two reads. If too high,
      # flag both (paste(Pair, other scar)).
      scar.indices <- get.dist.index(k, n.scars)
      high.scar <- min(scar.indices)
      low.scar <- max(scar.indices)
      if(log2(cell.scars$Reads[high.scar]/cell.scars$Reads[low.scar]) > log2.cutoff){
        scars.filter.1$Pair[scars.filter.1$Scar.id == cell.scars$Scar.id[low.scar]] <-
          paste(scars.filter.1$Pair[scars.filter.1$Scar.id == cell.scars$Scar.id[low.scar]],
                scars.filter.1$Scar.id[scars.filter.1$Scar.id == cell.scars$Scar.id[high.scar]],
                sep = ".ratio.")
        scars.filter.1$Pair[scars.filter.1$Scar.id == cell.scars$Scar.id[high.scar]] <-
          paste(scars.filter.1$Pair[scars.filter.1$Scar.id == cell.scars$Scar.id[high.scar]],
                scars.filter.1$Scar.id[scars.filter.1$Scar.id == cell.scars$Scar.id[low.scar]],
                sep = ".ratio.")
      }
    }
  }
}

# Assess and remove suspect scars ####
scars.assess <- scars.filter.1[scars.filter.1$Pair != "With", ]
scars.filter.2 <- scars.filter.1[, 1:7]
scars.assess$Incidence <-
  sapply(scars.assess$Sequence,
         function(x) sum(scars.filter.1$Sequence == x))
scars.assess$Pair.incidence <-
  sapply(scars.assess$Pair,
         function(x){
           pairs <- unlist(strsplit(x, "[.]"))[-1]
           return(min(scars.assess$Incidence[scars.assess$Scar.id %in% pairs]))}
  )
scars.assess$Min.pair.incidence <-
  apply(scars.assess[, c("Incidence", "Pair.incidence")], 1, min)

scars.assess.2 <- scars.assess[scars.assess$Min.pair.incidence > 1, ]

seq.freq <- data.frame(table(scars.assess.2$Sequence))
colnames(seq.freq)[1] <- "Sequence"
scars.assess.2 <- merge(scars.assess.2, seq.freq)
scars.assess.2$Crit.1 <- 
  !(scars.assess.2$Min.pair.incidence == scars.assess.2$Freq)

# Read in all >1 read sequences (all UMIs)
all.scars.g1 <- read.table("./Data/2017_10X_2/scar_2_reads_over1.txt",
                           sep = "\t", stringsAsFactors = F)
colnames(all.scars.g1)[-1] <- c("Barcode", "UMI", "Location", "Sequence")
all.scars.g1$V11 <- trimws(all.scars.g1$V1)

all.scars.g1$Reads <- sapply(all.scars.g1$V11,
                             function(x) unlist(strsplit(x, " "))[1])
all.scars.g1$CIGAR <- sapply(all.scars.g1$V11,
                             function(x) unlist(strsplit(x, " "))[2])
all.barcode.UMIs <-
  data.frame(table(all.scars.g1$Sequence, all.scars.g1$Barcode))
all.barcode.UMIs <- all.barcode.UMIs[all.barcode.UMIs$Freq > 0, ]
colnames(all.barcode.UMIs) <- c("Sequence", "Barcode", "UMIs")
scars.assess.2 <- merge(scars.assess.2, all.barcode.UMIs)
scars.assess.2$Crit.2 <- (scars.assess.2$UMIs > 1)

scars.assess.2.low.inc <- 
  scars.assess.2[scars.assess.2$Incidence <= scars.assess.2$Pair.incidence, ]

scars.assess.2.out <- 
  scars.assess.2.low.inc[!scars.assess.2.low.inc$Crit.1 & 
                           !scars.assess.2.low.inc$Crit.2, ]
scars.assess.2.maybe <- 
  scars.assess.2.low.inc[scars.assess.2.low.inc$Crit.1 | 
                           scars.assess.2.low.inc$Crit.2, ]

scars.assess.2.maybe$Crit.3 <- NA
# c.s <- 24
for(c.s in 1:nrow(scars.assess.2.maybe)){
  c.sequence <- scars.assess.2.maybe$Sequence[c.s]
  c.barcode <- scars.assess.2.maybe$Barcode[c.s]
  c.scar <- scars.assess.2.maybe$Scar.id[c.s]
  c.parent.sequence <- 
    scars.assess$Sequence[grepl(paste(c.scar, "\\.", sep = ""), scars.assess$Pair) |
                              grepl(paste(c.scar, "$", sep = ""), scars.assess$Pair)]
  
  validate.seq <- 
    all.scars.g1[all.scars.g1$Barcode == c.barcode &
                   all.scars.g1$Sequence %in% c.sequence, ]
  validate.parent <- 
    all.scars.g1[all.scars.g1$Barcode == c.barcode &
                   all.scars.g1$Sequence %in% c.parent.sequence, ]
  View(validate.seq)
  View(validate.parent)
  crit.3 <- 
    readline(prompt = paste("Scar ", c.s, " of ", nrow(scars.assess.2.maybe), 
                            ": amount of UMIs with HD > 1 to the parent sequence UMIs? ", sep = ""))
  scars.assess.2.maybe$Crit.3[c.s] <- crit.3
}
scars.assess.2.maybe$Out <-ifelse(scars.assess.2.maybe$Crit.3 ==0 , T, F)


sequencing.error.scars <-
  c(scars.assess.2.out$Scar.id, 
    scars.assess.2.maybe$Scar.id[scars.assess.2.maybe$Out])

scars.output <- scars.filter.1[!(scars.filter.1$Scar.id %in% sequencing.error.scars),
                               1:6]

# Write filtered results ####
# write.csv(scars.output, "./Data/2017_10X_2/Z2_1_preprocessed_scars_allreads.csv",
#           quote = F, row.names = F)

# Check scars per UMI ####
cell.scars <- data.frame(table(scars.filter.1$Barcode))
colnames(cell.scars)[1] <- "Barcode"
cell.scars$Mean.UMI.seq <-
  sapply(cell.scars$Barcode,
         function(x) {
           mean(table(all.scars.g1$UMI[all.scars.g1$Barcode == x]))})

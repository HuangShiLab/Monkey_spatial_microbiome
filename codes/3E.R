#######################################
# Monkey Study
#
# Mantel test, mocosal and content
#
# Author: HOU Shuwen
#######################################
library(FEAST)

# Content -> Mucosa
# Set working directory and load your own data
setwd("~/Downloads/monkey_data")
content_genus <- read.table("2b_species_abundance.txt", header = TRUE, row.names = 1)

setwd("~/Downloads/monkey_mucosa")
mucosa_genus <- read.table("species_abundance.txt", header = TRUE, row.names = 1)

# Transpose if needed so that samples are rows
content_genus <- as.matrix(content_genus)
mucosa_genus <- as.matrix(mucosa_genus)

# Find sample names common to both content and mucosa
common_samples <- intersect(colnames(content_genus), colnames(mucosa_genus))
content_genus <- content_genus[,common_samples , drop = FALSE]
mucosa_genus <- mucosa_genus[,common_samples , drop = FALSE]

# Merge source and sink
merged <- merge(content_genus, mucosa_genus, by = 0, all = T)
merged[is.na(merged)] <- 0
rownames(merged)<-merged$Row.names
merged<-merged[,-1]

# Normalize
merged <- ceiling((merged / colSums(merged))*10000)
merged <-t(merged)


# Grab all column names except the first one ("Row.names") that merge() created
sample_cols <- setdiff(colnames(merged), "Row.names")

# Construct metadata
metadata <- data.frame(
  Env    = ifelse(grepl("\\.x$", sample_cols), "content", "mucosa"),
  SourceSink  = ifelse(grepl("\\.x$", sample_cols), "source",  "sink"),
  id          = sub("\\.[xy]$", "", sample_cols),      # strip suffix to get the shared ID
  row.names   = sample_cols,
  stringsAsFactors = FALSE
)

# Estimate source proportions for each sink
FEAST_output <- FEAST(C = merged, metadata = metadata, different_sources_flag = 1,
                      dir_path ="~/Downloads/monkey_data", outfile="C_M")
output <- read.table("C_M_transmission.txt", header = TRUE, row.names = 1)

# Violin plot
plot <- ggplot(output, aes(x = "", y = transmission)) +
  geom_violin(trim = FALSE, fill = "#b996c5") +
  geom_boxplot(width = 0.1, fill = "black") +
  ylim(0, 1) + labs(x = NULL, y = "Transmission proportion", title="Content to Mucosa") +
  theme(panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 12),
        axis.line = element_line(size = 0.5, color = "black"))
ggsave("plots/3E.png", plot, width = 3, height = 5)


# Mucosa -> Blood
# Set working directory and load your own data
setwd("~/Downloads/monkey_mucosa")
mucosa_genus <- read.table("species_abundance.txt", header = TRUE, row.names = 1)

# Transpose if needed so that samples are rows
mucosa_genus <- as.matrix(mucosa_genus)

# Blood
blood_samples <- grep("B_", colnames(mucosa_genus), value = TRUE)

# Find sample names common to both content and mucosa
blood_genus <- mucosa_genus[,blood_samples , drop = FALSE]
mucosa_genus <- mucosa_genus[,common_samples , drop = FALSE]

# Merge source and sink
merged <- merge(mucosa_genus, blood_genus, by = 0, all = T)
merged[is.na(merged)] <- 0
rownames(merged)<-merged$Row.names
merged<-merged[,-1]

# Normalize
merged <- ceiling((merged / colSums(merged))*10000)
merged <- merged[rowSums(merged) > 0, , drop = FALSE]
merged <-t(merged)


# Construct metadata
metadata <- data.frame(
  Env        = c(rep("mucosa", length(common_samples)),
                 rep("blood",  length(blood_samples))),
  SourceSink = c(rep("Source", length(common_samples)),
                 rep("Sink",   length(blood_samples))),
  id         = c(rep(1, 54),
                 1:20),              
  row.names  = c(common_samples, blood_samples),
  stringsAsFactors = FALSE
)

metadata$id <- sub("_.*", "",        # keep text before first "_"
                       rownames(metadata))
metadata$id <- match(metadata$id, c("O1", "O2", "Y1", "Y2"))


# Estimate source proportions for each sink
FEAST_output <- FEAST(C = merged, metadata = metadata, different_sources_flag = 0,
                      dir_path ="~/Downloads/monkey_data", outfile="M_B")

# Draw a sankey bubble plot

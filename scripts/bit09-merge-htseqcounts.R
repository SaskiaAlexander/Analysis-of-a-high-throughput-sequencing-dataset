################################################################################
### Merge the count files after running htseq-count using this R script
### Download files from server on your local machine using FileZilla
################################################################################
### Get names count files
################################################################################
pathFiles <- "C:/Users/alexa/OneDrive/Howest/BIT 01 Linux operating systems/SF/BIT09/data/07_hisat2_trimmed_genome_tran_htseqcount"
patternFiles <- ".txt"
patternAfterSampleName <- "_filtered_sorted"
patternBeforeSampleName <- "counts_"
countFiles <- list.files(path = pathFiles,
                         pattern = patternFiles,
                         all.files = TRUE,
                         recursive = FALSE,
                         ignore.case = FALSE,
                         include.dirs = FALSE)
################################################################################
### Get counts from first file/sample in data frame
################################################################################
pathCountFile1 = paste(pathFiles, 
                      countFiles[1], 
                      sep = "/")
df.count <- read.table(pathCountFile1,
                       header = FALSE, 
                       stringsAsFactors = FALSE)
################################################################################
### Set column names
################################################################################
# Capture text using () and reprint using \\1
column2 <- gsub("(.*).txt", "\\1", countFiles[1])
column2 <- gsub(patternAfterSampleName, "", column2)
column2 <- gsub(patternBeforeSampleName, "", column2)
colnames(df.count) <- c("ID", column2)
################################################################################
### Add counts from other samples in data frame
################################################################################
for (i in 2:length(countFiles) ) {
  pathCountFile = paste(pathFiles, 
                        countFiles[i], 
                        sep = "/")
  df.next <- read.table(pathCountFile,
                        header = FALSE, 
                        stringsAsFactors = FALSE)
  column2 <- gsub("(.*).txt", "\\1", countFiles[i])
  column2 <- gsub(patternAfterSampleName, "", column2)
  column2 <- gsub(patternBeforeSampleName, "", column2)
  colnames(df.next) <- c("ID", column2)
  df.count <- merge(df.count, df.next, by = c("ID"))
}
# Show first 10 rows, first 5 columns
df.count[1:10,1:4]
################################################################################
### Calculate total amount counts per sample, get summary and write to file
################################################################################
### IDs as rownames
rownames(df.count) <- df.count$ID
df.count <- df.count[,-1]
### Calculate total amount counts per sample (added as last line)
df.count <- rbind(df.count, total.counts = colSums(df.count))
# tail(df.count[,c(1:3)])
### Get __lines
data.summary <- df.count[grep("__", rownames(df.count), 
                              perl = TRUE, 
                              invert = FALSE), ]
### Get last line with total.counts
lastLine <- df.count[grep("total.counts", rownames(df.count), 
                          perl = TRUE, 
                          invert = FALSE), ]
data.summary <- rbind(data.summary, lastLine)
### transpose (pivottable) -> t()
write.table(t(data.summary), 
            file = paste0(pathFiles,"/summary_counts.txt"), 
            sep = "\t",
            col.names = NA)
################################################################################
### Create data frame with all count values and save to file
### Remove __rows and last total.counts line
################################################################################
### invert = true -> everything but __ lines (remove __lines)
allCounts <- df.count[grep("__", rownames(df.count), 
                           perl = TRUE, 
                           invert = TRUE), ]
### every row but the last one (total counts row)
write.table(allCounts[-nrow(allCounts),], 
            file = paste0(pathFiles, "/all_counts.txt"), 
            sep = "\t",
            col.names = NA)
################################################################################
################################################################################
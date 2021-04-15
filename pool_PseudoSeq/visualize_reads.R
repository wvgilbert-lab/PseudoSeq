library(ggplot2)
library(tidyr)
library(dplyr)


#filenames = c("fname1","fname2",..."fnamen")

# Replace filenames here! Filenames is a list of input files,
# figures will be save to outfile, which should be a .png
filenames <- c("../../PUS1_+CMC.txt","../../PUS1_-CMC.txt",
               "../../TRUB1_+CMC.txt","../../TRUB1_-CMC.txt")
outfile <- "all_reads.pdf"

# Load the reads from all files into a list of dataframes,

data <- list()    
classes = c("character","numeric","numeric")

for (i in seq_along(filenames)) {
    data[[i]] <- read.csv(file = filenames[i],
                          colClasses=classes,
                          sep="\t",header=FALSE) %>%
        mutate(fname=filenames[i])
}

# Combine data into a single data frame, then update column names
data <- bind_rows(data)
colnames(data)<- c("Sequence","Position","Reads","fname")




# Count the number of reads mapped to each Sequence
totals <- data %>%
    group_by(Sequence,fname) %>%
    summarize(Total=sum(Reads)) 

# Use the sequence totals to normalize the reads at each position
data <- merge(data,totals,
              by=c("Sequence","fname")) %>%
    filter(Total > 100) %>%
    mutate(ReadFrac=Reads/Total) 


# Plot the read fractions.
ggplot(data) +
    facet_wrap(~fname,
               ncol=1,
               scales="fixed") +
    geom_point(aes(x=Position,
                    y=ReadFrac),
                size=0.5,
                alpha=0.5) +
    ylab("fraction of reads at position") +
    xlab("position") +
#    scale_y_log10()
    theme_classic() 

# Save the figure to outfile, scaling the figure size to
# the number of panels/input files
N.files <- length(filenames)
ggsave(outfile,
       height=N.files*1.3,width=4.5)

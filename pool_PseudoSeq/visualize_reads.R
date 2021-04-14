library(tidyverse)

classes=c("character","numeric","numeric")

data <- rbind(read.csv("../PUS1_+CMC.txt",
                       colClasses=classes,
                       sep="\t",header=FALSE) %>%
                 mutate(Treatment="CMC"),
              read.csv("../PUS1_-CMC.txt",
                       colClasses=classes,
                       sep="\t",header=FALSE) %>%
                 mutate(Treatment="mock"))

colnames(data)<- c("Sequence","Position","Reads","Treatment")

totals <- data %>%
    group_by(Sequence,Treatment) %>%
    summarize(Total=sum(Reads)) #%>%
    #filter(Treatment=="CMC") %>%
    #arrange(desc(Total))


sequences = c("ENST00000361624_U391","ENST00000361381_U396")

# Print a "browser view" for a given site
ggplot(data %>% filter(Sequence=="ENST00000361624_U391")) +
    facet_wrap(~Treatment,
               scales='fixed',
               ncol=1) +
    geom_vline(xintercept=69,
               color="#E69F00",
               linetype="dashed") +
    geom_bar(aes(x=Position,
                 y=Reads,
                 fill=Treatment),
             stat="identity") +
    scale_fill_manual(values=c("#000000", "#666666")) +
    xlab("position") +
    ylab("read 5' ends") +
    theme_classic() +
    theme(legend.position="none",
          strip.background = element_blank(),
          strip.text.x = element_blank())
ggsave("sample_site.pdf",
       height=3,width=4.5)


data <- merge(data,totals,
              by=c("Sequence","Treatment")) %>%
    filter(Total > 100) %>%
    mutate(ReadFrac=Reads/Total)


ggplot(data) +
    facet_wrap(~Treatment,
               ncol=1,
               scales="fixed") +
    geom_vline(xintercept=35,color="blue") +
    geom_vline(xintercept=55,color="green") +
    geom_vline(xintercept=66,color="red") +
    geom_vline(xintercept=95,color="green") +
    geom_vline(xintercept=115,color="blue") +
    geom_point(aes(x=Position,
                   y=ReadFrac),
               alpha=0.5) +
    scale_y_log10() +
    theme_bw()
ggsave("../plots/all_reads_log.png")

ggplot(data) +
    facet_wrap(~Treatment,
               ncol=1,
               scales="fixed") +
     geom_point(aes(x=Position,
                    y=ReadFrac),
                size=0.5,
                alpha=0.5) +
    ylab("fraction of reads at position") +
    xlab("position") +
    theme_classic() +
    theme(legend.position="none",
          strip.background = element_blank(),
          strip.text.x = element_blank())

 #   geom_vline(xintercept=35,color="blue") +
  #  geom_vline(xintercept=60,color="green") +
#    geom_vline(xintercept=66,color="red") +
 #   geom_vline(xintercept=95,color="green") +
  #  geom_vline(xintercept=115,color="blue") +
   
ggsave("../plots/all_reads.png",
       height=3,width=4.5)

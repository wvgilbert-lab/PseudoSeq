library(tidyverse)

# Load data from all the different analyses and conditions and
# merge them into a single dataframe

data <- rbind(read.csv("../lr_PUS1_+CMC.peaks.txt",sep="\t") %>%
                 mutate(Treatment="CMC"),
              read.csv("../lr_PUS1_-CMC.peaks.txt",sep="\t") %>%
                 mutate(Treatment="mock")) %>%
    filter(reads>10) %>%
    mutate(Analysis="lr_wide")

#data2 <- rbind(read.csv("../alt_PUS1_+CMC.peaks.txt",sep="\t") %>%
 #                mutate(Treatment="CMC"),
  #            read.csv("../alt_PUS1_-CMC.peaks.txt",sep="\t") %>%
   #              mutate(Treatment="mock")) %>%
    #filter(reads>10) %>%
    #mutate(Analysis="cov_only")

#data <- rbind(data,data2)

data2 <- rbind(read.csv("../PUS1_+CMC.peaks.txt",sep="\t") %>%
                 mutate(Treatment="CMC"),
              read.csv("../PUS1_-CMC.peaks.txt",sep="\t") %>%
                 mutate(Treatment="mock")) %>%
    filter(reads>10) %>%
    mutate(Analysis="original")

data <- rbind(data,data2)

data2 <- rbind(read.csv("../lrnarrow_PUS1_+CMC.peaks.txt",sep="\t") %>%
                 mutate(Treatment="CMC"),
              read.csv("../lrnarrow_PUS1_-CMC.peaks.txt",sep="\t") %>%
                 mutate(Treatment="mock")) %>%
    filter(reads>10) %>%
    mutate(Analysis="lr_narrow")

data <- rbind(data,data2)



# Visualize peak heights relative to library position
# comparing the +CMC condition across analyses
ggplot(data %>% filter(Treatment=="CMC")) +
    facet_wrap(~Analysis,
               ncol=1) +
    geom_point(aes(x=pos,
                   y=Z,
                   color=reads)) +
    scale_x_continuous(limits=c(0,120)) +
    scale_y_continuous(limits=c(0,100)) + 
    theme_bw()
ggsave("peaks_per_pos.png")

# Visualize the peak height distribution as a CDF
ggplot(data %>% filter(pos>=60 & pos<=95)) +
    stat_ecdf(aes(x=Z,
                  color=Treatment,
                  linetype=Analysis)) +
    coord_cartesian(xlim = c(-3, 15)) +
    theme_bw()
ggsave("peakdist.png",height=4,width=4.3)

# Visualize the peak heigh distribution as a log-scaled inverse CDF
ggplot(data) +
    # plot 1 - CDF
    stat_ecdf(aes(x=Z,
                  y=1-..y..,
                  color=Treatment,
                  linetype=Analysis)) +
    # Manually plot the cutoff
    geom_vline(xintercept=11,
               color='blue') +
    # log scale, change the background
    scale_y_log10() +
    theme_bw()
ggsave('cutoff.pdf')


# Compare median peak heights relative to library position across
# conditions and treatments

summary <- data %>%
    group_by(Treatment,Analysis,pos) %>%
    summarize(p100=max(Z),
              p95=quantile(Z,.95),
              p75=quantile(Z,.75),
              p50=median(Z),
              p25=quantile(Z,.25),
              p05=quantile(Z,.05))

ggplot(summary,
       aes(x=pos)) +
    facet_wrap(~Treatment,
               ncol=1) +
    geom_ribbon(aes(ymin=p25,ymax=p75,
                    fill=Analysis),
                alpha=0.4) +
    geom_line(aes(y=p50,
                  color=Analysis)) +
    geom_vline(xintercept=70) +
    theme_bw()
ggsave("peakdist_per_pos.pdf",
       height=4,width=6)


ggplot(summary,
       aes(x=pos)) +
    facet_wrap(~Treatment,
               ncol=1) +
    geom_ribbon(aes(ymin=p25,ymax=p75,
                    fill=Analysis),
                alpha=0.4) +
    geom_line(aes(y=p50,
                  color=Analysis)) +
    geom_vline(xintercept=70) +
    coord_cartesian(ylim=c(-2,5)) +
    theme_bw()
ggsave("peakdist_per_pos.filtered.pdf",
       height=4,width=6)

data <- data %>%
    select(chrom,pos,strand,Analysis,Z,Treatment) %>%
    spread(Treatment,Z) %>%
    arrange(desc(CMC)) %>%
    replace_na(list(CMC=0.0,
                    mock=0.0))



ggplot(data) +
    facet_wrap(~Analysis,
               nrow=1) +
    geom_point(aes(x=mock,
                   y=CMC)) +
    coord_cartesian(xlim = c(-3, 60),
                    ylim = c(-3, 60)) +
    theme_bw()
ggsave("scatter.png",height=3.5,width=7)


ggplot(data %>% filter(Analysis=="lr_narrow")) +
    geom_point(aes(x=mock,
                   y=CMC),
               alpha=0.7) +
    coord_cartesian(xlim = c(-3, 60),
                    ylim = c(-3, 60)) +
    xlab("site Z-score in -CMC library") +
    ylab("site Z-score in +CMC library") +
    
    theme_classic()
ggsave("scatter_panel.pdf",height=4,width=4)

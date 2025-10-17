---
jupyter:
  kernelspec:
    display_name: R
    language: R
    name: ir
  language_info:
    codemirror_mode: r
    file_extension: .r
    mimetype: text/x-r-source
    name: R
    pygments_lexer: r
    version: 4.3.3
  nbformat: 4
  nbformat_minor: 5
---

::: {#ae649819-e78d-476a-85d3-d2833e309def .cell .code execution_count="2"}
``` R
# install.packages("pals")
# install.packages("ggpubr")
# install.packages("igraph")
# install.packages("ggrepel")
```
:::

::: {#448aa600-0bb4-4660-a7fe-8f6c74e505f2 .cell .code execution_count="5"}
``` R
# LOAD LIBRARIES ----------------------------------------------------------

library(ggpubr)
library(tidyverse)
#library(igraph)
library(ggrepel)
```

::: {.output .stream .stderr}
    Loading required package: ggplot2

    ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ✔ lubridate 1.9.3     ✔ tibble    3.2.1
    ✔ purrr     1.0.2     ✔ tidyr     1.3.1
    ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ✖ dplyr::filter() masks stats::filter()
    ✖ dplyr::lag()    masks stats::lag()
    ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
:::
:::

::: {#0b8a2416-6083-4ff7-b36c-189f9a7bdc9b .cell .code}
``` R
# INITIALIZATION ----------------------------------------------------------

# Create a variable to our working directory for produced outputs
WORKDIR <- getwd()
FILESOURCE <- file.path(WORKDIR, "data")

set.seed(16341)

# This is a nice color palette.
cols <- pals::tableau20(20)
# Take a look at it.
scales::show_col(cols)

popCols <- c("#1F77B4", "#FFBFD4", "#82CBFF", "#2CA02C", "#930000", "#000000")
names(popCols) <- c('African','American','EastAsian','European','SouthAsian','Middle East')
ref_pop <- c("European", "East Asian", "American", "South Asian", "African")   
```
:::

::: {#8233dd2e-6ab9-40d9-9d59-5254222438bd .cell .code}
``` R
# FUNCTIONS ---------------------------------------------------------------

## UAE onto G1K projection plot
panelPCA.main <- function(plot.df, metric.df, p.axes, anc.cutoff = 0.9, l.pos = "none") {
  
  col.var <- "Main.ancestry"
  fill.var <- "Main.ancestry"
  
  ## 1. G1K as background
  g1k.df <- dplyr::filter(plot.df, grouping == "1KGenome")
  ## 2. all the ERGP samples
  ergp.df <- dplyr::filter(plot.df, grouping == "Emirati")
  ## 3. all Emiratis with main ancestry exceeding 90%
  ## 2. all the ERGP samples
  uae_prime.df <- dplyr::filter(ergp.df, European > eval(anc.cutoff) | 
                                  SouthAsian >= eval(anc.cutoff) | 
                                  EastAsian >= eval(anc.cutoff) | 
                                  African >= eval(anc.cutoff) | 
                                  American >= eval(anc.cutoff))
  
  
  tmp.return <- ggplot(data=plot.df, aes(x=get(p.axes[1]), y=get(p.axes[2]), fill=get(fill.var), colour=get(col.var))) +
    geom_point(data=g1k.df, size=5, alpha = 0.65) + 
    geom_point(data=ergp.df, shape=21, size=4, colour="black", alpha = 0.85) +
    geom_point(data=uae_prime.df, shape=23, size=7.5, colour="white", alpha = 0.85) + 
    scale_shape_manual(values=c(21, 22)) +
    scale_fill_manual(values = popCols, name = "Main Ancestry",
                      labels = c("European","South Asian","African","American","East Asian"), 
                      breaks = c("European","SouthAsian","African","American","EastAsian")) +
    scale_colour_manual(values = popCols, name = "Main Ancestry",
                        labels = c("European","South Asian","African","American","East Asian"), 
                        breaks = c("European","SouthAsian","African","American","EastAsian")) +
    xlim(metric.df[eval(p.axes[1]), "axis.min"], metric.df[eval(p.axes[1]), "axis.max"]) +
    ylim(metric.df[eval(p.axes[2]), "axis.min"], metric.df[eval(p.axes[2]), "axis.max"]) +
    xlab(metric.df[eval(p.axes[1]), "var.label"]) +
    ylab(metric.df[eval(p.axes[2]), "var.label"]) + 
    theme_bw() +
    theme(legend.position=eval(l.pos),
          panel.border = element_blank(),            # Remove panel border
          panel.background = element_blank(),        # Remove panel background
          plot.background = element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x = element_text(size =20),
          axis.title.y = element_text(size =20),
          legend.title = element_text(size =20),
          legend.text = element_text(size =20))
  
  return(tmp.return)
  
}

## G1k onto UAE projection
panelPCA.main_UAE <- function(plot.df, metric.df, p.axes, anc.cutoff = 0.9, l.pos = "none") {
  
  col.var <- "Main.ancestry"
  fill.var <- "Main.ancestry"
  
  ## 1. G1K as background
  g1k.df <- dplyr::filter(plot.df, grouping == "1KGenome")
  ## 2. all the ERGP samples
  ergp.df <- dplyr::filter(plot.df, grouping == "Emirati")
  ## 3. all Emiratis with main ancestry exceeding 90%
  ## 2. all the ERGP samples
  uae_prime.df <- dplyr::filter(ergp.df, European > eval(anc.cutoff) | 
                                  SouthAsian >= eval(anc.cutoff) | 
                                  EastAsian >= eval(anc.cutoff) | 
                                  African >= eval(anc.cutoff) | 
                                  American >= eval(anc.cutoff))
  
  
  tmp.return <- ggplot(data=plot.df, aes(x=get(p.axes[1]), y=get(p.axes[2]), fill=get(fill.var), colour=get(col.var))) +
    geom_point(data=ergp.df, shape=21, size=2, alpha = 0.95) +
    geom_point(data=g1k.df, shape=21, colour="black", size=5, alpha = 0.85) + 
    geom_point(data=uae_prime.df, shape=23, size=7.5, colour="white", alpha = 0.85) + 
    scale_shape_manual(values=c(21, 22)) +
    scale_fill_manual(values = popCols, name = "Main Ancestry",
                      labels = c("European","South Asian","African","American","East Asian"), 
                      breaks = c("European","SouthAsian","African","American","EastAsian")) +
    scale_colour_manual(values = popCols, name = "Main Ancestry",
                        labels = c("European","South Asian","African","American","East Asian"), 
                        breaks = c("European","SouthAsian","African","American","EastAsian")) +
    xlim(metric.df[eval(p.axes[1]), "axis.min"], metric.df[eval(p.axes[1]), "axis.max"]) +
    ylim(metric.df[eval(p.axes[2]), "axis.min"], metric.df[eval(p.axes[2]), "axis.max"]) +
    xlab(metric.df[eval(p.axes[1]), "var.label"]) +
    ylab(metric.df[eval(p.axes[2]), "var.label"]) + 
    theme_bw() +
    theme(legend.position=eval(l.pos), 
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x = element_text(size =20),
          axis.title.y = element_text(size =20),
          #axis.line=element_blank(),
          panel.border = element_blank(),            # Remove panel border
          panel.background = element_blank(),        # Remove panel background
          plot.background = element_blank(),
          legend.title = element_text(size =20),
          legend.text = element_text(size =20))
  
  return(tmp.return)
  
}
```
:::

::: {#6ad4804a-dcc0-4daf-a2ef-2b53961af4a7 .cell .code}
``` R
# Admixture Bar-Plot ILLUMINA ---------------------------------------------

## load dataframe
ancestry_order <- c("European", "SouthAsian", "American", "African", "EastAsian")

admix.df <- readRDS(file.path("../data/admix_dummy_data.RDS"))


PG.df <- admix.df  %>%
  mutate(Main.ancestry = sapply(main_pop, function(elemelon) str_replace_all(elemelon, "\\s+", ""))) %>%
  rename(
    SouthAsian = "South Asian",
    EastAsian = "East Asian"
  )


df_long <- PG.df %>%
  mutate(Main.ancestry = factor(Main.ancestry, levels = ancestry_order)) %>% # Order by Main.ancestry
  arrange(
    Main.ancestry, # Order by Main.ancestry groups
    desc(case_when( # Order within each group by the respective ancestry column
      Main.ancestry == "European" ~ European,
      Main.ancestry == "SouthAsian" ~ SouthAsian,
      Main.ancestry == "African" ~ African,
      Main.ancestry == "American" ~ American,
      Main.ancestry == "EastAsian" ~ EastAsian
    ))
  ) %>%
  dplyr::select(sampleID, Main.ancestry, European, SouthAsian, American, African, EastAsian) %>%
  mutate(Prefix = factor(sampleID, levels = sampleID)) %>%
  pivot_longer(cols = European:EastAsian, names_to = "Ancestry", values_to = "value") %>%
  mutate(Ancestry = factor(Ancestry, levels = ancestry_order))


## Sanity ##
popCols <- c("#1F77B4", "#FFBFD4", "#82CBFF", "#2CA02C", "#930000", "#000000")
names(popCols) <- c('African','American','EastAsian','European','SouthAsian','Middle East')
## Sanity ##
```
:::

::: {#aeca318c-e0ed-40d0-a74e-462577a7cf7e .cell .code}
``` R
######## AD plot EAS only
df_long.sa <- df_long %>%
  dplyr::filter(Main.ancestry == "EastAsian")

ad_plot.sa <- ggplot(df_long.sa, aes(x=Prefix, y=value, fill=Ancestry)) + 
  geom_bar(stat = "identity", width=0.9) + 
  scale_fill_manual(values = popCols, name = "Admixture",
                    labels = c("European","South Asian","African","American","East Asian"), 
                    breaks = c("European","SouthAsian","African","American","EastAsian")) +
  scale_y_continuous(expand = c(0,0))+ ylab("Admixture proportion") +
  theme_bw() +
  facet_grid(. ~ Main.ancestry, scales = "free", space = "free", 
             labeller = labeller(Main.ancestry = c(
               "European" = "EUR",
               "SouthAsian" = "SAS",
               "African" = "AFR",
               "American" = "AM",
               "EastAsian" = "EAS"
             ))) +
  theme(
    legend.position = "none", # Top left corner inside the plot
    legend.justification = c(0, 1), # Align the legend to the top left
    legend.direction = "vertical", # Vertical legend
    legend.box.background = element_rect(color = "black", fill = "white", size = 1), # Black frame with white background
    legend.title = element_text(size = 20), 
    legend.text = element_text(size = 16), 
    #axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
    panel.background = element_blank(),        # Remove panel background
    plot.background = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    strip.text = element_text(size = 20),
    strip.background = element_rect(fill = NA, color = NA), 
    panel.border = element_blank()
  ) 
```
:::

::: {#5e056f62-ce2a-4b57-9eec-62e0510aefde .cell .code}
``` R
### no EAS
df_long.NOsa <- df_long %>%
  dplyr::filter(Main.ancestry != "EastAsian")

ad_plot.NOsa <- ggplot(df_long.NOsa, aes(x=Prefix, y=value, fill=Ancestry)) + 
  geom_bar(stat = "identity", width=0.9) + 
  scale_fill_manual(values = popCols, name = "Admixture",
                    labels = c("European","South Asian","African","American","East Asian"), 
                    breaks = c("European","SouthAsian","African","American","EastAsian")) +
  scale_y_continuous(expand = c(0,0))+ ylab("Admixture proportion") +
  theme_bw() +
  facet_grid(. ~ Main.ancestry, scales = "free", space = "free", 
             labeller = labeller(Main.ancestry = c(
               "European" = "EUR",
               "SouthAsian" = "SAS",
               "African" = "AFR",
               "American" = "AM",
               "EastAsian" = "EAS"
             ))) +
  theme(
    legend.position = c(0.05, 0.95), # Top left corner inside the plot
    legend.justification = c(0, 1), # Align the legend to the top left
    legend.direction = "vertical", # Vertical legend
    legend.box.background = element_rect(color = "black", fill = "white", size = 1), # Black frame with white background
    legend.title = element_text(size = 20), 
    legend.text = element_text(size = 16), 
    #axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
    panel.background = element_blank(),        # Remove panel background
    plot.background = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 20),
    strip.text = element_text(size = 20),
    strip.background = element_rect(fill = NA, color = NA), 
    panel.border = element_blank()
  ) 
```
:::

::: {#e66e792f-79f4-411c-a40a-f2b108b3bf99 .cell .code}
``` R
# Admixture PCA-Plot ILLUMINA ---------------------------------------------

df <- readRDS(file = file.path("../data/PCA_dummy.RDS"))
plot.df <- df[[1]]
metric.df <- df[[2]]


plot.df <- plot.df  %>%
  mutate(Main.ancestry = sapply(main_pop, function(elemelon) str_replace_all(elemelon, "\\s+", ""))) %>%
  mutate(supPopNameClean = sapply(supPopNameClean, function(elemelon) str_replace_all(elemelon, "\\s+", ""))) %>%
  rename(
    SouthAsian = "South Asian",
    EastAsian = "East Asian"
  )


pca12.main <- panelPCA.main(plot.df, metric.df, c("PC1","PC2"), anc.cutoff = 0.75)
pca23.main <- panelPCA.main(plot.df, metric.df, c("PC2","PC3"), anc.cutoff = 0.75)
pca34.main <- panelPCA.main(plot.df, metric.df, c("PC3","PC4"), anc.cutoff = 0.75)

#

# Admixture UAE PCA-Plot ILLUMINA -----------------------------------------

df <- readRDS(file = file.path("../data/PCA_dummy.RDS"))

plot.df <- df[[1]]  %>%
  mutate(Main.ancestry = sapply(main_pop, function(elemelon) str_replace_all(elemelon, "\\s+", ""))) %>%
  mutate(supPopNameClean = sapply(supPopNameClean, function(elemelon) str_replace_all(elemelon, "\\s+", ""))) %>%
  rename(
    SouthAsian = "South Asian",
    EastAsian = "East Asian"
  )
metric.df <- df[[2]]


pca12.uae <- panelPCA.main_UAE(plot.df, metric.df, c("PC1","PC2"), anc.cutoff = 0.75)
pca23.uae <- panelPCA.main_UAE(plot.df, metric.df, c("PC2","PC3"), anc.cutoff = 0.75)
pca34.uae <- panelPCA.main_UAE(plot.df, metric.df, c("PC3","PC4"), anc.cutoff = 0.75)


#
```
:::

::: {#ff131dcf-6f95-408a-84d2-c9565ee5bb0b .cell .code}
``` R
# MAIN ANCESTRY PREVALENCE CHART ------------------------------------------

ancestry_order <- c("European", "SouthAsian", "American", "African", "EastAsian")

admix.df <- readRDS(file.path("../data/admix_dummy_data.RDS"))


PG.df <- admix.df  %>%
  mutate(Main.ancestry = sapply(main_pop, function(elemelon) str_replace_all(elemelon, "\\s+", ""))) %>%
  rename(
    SouthAsian = "South Asian",
    EastAsian = "East Asian"
  ) %>%
  dplyr::mutate(European = European * 100,
                SouthAsian = SouthAsian * 100,
                African = African * 100,
                American = American * 100,
                EastAsian = EastAsian * 100)

df_cnt <- PG.df %>%
  dplyr::mutate(Main.ancestry = gsub(" ", "", Main.ancestry)) %>%
  group_by(Main.ancestry) %>%
  summarise(
    TotalSamples = n()
  ) %>%
  arrange(desc(TotalSamples)) %>%
  dplyr::mutate(Main.ancestry = factor(Main.ancestry, levels = Main.ancestry)) %>%
  mutate(prop = TotalSamples / sum(TotalSamples) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop ) %>%
  dplyr::mutate(group = "Group1", .before = Main.ancestry)

### make count frame for cutoff main ancestry values
anc.cutoff <- 90
uae_prime.df <- dplyr::filter(PG.df, European >= eval(anc.cutoff) | 
                                SouthAsian >= eval(anc.cutoff) | 
                                EastAsian >= eval(anc.cutoff) | 
                                African >= eval(anc.cutoff) | 
                                American >= eval(anc.cutoff)) %>%
  dplyr::mutate(Main.ancestry = gsub(" ", "", Main.ancestry)) %>%
  group_by(Main.ancestry) %>%
  summarise(
    TotalSamples = n()
  ) %>%
  arrange(desc(TotalSamples)) %>%
  dplyr::mutate(Main.ancestry = factor(Main.ancestry, levels = Main.ancestry)) %>%
  mutate(prop = TotalSamples / sum(TotalSamples) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop ) %>%
  dplyr::mutate(group = "Group2", .before = Main.ancestry)


### Put humpdy-Dumpdy back together again
sample_data <- rbind.data.frame(df_cnt, uae_prime.df) %>%
  arrange(Main.ancestry) %>%
  group_by(group) %>%
  mutate(total_prop = prop / sum(prop)) %>%
  mutate(nested_group = ifelse(group == "Group1", "Outer", "Inner"))


# Nested pie chart
ancestry_plot <- ggplot(sample_data, aes(x = nested_group, y = total_prop, fill = Main.ancestry)) +
  # Set alpha conditionally for outer and inner rings
  geom_bar(
    aes(alpha = ifelse(nested_group == "Outer", 0.75, 1)), # Outer ring brighter (alpha = 1), Inner ring dimmer (alpha = 0.7)
    stat = "identity", width = 0.95, linewidth = 0.25, color = "white"
  ) +
  scale_fill_manual(
    values = popCols,
    name = "Admixture",
    labels = c("European", "South Asian", "American", "African", "East Asian"), 
    breaks = c("European", "SouthAsian", "American", "African", "EastAsian")
  ) +
  scale_alpha_identity() + # Use identity scale to apply custom alpha directly
  coord_polar(theta = "y", start = 0) +
  theme_void() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),            # Remove panel border
    panel.background = element_blank(),        # Remove panel background
    plot.background = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  # Add repel labels with dynamic font color
  geom_label_repel(
    aes(label = TotalSamples, color = Main.ancestry), # Map color to Main.ancestry
    fontface = "bold", # Bold font for labels
    fill = "white", # Background color of the label
    box.padding = 0.5, # Padding around the label
    segment.color = "grey50", # Line connecting label to the segment
    segment.size = 0.75, # Thickness of the connecting line
    label.size = 0.25, # Border thickness around the label
    size = 5.5,
    show.legend = FALSE # Disable legend for label colors
  ) +
  # Ensure the text color matches the fill colors
  scale_color_manual(values = popCols, guide = "none")
```
:::

::: {#742afe38-198b-4cf9-9f7f-f608727c7528 .cell .code}
``` R
# MAIN ANCESTRY METRICS ---------------------------------------------------

admix.df <- readRDS(file.path("../data/admix_dummy_data.RDS"))

data <- admix.df %>%
  mutate(Main.ancestry = sapply(main_pop, function(elemelon) str_replace_all(elemelon, "\\s+", ""))) %>%
  rename(
    SouthAsian = "South Asian",
    EastAsian = "East Asian"
  ) %>%
  dplyr::select(sampleID, Main.ancestry, European, SouthAsian, American, African, EastAsian)


long_data <- data %>%
  pivot_longer(cols = European:EastAsian, names_to = "Component", values_to = "Proportion") %>%
  mutate(
    Subset = ifelse(Component == Main.ancestry, paste(eval(Component), "Main", sep = "_"), "Overall")
  )

### subset component distribution for main ancestry
long_data.sub <- long_data %>%
  dplyr::filter(Subset != "Overall") %>%
  dplyr::mutate(Component = Subset)

plotframe <- rbind.data.frame(long_data, long_data.sub) %>%
  dplyr::mutate(col.alpha = ifelse(grepl("Main", Component), 0.75, 1)) %>%
  dplyr::mutate(Component = factor(Component, levels = c("European", "European_Main", "SouthAsian", "SouthAsian_Main",
                                                         "American", "American_Main", "African", "African_Main",
                                                         "EastAsian", "EastAsian_Main")))
## SANITY CHECK
table(plotframe$Component, plotframe$col.alpha)

## adjust color palette
unique(plotframe$Component)

ancCols <- c("#2CA02C","#2CA02C","#930000","#930000","#FFBFD4","#FFBFD4","#1F77B4","#1F77B4", "#82CBFF", "#82CBFF")
names(ancCols) <- c("European", "European_Main", "SouthAsian", "SouthAsian_Main",
                    "American", "American_Main", "African", "African_Main",
                    "EastAsian", "EastAsian_Main")

plotframe <- plotframe %>%
  dplyr::mutate(facets = gsub("_Main", "", Component)) %>%
  dplyr::mutate(facets = factor(facets, levels = c("European", "SouthAsian", "American" , "African", "EastAsian")))

# Modify the plot
violin.plot <- ggplot(plotframe, aes(x = Component, y = Proportion, fill = Component)) +
  # Violin plot with dynamic alpha for "main" labels
  geom_violin(aes(alpha = ifelse(grepl("Main", Component, ignore.case = TRUE), 0.75, 1)),
              color = "black", scale = "width", trim = TRUE) +
  # Mean point
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "#FF7F0E") +
  # Median line
  stat_summary(fun = median, geom = "point", shape = 95, size = 17, color = "#7F7F7F") +
  # Scatter points for outliers exceeding the 95th quantile
  # geom_jitter(
  #   data = plotframe %>%
  #     group_by(Component) %>%
  #     mutate(Q95 = quantile(Proportion, 0.99)) %>%
  #     filter(Proportion > Q95),
  #   aes(x = Component, y = Proportion),
  #   color = "lightgrey", alpha = 0.15, size = 1, width = 0.2, inherit.aes = FALSE
  # ) +
  # Update fill scale and x-axis labels
  scale_fill_manual(
    values = ancCols,
    name = "Admixture",
    labels = c("European", "European Main", "South Asian", "South Asian Main",
               "American", "American Main", "African", "African Main",
               "East Asian", "East Asian Main"), 
    breaks = c("European", "European_main", "SouthAsian", "SouthAsian_main",
               "American", "American_main", "African", "African_main",
               "EastAsian", "EastAsian_main")
  ) +
  scale_alpha_identity() +  # Ensure alpha values are applied directly
  theme_minimal() +
  facet_grid(. ~ facets, scales = "free", space = "free", 
             labeller = labeller(facets = c(
               "European" = "EUR", 
               "SouthAsian" = "SAS", 
               "American" = "AM", 
               "African" = "AFR", 
               "EastAsian" = "EAS"
             ))) +
  theme(
    legend.position = "none",
    legend.justification = c(0, 1),
    legend.direction = "vertical",
    legend.box.background = element_rect(color = "black", fill = "white", size = 1),
    legend.title = element_text(size = 20), 
    legend.text = element_text(size = 16), 
    panel.background = element_blank(),       
    plot.background = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 16),
    #axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 20),
    strip.text = element_text(size = 18),
    strip.background = element_rect(fill = NA, color = NA), 
    panel.border = element_blank(),
    panel.grid.minor = element_line(linewidth = 0.25),
    panel.grid.major = element_line(color = "gray80", linewidth = 0.5)
  ) 






#
```
:::

::: {#daf6667b-565f-4629-bab2-564cf18c32f1 .cell .code}
``` R
# MT HAPLOGROUP PLOT ------------------------------------------------------

df <- readRDS(file.path("../data/MTHaplo_dummy.RDS"))

le.filter <- data.frame(table(df$Macro_group)) %>%
  dplyr::filter(Freq >= 100)

ancestry.order.MT <- c("European", "South Asian", "American", "African", "East Asian")

df.plt <- df %>%
  dplyr::filter(Macro_group %in% le.filter$Var1) %>% 
  dplyr::mutate(Main.ancestry = factor(main_pop, levels = ancestry.order.MT))


popCols.mt <- c("#1F77B4", "#FFBFD4", "#82CBFF", "#2CA02C", "#930000", "#000000")
names(popCols.mt) <- c('African','American','East Asian','European','South Asian','Middle East')
## Sanity ##
table(df.plt$Main.ancestry)


mt.plot <- ggplot(df.plt, aes(x=Macro_group, y=Main.ancestry, fill=Main.ancestry)) + 
  geom_bar(stat = "identity", width=0.9) + 
  scale_fill_manual(values = popCols.mt, name = "Admixture",
                    labels = c("European","South Asian","African","American","East Asian"), 
                    breaks = c("European","South Asian","African","American","East Asian")) +
  theme_bw() +
  labs(
    title = "Violin Plots of Admixture Components",
    x = "Admixture Component",
    y = "Ancestry Proportion",
    fill = "Admixture Component"
  ) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1, size = 14),
    legend.position = "none",
    title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 20),
    strip.text = element_text(size = 20),
    strip.background = element_rect(fill = "white", color = NA), 
    panel.border = element_blank(),            
    panel.background = element_blank(),        
    plot.background = element_blank()
  ) 


#### ratio
df.plt_long <- df.plt %>%
  group_by(Macro_group, Main.ancestry) %>%
  summarise(Count = n(), .groups = "drop") %>%  # Count occurrences
  mutate(Log_Count = log10(Count)) %>%
  group_by(Macro_group) %>%
  mutate(ratio = Count/sum(Count))


mt.ratio.plot <- ggplot(df.plt_long, aes(x = Macro_group, y = ratio, fill = Main.ancestry)) + 
  # Stacked bar plot
  geom_bar(stat = "identity", width = 0.9) + 
  # Custom fill colors
  scale_fill_manual(
    values = popCols.mt,
    name = "Admixture",
    labels = c("European", "South Asian", "African", "American", "East Asian"), 
    breaks = c("European", "South Asian", "African", "American", "East Asian")
  ) +
  # Optional log10-transformation for y-axis
  #scale_y_continuous(trans = "log10", name = "Log10(Count)") + 
  # Theme adjustments
  theme_bw() +
  labs(
    title = "Stacked Bar Plot of Main Ancestry Components",
    x = "Macro Group",
    y = "Count",
    fill = "Main Ancestry"
  ) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1, size = 14),
    axis.text.y = element_text(size = 12),
    legend.position = "bottom",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    title = element_text(size = 16),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text = element_text(size = 14),
    strip.background = element_rect(fill = "white", color = NA), 
    panel.border = element_blank(),            
    panel.background = element_blank(),        
    plot.background = element_blank()
  )
```
:::

::: {#f7972166-cf8a-409d-89c7-142fc63a171c .cell .code}
``` R
# PLOT PANEL --------------------------------------------------------------


# Adjust margins for individual plots
ad_plot.NOsa <- ad_plot.NOsa + theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
ad_plot.sa <- ad_plot.sa + theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.01), "cm"))


Fig3.plt <- ggarrange(
  # First row with 90/10 width ratio
  ggarrange(
    ad_plot.NOsa, ad_plot.sa, 
    ncol = 2, nrow = 1, 
    widths = c(95, 5), 
    labels = c("A)"), 
    font.label = list(size = 24, color = "black", face = "bold")
  ), 
  # Second row
  ggarrange(
    pca12.main, pca12.uae, ancestry_plot, 
    ncol = 3, nrow = 1, 
    labels = c("B)", "C)", "D)"), 
    font.label = list(size = 24, color = "black", face = "bold")
  ), 
  # Third row
  ggarrange(
    mt.plot, violin.plot, 
    ncol = 2, nrow = 1, 
    widths = c(1, 2), 
    labels = c("E)", "F)"), 
    font.label = list(size = 24, color = "black", face = "bold")
  ), 
  # Combine rows
  ncol = 1, nrow = 3, 
  #  labels = c("A)", "B)", "C)"), 
  font.label = list(size = 24, color = "black", face = "bold",
                    panel.spacing = unit(0.1, "cm") )
)


ggsave(file.path("../results", "Figure_3_Panel.png"), plot = Fig3.plt, width = 13, height = 12, dpi = "retina")
```
:::

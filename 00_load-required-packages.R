###############################################################################
## ## Project: Collapsibility 
## Program 0: load required packages 
# Goal:
#### install and load required packages 

## Created by: Ya-Hui Yu (Jan 29, 2025)
################################################################################

pacman::p_load(
  rio,
  here,
  skimr,
  tidyverse,
  lmtest,
  sandwich,
  broom,
  parallel, 
  SuperLearner,
  tmle,
  #ggpubr,
  patchwork
)

thm <- theme_classic() +
  theme(
    legend.position = "top",
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.key = element_rect(fill = "transparent", colour = NA)
  )
theme_set(thm)

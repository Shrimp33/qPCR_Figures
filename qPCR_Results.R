library("readxl")

setwd("F:")
# Get all files that end with .xlsx
# files <- list.files(pattern = ".xlsx")

# get all sheets in "2023-12-05_ALI_mColon__huCYP_huALDH_huCYP24A1.xlsx"
sheet_names <- excel_sheets("2023-12-05_ALI_mColon__huCYP_huALDH_huCYP24A1.xlsx")
# Keep names that contain "REP"
sheet_names <- sheet_names[grep("REP", sheet_names)]
sheet_names

# # > sheet_names
# [1] "m_Msi1 REP"    "m_Lgr5 REP"    "m_Bmi1 REP"    "m_Cyp24a1 REP"
# [5] "m_Bcl2 REP"    "m_Ki67 REP"    "m_Caspn3 REP"

# Groups for each sample in rows
# NO DSS
# DSS
# DSS CYP
# DSS 1981
# DSS ALDH
# DSS MAC


m_MsiREP_2_neg_del_Ct <- read_excel("2023-12-05_ALI_mColon__huCYP_huALDH_huCYP24A1.xlsx", sheet = "m_Msi1 REP", range = "AI10:AQ16")
m_Lgr5REP_2_neg_del_Ct <- read_excel("2023-12-05_ALI_mColon__huCYP_huALDH_huCYP24A1.xlsx", sheet = "m_Lgr5 REP", range = "AI10:AQ16")
m_Bmi1REP_2_neg_del_Ct <- read_excel("2023-12-05_ALI_mColon__huCYP_huALDH_huCYP24A1.xlsx", sheet = "m_Bmi1 REP", range = "AI10:AQ16")
m_Cyp24a1REP_2_neg_del_Ct <- read_excel("2023-12-05_ALI_mColon__huCYP_huALDH_huCYP24A1.xlsx", sheet = "m_Cyp24a1 REP", range = "AI10:AQ16")
m_Bcl2REP_2_neg_del_Ct <- read_excel("2023-12-05_ALI_mColon__huCYP_huALDH_huCYP24A1.xlsx", sheet = "m_Bcl2 REP", range = "AI10:AQ16")
m_Ki67REP_2_neg_del_Ct <- read_excel("2023-12-05_ALI_mColon__huCYP_huALDH_huCYP24A1.xlsx", sheet = "m_Ki67 REP", range = "AI10:AQ16")
m_Caspn3REP_2_neg_del_Ct <- read_excel("2023-12-05_ALI_mColon__huCYP_huALDH_huCYP24A1.xlsx", sheet = "m_Caspn3 REP", range = "AI10:AQ16")

# For each of neg_del_Ct, make them dfs then add rownames
m_MsiREP_2_neg_del_Ct <- as.data.frame(m_MsiREP_2_neg_del_Ct)
m_Lgr5REP_2_neg_del_Ct <- as.data.frame(m_Lgr5REP_2_neg_del_Ct)
m_Bmi1REP_2_neg_del_Ct <- as.data.frame(m_Bmi1REP_2_neg_del_Ct)
m_Cyp24a1REP_2_neg_del_Ct <- as.data.frame(m_Cyp24a1REP_2_neg_del_Ct)
m_Bcl2REP_2_neg_del_Ct <- as.data.frame(m_Bcl2REP_2_neg_del_Ct)
m_Ki67REP_2_neg_del_Ct <- as.data.frame(m_Ki67REP_2_neg_del_Ct)
m_Caspn3REP_2_neg_del_Ct <- as.data.frame(m_Caspn3REP_2_neg_del_Ct)

rownames(m_MsiREP_2_neg_del_Ct) <- c("NO DSS", "DSS", "DSS CYP", "DSS 1981", "DSS ALDH", "DSS MAC")
rownames(m_Lgr5REP_2_neg_del_Ct) <- c("NO DSS", "DSS", "DSS CYP", "DSS 1981", "DSS ALDH", "DSS MAC")
rownames(m_Bmi1REP_2_neg_del_Ct) <- c("NO DSS", "DSS", "DSS CYP", "DSS 1981", "DSS ALDH", "DSS MAC")
rownames(m_Cyp24a1REP_2_neg_del_Ct) <- c("NO DSS", "DSS", "DSS CYP", "DSS 1981", "DSS ALDH", "DSS MAC")
rownames(m_Bcl2REP_2_neg_del_Ct) <- c("NO DSS", "DSS", "DSS CYP", "DSS 1981", "DSS ALDH", "DSS MAC")
rownames(m_Ki67REP_2_neg_del_Ct) <- c("NO DSS", "DSS", "DSS CYP", "DSS 1981", "DSS ALDH", "DSS MAC")
rownames(m_Caspn3REP_2_neg_del_Ct) <- c("NO DSS", "DSS", "DSS CYP", "DSS 1981", "DSS ALDH", "DSS MAC")


m_MsiREP_GAPDH <- read_excel("2023-12-05_ALI_mColon__huCYP_huALDH_huCYP24A1.xlsx", sheet = "m_Msi1 REP", range = "B10:K16")
m_Lgr5REP_GAPDH <- read_excel("2023-12-05_ALI_mColon__huCYP_huALDH_huCYP24A1.xlsx", sheet = "m_Lgr5 REP", range = "B10:K16")
m_Bmi1REP_GAPDH <- read_excel("2023-12-05_ALI_mColon__huCYP_huALDH_huCYP24A1.xlsx", sheet = "m_Bmi1 REP", range = "B10:K16")
m_Cyp24a1REP_GAPDH <- read_excel("2023-12-05_ALI_mColon__huCYP_huALDH_huCYP24A1.xlsx", sheet = "m_Cyp24a1 REP", range = "B10:K16")
m_Bcl2REP_GAPDH <- read_excel("2023-12-05_ALI_mColon__huCYP_huALDH_huCYP24A1.xlsx", sheet = "m_Bcl2 REP", range = "B10:K16")
m_Ki67REP_GAPDH <- read_excel("2023-12-05_ALI_mColon__huCYP_huALDH_huCYP24A1.xlsx", sheet = "m_Ki67 REP", range = "B10:K16")
m_Caspn3REP_GAPDH <- read_excel("2023-12-05_ALI_mColon__huCYP_huALDH_huCYP24A1.xlsx", sheet = "m_Caspn3 REP", range = "B10:K16")

m_Msi1REP_Ct <- read_excel("2023-12-05_ALI_mColon__huCYP_huALDH_huCYP24A1.xlsx", sheet = "m_Msi1 REP", range = "N10:W16")
m_Lgr5REP_Ct <- read_excel("2023-12-05_ALI_mColon__huCYP_huALDH_huCYP24A1.xlsx", sheet = "m_Lgr5 REP", range = "N10:W16")
m_Bmi1REP_Ct <- read_excel("2023-12-05_ALI_mColon__huCYP_huALDH_huCYP24A1.xlsx", sheet = "m_Bmi1 REP", range = "N10:W16")
m_Cyp24a1REP_Ct <- read_excel("2023-12-05_ALI_mColon__huCYP_huALDH_huCYP24A1.xlsx", sheet = "m_Cyp24a1 REP", range = "N10:W16")
m_Bcl2REP_Ct <- read_excel("2023-12-05_ALI_mColon__huCYP_huALDH_huCYP24A1.xlsx", sheet = "m_Bcl2 REP", range = "N10:W16")
m_Ki67REP_Ct <- read_excel("2023-12-05_ALI_mColon__huCYP_huALDH_huCYP24A1.xlsx", sheet = "m_Ki67 REP", range = "N10:W16")
m_Caspn3REP_Ct <- read_excel("2023-12-05_ALI_mColon__huCYP_huALDH_huCYP24A1.xlsx", sheet = "m_Caspn3 REP", range = "N10:W16")


library(tidyr)
library(dplyr)
library(tibble)

Msi1_exprs <- m_MsiREP_2_neg_del_Ct |>
tibble::rownames_to_column("Group") |> 
  tidyr::pivot_longer(cols=!Group, names_to = "Samples", values_to = "neg_del_Ct")

Lgr5_exprs <- m_Lgr5REP_2_neg_del_Ct |>
tibble::rownames_to_column("Group") |> 
  tidyr::pivot_longer(cols=!Group, names_to = "Samples", values_to = "neg_del_Ct")

Bmi1_exprs <- m_Bmi1REP_2_neg_del_Ct |>
tibble::rownames_to_column("Group") |> 
  tidyr::pivot_longer(cols=!Group, names_to = "Samples", values_to = "neg_del_Ct")

Cyp24a1_exprs <- m_Cyp24a1REP_2_neg_del_Ct |>
tibble::rownames_to_column("Group") |> 
  tidyr::pivot_longer(cols=!Group, names_to = "Samples", values_to = "neg_del_Ct")

Bcl2_exprs <- m_Bcl2REP_2_neg_del_Ct |>
tibble::rownames_to_column("Group") |> 
  tidyr::pivot_longer(cols=!Group, names_to = "Samples", values_to = "neg_del_Ct")

Ki67_exprs <- m_Ki67REP_2_neg_del_Ct |>
tibble::rownames_to_column("Group") |> 
  tidyr::pivot_longer(cols=!Group, names_to = "Samples", values_to = "neg_del_Ct")

Caspn3_exprs <- m_Caspn3REP_2_neg_del_Ct |>
tibble::rownames_to_column("Group") |> 
  tidyr::pivot_longer(cols=!Group, names_to = "Samples", values_to = "neg_del_Ct")

# Remove values equal to NA
Msi1_exprs <- Msi1_exprs[!is.na(Msi1_exprs$neg_del_Ct),]
Lgr5_exprs <- Lgr5_exprs[!is.na(Lgr5_exprs$neg_del_Ct),]
Bmi1_exprs <- Bmi1_exprs[!is.na(Bmi1_exprs$neg_del_Ct),]
Cyp24a1_exprs <- Cyp24a1_exprs[!is.na(Cyp24a1_exprs$neg_del_Ct),]
Bcl2_exprs <- Bcl2_exprs[!is.na(Bcl2_exprs$neg_del_Ct),]
Ki67_exprs <- Ki67_exprs[!is.na(Ki67_exprs$neg_del_Ct),]
Caspn3_exprs <- Caspn3_exprs[!is.na(Caspn3_exprs$neg_del_Ct),]

library(ggpubr)
library(ggplot2)

# # This is the comparison for Msi1

# compare_means(neg_del_Ct ~ Group,  data = Msi1_exprs, method="t.test")
# # Compare DSS to all and DSS MAC to all
# my_comparisons <- list(c("DSS", "DSS CYP"), c("DSS", "DSS 1981"), c("DSS", "DSS ALDH"), c("DSS", "DSS MAC"), c("DSS MAC", "DSS CYP"), c("DSS MAC", "DSS 1981"), c("DSS MAC", "DSS ALDH"), c("DSS", "NO DSS"))
# ggboxplot(Msi1_exprs, x = "Group", y = "neg_del_Ct",
#           color = "Group", palette = "jco")+ 
#   stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value     # Add global p-value

# ggsave("Msi1_relative_Expression.png")

# # Now generate graph for Lgr5
# compare_means(neg_del_Ct ~ Group,  data = Lgr5_exprs, method="t.test")
# # Compare DSS to all and DSS MAC to all
# my_comparisons <- list(c("DSS", "DSS CYP"), c("DSS", "DSS 1981"), c("DSS", "DSS ALDH"), c("DSS", "DSS MAC"), c("DSS MAC", "DSS CYP"), c("DSS MAC", "DSS 1981"), c("DSS MAC", "DSS ALDH"), c("DSS", "NO DSS"))
# ggboxplot(Lgr5_exprs, x = "Group", y = "neg_del_Ct",
#           color = "Group", palette = "jco")+ 
#   stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value     # Add global p-value

# ggsave("Lgr5_relative_Expression.png")

# # Now generate graph for Bmi1
# compare_means(neg_del_Ct ~ Group,  data = Bmi1_exprs, method="t.test")
# # Compare DSS to all and DSS MAC to all
# my_comparisons <- list(c("DSS", "DSS CYP"), c("DSS", "DSS 1981"), c("DSS", "DSS ALDH"), c("DSS", "DSS MAC"), c("DSS MAC", "DSS CYP"), c("DSS MAC", "DSS 1981"), c("DSS MAC", "DSS ALDH"), c("DSS", "NO DSS"))
# ggboxplot(Bmi1_exprs, x = "Group", y = "neg_del_Ct",
#           color = "Group", palette = "jco")+ 
#   stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value     # Add global p-value

# ggsave("Bmi1_relative_Expression.png")

# # Now generate graph for Cyp24a1
# compare_means(neg_del_Ct ~ Group,  data = Cyp24a1_exprs, method="t.test")
# # Compare DSS to all and DSS MAC to all
# my_comparisons <- list(c("DSS", "DSS CYP"), c("DSS", "DSS 1981"), c("DSS", "DSS ALDH"), c("DSS", "DSS MAC"), c("DSS MAC", "DSS CYP"), c("DSS MAC", "DSS 1981"), c("DSS MAC", "DSS ALDH"), c("DSS", "NO DSS"))
# ggboxplot(Cyp24a1_exprs, x = "Group", y = "neg_del_Ct",
#           color = "Group", palette = "jco")+ 
#   stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value     # Add global p-value

# ggsave("Cyp24a1_relative_Expression.png")

# # Now generate graph for Bcl2
# compare_means(neg_del_Ct ~ Group,  data = Bcl2_exprs, method="t.test")
# # Compare DSS to all and DSS MAC to all
# my_comparisons <- list(c("DSS", "DSS CYP"), c("DSS", "DSS 1981"), c("DSS", "DSS ALDH"), c("DSS", "DSS MAC"), c("DSS MAC", "DSS CYP"), c("DSS MAC", "DSS 1981"), c("DSS MAC", "DSS ALDH"), c("DSS", "NO DSS"))
# ggboxplot(Bcl2_exprs, x = "Group", y = "neg_del_Ct",
#           color = "Group", palette = "jco")+ 
#   stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value     # Add global p-value

# ggsave("Bcl2_relative_Expression.png")

# # Now generate graph for Ki67
# compare_means(neg_del_Ct ~ Group,  data = Ki67_exprs, method="t.test")
# # Compare DSS to all and DSS MAC to all
# my_comparisons <- list(c("DSS", "DSS CYP"), c("DSS", "DSS 1981"), c("DSS", "DSS ALDH"), c("DSS", "DSS MAC"), c("DSS MAC", "DSS CYP"), c("DSS MAC", "DSS 1981"), c("DSS MAC", "DSS ALDH"), c("DSS", "NO DSS"))
# ggboxplot(Ki67_exprs, x = "Group", y = "neg_del_Ct",
#           color = "Group", palette = "jco")+ 
#   stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value     # Add global p-value

# ggsave("Ki67_relative_Expression.png")

# # Now generate graph for Caspn3
# compare_means(neg_del_Ct ~ Group,  data = Caspn3_exprs, method="t.test")
# # Compare DSS to all and DSS MAC to all
# my_comparisons <- list(c("DSS", "DSS CYP"), c("DSS", "DSS 1981"), c("DSS", "DSS ALDH"), c("DSS", "DSS MAC"), c("DSS MAC", "DSS CYP"), c("DSS MAC", "DSS 1981"), c("DSS MAC", "DSS ALDH"), c("DSS", "NO DSS"))
# ggboxplot(Caspn3_exprs, x = "Group", y = "neg_del_Ct",
#           color = "Group", palette = "jco")+ 
#   stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value     # Add global p-value

# ggsave("Caspn3_relative_Expression.png")

listOfGEI <- list("Msi1"=Msi1_exprs, "Lgr5"=Lgr5_exprs, "Bmi1"=Bmi1_exprs, "Cyp24a1"=Cyp24a1_exprs, "Bcl2"=Bcl2_exprs, "Ki67"=Ki67_exprs, "Casp3"=Caspn3_exprs)

library("stringr")

backPaste <- function(x, y) {
  paste0(y, x)
}

mapGroup <- function(v, map_ref){
  # v is a vector containing values we want to map
  # map_ref is a matrix with the first column as the key and the second column as the value we want to map to
  l <- v
  for (i in seq_along(l)){
    l[i] <- map_ref[match(l[i], map_ref[,1]), 2]
  }
  return(l)
}

teqRep <-  matrix(c(c(1, 2, 3, 4, 5, 6, 7, 8, 9), c(rep(1, 3), rep(2, 3), rep(3, 3))), ncol = 2)

# For each element in listOfGEI add color column defined by sapply(as.integer((as.integer(sapply(neg_del_Ct, str_sub, -1, -1)) - 1)/ 3), backPaste, "_color")
for (i in seq_along(listOfGEI)) {
  listOfGEI[[i]]$color <- sapply(mapGroup(as.integer(sapply(listOfGEI[[i]]$Samples, str_sub, -1, -1)), teqRep), backPaste, "color_")
}


# sapply(as.integer((as.integer(sapply(Msi1_exprs$neg_del_Ct, str_sub, -1, -1)) - 1)/ 3), backPaste, "color_")

for (i in seq_along(listOfGEI)) {
  ggplot(as.data.frame(listOfGEI[[i]]), aes(Group, neg_del_Ct)) +
    geom_boxplot(fill='white', aes(group=Group), outliers =FALSE, position=position_dodge2(width = 1)) +
    geom_jitter(aes(color=color), size=3) +
    scale_color_manual(values = c( "#0ed1b1","#C4961A", "steelblue")) +
    ggtitle(paste0(names(listOfGEI)[i]), " Relative Expression") +
    guides(fill="none")
  ggsave(paste0(names(listOfGEI)[i], "_plot.png"))

  
}

i <- 1

for(i in seq_along(listOfGEI)) {
  copy_listOfGEI <- listOfGEI[[i]]
  # Combine rows which have the same group and color, take the mean of the neg_del_Ct
  copy_listOfGEI <- copy_listOfGEI %>%
    group_by(Group, color) %>%
    summarise(neg_del_Ct = mean(neg_del_Ct))
    res.anova <- aov(neg_del_Ct ~ Group, data=copy_listOfGEI)
    # print(summary(res.anova))
    posthoc <- TukeyHSD(res.anova)
    print((paste0(names(listOfGEI)[i], "(Mouse Level): ")))
    print(filter(as.data.frame(posthoc$Group), `p adj` < 0.05))
    print((paste0(names(listOfGEI)[i], "(Group Level): ")))
    res.anova <- aov(neg_del_Ct ~ Group, data=listOfGEI[[i]])
    # print(summary(res.anova))
    posthoc <- TukeyHSD(res.anova)
    print(filter(as.data.frame(posthoc$Group), `p adj` < 0.05))
}

#listOfGEI[[i]]$color <- sapply(mapGroup(as.integer(sapply(listOfGEI[[i]]$Samples, str_sub, -1, -1)), teqRep), backPaste, "color_")
# Asign colors to orignal exprs
Msi1_exprs$color <- sapply(mapGroup(as.integer(sapply(Msi1_exprs$Samples, str_sub, -1, -1)), teqRep), backPaste, "color_")
Lgr5_exprs$color <- sapply(mapGroup(as.integer(sapply(Lgr5_exprs$Samples, str_sub, -1, -1)), teqRep), backPaste, "color_")
Bmi1_exprs$color <- sapply(mapGroup(as.integer(sapply(Bmi1_exprs$Samples, str_sub, -1, -1)), teqRep), backPaste, "color_")
Cyp24a1_exprs$color <- sapply(mapGroup(as.integer(sapply(Cyp24a1_exprs$Samples, str_sub, -1, -1)), teqRep), backPaste, "color_")
Bcl2_exprs$color <- sapply(mapGroup(as.integer(sapply(Bcl2_exprs$Samples, str_sub, -1, -1)), teqRep), backPaste, "color_")
Ki67_exprs$color <- sapply(mapGroup(as.integer(sapply(Ki67_exprs$Samples, str_sub, -1, -1)), teqRep), backPaste, "color_")
Caspn3_exprs$color <- sapply(mapGroup(as.integer(sapply(Caspn3_exprs$Samples, str_sub, -1, -1)), teqRep), backPaste, "color_")


# Remove outliers from the data
# Remove DSS ALDH color2 with the higest value in Lgr5
Lgr5_exprs <- Lgr5_exprs[!(Lgr5_exprs$neg_del_Ct == max(Lgr5_exprs$neg_del_Ct) & Lgr5_exprs$color == "color_2"),]

# All samples from color_2 in DSS in Cyp24a1, color_1 in DSS CYP in Cyp24a1, color_2 in  NO DSS in Cyp24a1 will be removed
Cyp24a1_exprs <- Cyp24a1_exprs[!(Cyp24a1_exprs$color == "color_2" & Cyp24a1_exprs$Group == "DSS"),]
Cyp24a1_exprs <- Cyp24a1_exprs[!(Cyp24a1_exprs$color == "color_1" & Cyp24a1_exprs$Group == "DSS CYP"),]
Cyp24a1_exprs <- Cyp24a1_exprs[!(Cyp24a1_exprs$color == "color_2" & Cyp24a1_exprs$Group == "NO DSS"),]

# We will remove DSS color 2 as outlier in Bmi1
Bmi1_exprs <- Bmi1_exprs[!(Bmi1_exprs$color == "color_2" & Bmi1_exprs$Group == "DSS"),]

# We will remove the highest value of color_2 in DSS 1981 in Bcl2, the highest value of color_2 in DSS MAC in Bcl2, the highest value of color2 in no dss in Bcl2 and all of color_3 in DSS MAC in Bcl2
Bcl2_exprs <- Bcl2_exprs[!(Bcl2_exprs$color == "color_2" & Bcl2_exprs$Group == "DSS 1981"),]
Bcl2_exprs <- Bcl2_exprs[!(Bcl2_exprs$color == "color_2" & Bcl2_exprs$Group == "DSS MAC"),]
Bcl2_exprs <- Bcl2_exprs[!(Bcl2_exprs$color == "color_2" & Bcl2_exprs$Group == "NO DSS"),]
Bcl2_exprs <- Bcl2_exprs[!(Bcl2_exprs$color == "color_3" & Bcl2_exprs$Group == "DSS MAC"),]

# We will remove Color 3 in DSS Mac in Ki67 and the higest values of color 2 in NO DSS in Ki67
Ki67_exprs <- Ki67_exprs[!(Ki67_exprs$color == "color_3" & Ki67_exprs$Group == "DSS MAC"),]
Ki67_exprs <- Ki67_exprs[!(Ki67_exprs$color == "color_2" & Ki67_exprs$Group == "NO DSS"),]

# Remove color_2 in DSS from Caspn3
Caspn3_exprs <- Caspn3_exprs[!(Caspn3_exprs$color == "color_2" & Caspn3_exprs$Group == "DSS"),]


listOfGEI <- list("Msi1"=Msi1_exprs, "Lgr5"=Lgr5_exprs, "Bmi1"=Bmi1_exprs, "Cyp24a1"=Cyp24a1_exprs, "Bcl2"=Bcl2_exprs, "Ki67"=Ki67_exprs, "Casp3"=Caspn3_exprs)

for (i in seq_along(listOfGEI)) {
  ggplot(as.data.frame(listOfGEI[[i]]), aes(Group, neg_del_Ct)) +
    geom_boxplot(fill='white', aes(group=Group), outliers =FALSE, position=position_dodge2(width = 1)) +
    geom_jitter(aes(color=color), size=3) +
    scale_color_manual(values = c( "#0ed1b1","#C4961A", "steelblue")) +
    ggtitle(paste0(names(listOfGEI)[i]), " Relative Expression") +
    guides(fill="none")
  ggsave(paste0(names(listOfGEI)[i], "_plot_with_outliers_removed.png"))

  
}

for(i in seq_along(listOfGEI)) {
  copy_listOfGEI <- listOfGEI[[i]]
  # Combine rows which have the same group and color, take the mean of the neg_del_Ct
  copy_listOfGEI <- copy_listOfGEI %>%
    group_by(Group, color) %>%
    summarise(neg_del_Ct = mean(neg_del_Ct))
    res.anova <- aov(neg_del_Ct ~ Group, data=copy_listOfGEI)
    # print(summary(res.anova))
    posthoc <- TukeyHSD(res.anova)
    print((paste0(names(listOfGEI)[i], "(Mouse Level): ")))
    print(filter(as.data.frame(posthoc$Group), `p adj` < 0.05))
    print((paste0(names(listOfGEI)[i], "(Group Level): ")))
    res.anova <- aov(neg_del_Ct ~ Group, data=listOfGEI[[i]])
    # print(summary(res.anova))
    posthoc <- TukeyHSD(res.anova)
    print(filter(as.data.frame(posthoc$Group), `p adj` < 0.05))
}

library(rstatix)
library(ggpubr)


for (i in seq_along(listOfGEI)) {
  res.anova <- aov(neg_del_Ct ~ Group, data=listOfGEI[[i]])
  test <- tukey_hsd(res.anova)
  test <- test |> filter(p.adj < 0.10)
  test <- test |> add_y_position(data=listOfGEI[[i]], formula = neg_del_Ct ~ Group, fun = "max")
  ggplot(as.data.frame(listOfGEI[[i]]), aes(Group, neg_del_Ct)) +
    geom_boxplot(fill='white', aes(group=Group), outliers =FALSE, position=position_dodge2(width = 1)) +
    geom_jitter(aes(color=color), size=3) +
    scale_color_manual(values = c( "#0ed1b1","#C4961A", "steelblue")) +
    ggtitle(paste0(names(listOfGEI)[i]), " Relative Expression") +
    labs(subtitle =paste0("ANOVA Signifance: ", summary(res.anova)[[1]][["Pr(>F)"]][1]), caption = "* = p < 0.05, ** = p < 0.01, *** = p < 0.001") +
    guides(fill="none") +
    stat_pvalue_manual(test, label = "p.adj", tip.length = 0.01)
  ggsave(paste0(names(listOfGEI)[i], "_plot_with_outliers_removed_and_signifance.png"))
  
}
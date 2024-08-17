library("readxl")

setwd("D:/")
# Get all files that end with .xlsx
# files <- list.files(pattern = ".xlsx")

# get all sheets in "qPCR_Results.xlsx"
sheet_names <- excel_sheets("qPCR_Results.xlsx")
# Keep names that contain "REP"
sheet_names <- sheet_names[grep("REP", sheet_names)]

# > sheet_names
#  [1] "m_Msi1 REP"      "m_Lgr5 REP"      "m_Bmi1 REP"      "m_Cyp24a1 REP"
#  [5] "m_Bcl2 REP"      "m_Ki67 REP"      "m_Caspn3 REP"    "m_Fgfbp1 REP"
#  [9] "m_Mex3a REP"     "m_Mex3a REP (2)"

m_MsiREP_2_neg_del_Ct_GAP <- read_excel("qPCR_Results.xlsx", sheet = "m_Msi1 REP", range = "G18:J36")
m_Lgr5REP_2_neg_del_Ct_GAP <- read_excel("qPCR_Results.xlsx", sheet = "m_Lgr5 REP", range = "G18:J36")
m_Bmi1REP_2_neg_del_Ct_GAP <- read_excel("qPCR_Results.xlsx", sheet = "m_Bmi1 REP", range = "G18:J36")
m_Cyp24a1REP_2_neg_del_Ct_GAP <- read_excel("qPCR_Results.xlsx", sheet = "m_Cyp24a1 REP", range = "G18:J36")
m_Bcl2REP_2_neg_del_Ct_GAP <- read_excel("qPCR_Results.xlsx", sheet = "m_Bcl2 REP", range = "G18:J36")
m_Ki67REP_2_neg_del_Ct_GAP <- read_excel("qPCR_Results.xlsx", sheet = "m_Ki67 REP", range = "G18:J36")
m_Caspn3REP_2_neg_del_Ct_GAP <- read_excel("qPCR_Results.xlsx", sheet = "m_Caspn3 REP", range = "G18:J36")
m_Fgfbp1REP_2_neg_del_Ct_GAP <- read_excel("qPCR_Results.xlsx", sheet = "m_Fgfbp1 REP", range = "G18:J36")
m_Mex3aREP_2_neg_del_Ct_GAP <- read_excel("qPCR_Results.xlsx", sheet = "m_Mex3a REP", range = "G18:J36")
m_Mex3aREP_2_neg_del_Ct_GAP_2 <- read_excel("qPCR_Results.xlsx", sheet = "m_Mex3a REP (2)", range = "G18:J36")

GAP <- list("Msi1" = m_MsiREP_2_neg_del_Ct_GAP,
            "Lgr5" = m_Lgr5REP_2_neg_del_Ct_GAP,
            "Bmi1" = m_Bmi1REP_2_neg_del_Ct_GAP,
            "Cyp24a1" = m_Cyp24a1REP_2_neg_del_Ct_GAP,
            "Bcl2" = m_Bcl2REP_2_neg_del_Ct_GAP,
            "Ki67" = m_Ki67REP_2_neg_del_Ct_GAP,
            "Casp3" = m_Caspn3REP_2_neg_del_Ct_GAP,
            "Fgfbp1" = m_Fgfbp1REP_2_neg_del_Ct_GAP,
            "Mex3a" = m_Mex3aREP_2_neg_del_Ct_GAP,
            "Mex3a2" = m_Mex3aREP_2_neg_del_Ct_GAP_2)

# [1] "Msi1 Lgr5"
# [1] "Msi1 Bmi1"
# [1] "Lgr5 Msi1"

# [1] "Cyp24a1

# [1] "Bcl2 Ki67"
# [1] "Bcl2 Caspn3"

# [1] "Fgfbp1 Mex3a"


# For each sheet plot as scatter plot with each group in different x value and color for technical replicates which are the row names
library(ggplot2)
library(tidyr)

for (gene in names(GAP)) {
  x <- as.data.frame(GAP[[gene]])
  rownames(x) <- x$`...1`
  x$`...1` <- NULL
  GAP[[gene]] <- x
}

# For each gene
for (gene in names(GAP)) {
  x <- GAP[[gene]] |>
tibble::rownames_to_column("Group") |> 
  tidyr::pivot_longer(cols=!Group, names_to = "Samples", values_to = "neg_del_Ct")


  # remove NA values
    x <- x[!is.na(x$neg_del_Ct),]
    # calculate the mean and standard deviation for each group
    x <- x |>
      dplyr::group_by(Group) |>
      dplyr::mutate(mean_var = mean(neg_del_Ct), sd_var = sd(neg_del_Ct)) |>
      dplyr::ungroup()
    print(x)
  # Have the x values tiles and spread out
  # Also make a line on the mean of each group with a bar for standard deviation 
    ggplot(x, aes(x = Group, y = neg_del_Ct, color = Samples)) +
    geom_point() +
    labs(title = paste(gene, "Experiment GAP Expression")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    # Calculte the mean and 0.5 * sd for each group
    geom_errorbar(aes(ymin=mean_var-(1*sd_var), ymax=mean_var+(1*sd_var)), width=.2, color = "black") +
    # plot the mean of each group
    geom_line(aes(y = mean_var), color = "black")
    ggsave(paste0(gene, "_GAP_Expression.png"))
}

# 
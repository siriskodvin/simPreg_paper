

# Example use of the simPreg package


library(simPreg)


# ---- 1. Simulate pregnancy outcome proportions ------------------------------

# Run simPregProp() using package-provided default hazard rate and hazard ratio
# vectors
df_prop <- simPregProp()

head(df_prop)


# ---- 2. Sample pregnancy data -----------------------------------------------

# Sample a dataset with grouped counts (expand = FALSE)
set.seed(9)
df_samp_grp <- simPregSamp(df = df_prop,
                           n = 100000,
                           expand = FALSE)
head(df_samp_grp)

# Sample a dataset with one row per pregnancy (expand = TRUE)
set.seed(9)
df_samp_one <- simPregSamp(df = df_prop,
                           n = 100000,
                           expand = TRUE)
df_samp_one[1:10,]


# ---- 3. Specify a hazard ratio for preterm birth ----------------------------

# Set a custom hazard ratio for preterm live birth (before week 37/day 259)
# following exposure
hr_preterm <- c(rep(5, 258), rep(1, 301-258))
df_prop_pt <- simPregProp(hr.spont.livebirth = hr_preterm,
                          hr.nonspont.livebirth = hr_preterm)

# Sample a dataset and compare the proportion of preterm births among exposed
# and unexposed pregnancies
set.seed(10)
df_samp_pt <- simPregSamp(df = df_prop_pt,
                          n = 100000,
                          expand = TRUE)
prop_preterm_exposed <- mean(df_samp_pt$GA[!is.na(df_samp_pt$ExpGA)] <= 258)
prop_preterm_unexposed <- mean(df_samp_pt$GA[is.na(df_samp_pt$ExpGA)] <= 258)
cat("Proportion of preterm births among exposed pregnancies:",
    round(prop_preterm_exposed, 3), "\n")
cat("Proportion of preterm births among unexposed pregnancies:",
    round(prop_preterm_unexposed, 3), "\n")


# ---- 4. Specify a trimester-specific hazard ratio ---------------------------

# Set a custom hazard ratio for late miscarriage/stillbirth following exposure
# during varying trimesters
hr_fetal_death <- c(rep(3, 90), rep(2, 105), rep(1, (301 - (90 + 105))))
df_prop_fd <- simPregProp(hr.late.miscarriage.stillbirth = hr_fetal_death)

# Sample a dataset and compare late miscarriage/stillbirth by trimester of
# exposure in the no-effect (default) and exposure-effect datasets
set.seed(11)
df_samp_fd <- simPregSamp(df = df_prop_fd,
                          n = 100000,
                          expand = TRUE)
trimester_breaks <- c(-Inf, 90, 195, Inf)
trimester_labels <- c("1st trimester", "2nd trimester", "3rd trimester")
no_effect_fetal_death <- df_samp_one[df_samp_one$Outcome ==
                                       "late_miscarriage_stillbirth" &
                                       !is.na(df_samp_one$ExpGA),]
exposure_effect_fetal_death <- df_samp_fd[df_samp_fd$Outcome ==
                                            "late_miscarriage_stillbirth" &
                                            !is.na(df_samp_fd$ExpGA),]
cat("\nLate miscarriage/stillbirth by trimester of exposure",
    "in the no-effect (default) dataset:\n")
print(table(cut(no_effect_fetal_death$ExpGA, breaks = trimester_breaks,
                labels = trimester_labels)))
cat("\nLate miscarriage/stillbirth by trimester of exposure",
    "in the exposure-effect dataset:\n")
print(table(cut(exposure_effect_fetal_death$ExpGA, breaks = trimester_breaks,
                labels = trimester_labels)))

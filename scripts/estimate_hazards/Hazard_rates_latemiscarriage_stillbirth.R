
library(ggplot2)
library(survival)
library(dplyr)



#### READ EXTRACT FROM THE MEDICAL BIRTH REGISTRY AND SUBSET ####

# Read
df_mbrn <- readRDS(file = "df_mbrn.rds")

# Subset to relevant time period
df_mbrn <- subset(df_mbrn, subset = FDATO <= as.Date("2019-12-31"))

# Randomly sample one pregnancy per mother
length(unique(df_mbrn$W20_1290_LOPENR_PERSON))
set.seed(424)
df_mbrn <- df_mbrn[unlist(lapply(split(seq_len(nrow(df_mbrn)), df_mbrn$W20_1290_LOPENR_PERSON),
                                 function(x){ifelse(length(x)>1, sample(x, 1), x)})),]



#### DESCRIPTIVES ####

# Number or late miscarriages/stillbirths and live births
table(df_mbrn$sb)

# Number of late miscarriages/stillbirths by week
table(df_mbrn$ga_wk[which(df_mbrn$sb == 1)])

# Adjust for plotting
df_mbrn$early <- ifelse(df_mbrn$ga_wk < 22, 1, 0)
tab <- tapply(df_mbrn$sb == 1, df_mbrn$ga_wk, sum)
df_mbrn$sb_adj <- df_mbrn$sb
df_mbrn$sb_adj[which(df_mbrn$ga_wk %in% as.numeric(names(tab)[which(tab < 5)]))] <- 0

# Plot late miscarriage/stillbirth distribution by gestational age
pl_sb <- ggplot(subset(df_mbrn, subset = sb_adj == 1), aes(x = ga_wk, color = factor(early), fill = factor(early))) +
  geom_bar(position = position_dodge(width = 0.8), linewidth = 0.1) +
  scale_color_manual(values = c("1" = "black", "0" = "black"),
                     guide = "none") +
  scale_fill_manual(values = c("1" = "#CDDCE1", "0" = "#6796A7"),
                    breaks = c("1", "0"),
                    labels = c("Lower reporting completeness", "Higher reporting completeness")) +
  scale_x_continuous(breaks = seq(from = 0, to = max(df_mbrn$ga_wk), by = 5)) +
  labs(title = "Late miscarriage/stillbirth distribution", x = "Gestational age (weeks)", y = NULL, fill = NULL, alpha = NULL) +
  guides(fill = guide_legend(override.aes = list(linewidth = 4))) +
  theme_minimal(base_size = 13)

# Save plot as pdf
ggsave(file = "pl_sb.pdf",
       plot = pl_sb, width = 6, height = 2.1, dpi = 600)



#### ESTIMATE LATE MISCARRIAGE/STILLBIRTH HAZARDS ####

# Survival analysis of GA of late miscarriage/stillbirth
km_ga_sb <- survfit(Surv(ga, sb) ~ 1, data = df_mbrn)

# Construct estimate dataframe, include frequencies
df_est <- data.frame(GA = seq(from = 0, to = max(km_ga_sb$time)))
tmp <- as.data.frame(table(df_mbrn$ga[which(df_mbrn$sb == 1)]))
tmp$Var1 <- as.numeric(as.character(tmp$Var1))
df_est <- left_join(df_est, tmp, by = c("GA" = "Var1"))
df_est$Freq[which(is.na(df_est$Freq))] <- 0

# Get cumulative hazard estimates
df_est <- left_join(df_est, data.frame(GA = km_ga_sb$time, cumhaz = km_ga_sb$cumhaz),
                    by = c("GA" = "GA"))

# Approximate missing cumulative hazard estimates
df_est$cumhaz_appr <- df_est$cumhaz
df_est$cumhaz_appr[which(df_est$GA < 12*7)] <- 0
df_est$cumhaz_appr <- approx(x = df_est$GA[which(!is.na(df_est$cumhaz_appr))],
                             y = df_est$cumhaz_appr[which(!is.na(df_est$cumhaz_appr))],
                             xout = df_est$GA)$y

# Calculate hazards
df_est$haz <- c(df_est$cumhaz_appr[1], diff(df_est$cumhaz_appr))

# Subset
haz.late.miscarriage.stillbirth.default <- df_est$haz[2:302]

# Save data as .rds and .rda
saveRDS(haz.late.miscarriage.stillbirth.default, file = "haz.late.miscarriage.stillbirth.default.rds")
save(haz.late.miscarriage.stillbirth.default, file = "haz.late.miscarriage.stillbirth.default.rda")



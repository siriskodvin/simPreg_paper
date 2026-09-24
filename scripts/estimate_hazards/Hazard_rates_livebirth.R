
library(ggplot2)
library(survival)
library(dplyr)



#### READ EXTRACT FROM THE MEDICAL BIRTH REGISTRY AND SUBSET ####

# Read
df_mbrn <- readRDS(file = "df_mbrn.rds")

# Subset to relevant time period
df_mbrn <- subset(df_mbrn, subset = (FDATO >= as.Date("2016-01-01") & FDATO <= as.Date("2019-12-31")))

# Exlude stillbirths
df_mbrn <- subset(df_mbrn, subset = sb == 0)

# Randomly sample one pregnancy per mother
set.seed(9)
df_mbrn <- df_mbrn[unlist(lapply(split(seq_len(nrow(df_mbrn)), df_mbrn$W20_1290_LOPENR_PERSON),
                            function(x){ifelse(length(x)>1, sample(x, 1), x)})),]

# Define spontaneous and non-spontaneous births
df_mbrn$birth_spont <- ifelse(df_mbrn$FSTART == 1, 1, 0)
df_mbrn$birth_nonspont <- ifelse(df_mbrn$FSTART == 1, 0, 1)



#### DESCRIPTIVES ####

# Code birth mode factor
df_mbrn$birth_mode <- factor(ifelse(df_mbrn$birth_spont == 1, "Spontaneous", "Non-spontaneous"),
                             levels = c("Spontaneous", "Non-spontaneous"))

# Number of spontaneous and non-spontaneous births
table(df_mbrn$birth_mode)

# Adjust for plotting
tab <- as.data.frame(table(df_mbrn$ga, df_mbrn$birth_mode))
tab$Var1 <- as.numeric(as.character(tab$Var1))
tab$Var2 <- as.character(tab$Var2)
colnames(tab) <- c("ga", "birth_mode", "Freq")
tab <- subset(tab, subset = Freq < 5)
df_mbrn$ga_adj <- df_mbrn$ga
df_mbrn$ga_adj[which(paste0(df_mbrn$ga, df_mbrn$birth_mode) %in%
                       paste0(tab$ga, tab$birth_mode))] <- NA

# Plot live birth distribution by gestational age and birth mode
pl_lb <- ggplot(df_mbrn, aes(x = ga_adj, fill = birth_mode)) +
  geom_bar(position = position_dodge(width = 1)) +
  scale_fill_manual(values = c("Spontaneous" = "#6A6B88", "Non-spontaneous" = "#EE756A")) +
  scale_x_continuous(breaks = seq(from = floor(min(df_mbrn$ga)/10)*10, to = max(df_mbrn$ga), by = 10)) +
  labs(title = "Live birth distribution", x = "Gestational age (days)", y = NULL, fill = NULL) +  guides(fill = guide_legend(override.aes = list(linewidth = 4))) +
  theme_minimal(base_size = 13)

# Save plot as pdf
ggsave(file = "pl_lb.pdf",
       plot = pl_lb, width = 7.9, height = 2.1, dpi = 600)



#### ESTIMATE SPONTANEOUS LIVE BIRTH HAZARDS ####

# Survival analysis of GA of spontaneous live birth
km_ga_lb_spont <- survfit(Surv(ga, birth_spont) ~ 1, data = df_mbrn)

# Construct estimate dataframe, include frequencies
df_est <- data.frame(GA = seq(from = 0, to = max(km_ga_lb_spont$time)))
tmp <- as.data.frame(table(df_mbrn$ga[which(df_mbrn$birth_spont == 1)]))
tmp$Var1 <- as.numeric(as.character(tmp$Var1))
df_est <- left_join(df_est, tmp, by = c("GA" = "Var1"))
df_est$Freq[which(is.na(df_est$Freq))] <- 0

# Get cumulative hazard estimates
df_est <- left_join(df_est, data.frame(GA = km_ga_lb_spont$time, cumhaz = km_ga_lb_spont$cumhaz),
                    by = c("GA" = "GA"))

# Approximate missing cumulative hazard estimates
df_est$cumhaz_appr <- df_est$cumhaz
df_est$cumhaz_appr[which(df_est$GA < 22*7)] <- 0
df_est$cumhaz_appr <- approx(x = df_est$GA[which(!is.na(df_est$cumhaz_appr))],
                             y = df_est$cumhaz_appr[which(!is.na(df_est$cumhaz_appr))],
                             xout = df_est$GA)$y

# Calculate hazards
df_est$haz <- c(df_est$cumhaz_appr[1], diff(df_est$cumhaz_appr))

# Subset
haz.spont.livebirth.default <- df_est$haz[2:302]

# Save data as .rds and .rda
saveRDS(haz.spont.livebirth.default, file = "haz.spont.livebirth.default.rds")
save(haz.spont.livebirth.default, file = "haz.spont.livebirth.default.rda")



#### ESTIMATE NON-SPONTANEOUS LIVE BIRTH HAZARDS ####

# Survival analysis of GA of non-spontaneous live birth
km_ga_lb_nonspont <- survfit(Surv(ga, birth_nonspont) ~ 1, data = df_mbrn)

# Construct estimate dataframe, include N for each GA
df_est <- data.frame(GA = seq(from = 0, to = max(km_ga_lb_nonspont$time)))
tmp <- as.data.frame(table(df_mbrn$ga[which(df_mbrn$birth_nonspont == 1)]))
tmp$Var1 <- as.numeric(as.character(tmp$Var1))
df_est <- left_join(df_est, tmp, by = c("GA" = "Var1"))
df_est$Freq[which(is.na(df_est$Freq))] <- 0

# Get cumulative hazard estimates
df_est <- left_join(df_est, data.frame(GA = km_ga_lb_nonspont$time, cumhaz = km_ga_lb_nonspont$cumhaz),
                    by = c("GA" = "GA"))

# Approximate missing survival estimates
df_est$cumhaz_appr <- df_est$cumhaz
df_est$cumhaz_appr[which(df_est$GA < 22*7)] <- 0
df_est$cumhaz_appr <- approx(x = df_est$GA[which(!is.na(df_est$cumhaz_appr))],
                             y = df_est$cumhaz_appr[which(!is.na(df_est$cumhaz_appr))],
                             xout = df_est$GA)$y

# Calculate hazards
df_est$haz <- c(df_est$cumhaz_appr[1], diff(df_est$cumhaz_appr))

# Subset
haz.nonspont.livebirth.wo1 <- df_est$haz[2:302]

# Set last hazard to 1
haz.nonspont.livebirth.default[length(haz.nonspont.livebirth.default)] <- 1

# Save data as .rds and .rda
saveRDS(haz.nonspont.livebirth.default, file = "haz.nonspont.livebirth.default.rds")
save(haz.nonspont.livebirth.default, file = "haz.nonspont.livebirth.default.rda")




library(ggplot2)
library(survival)
library(dplyr)



#### READ EXTRACT FROM THE MEDICAL BIRTH REGISTRY AND SUBSET ####

# Read
df_mbrn <- readRDS(file = "df_mbrn.rds")

# Subset to relevant time period
df_mbrn <- subset(df_mbrn, subset = (FDATO >= as.Date("2020-01-01") & FDATO <= as.Date("2023-12-31")))

# Randomly sample one pregnancy per mother
length(unique(df_mbrn$W20_1290_LOPENR_PERSON))
set.seed(150)
df_mbrn <- df_mbrn[unlist(lapply(split(seq_len(nrow(df_mbrn)), df_mbrn$W20_1290_LOPENR_PERSON),
                                 function(x){ifelse(length(x)>1, sample(x, 1), x)})),]



#### READ EXTRACT FROM SYSVAK AND SUBSET ####

# Read SYSVAK and clean date variable
sysvak <- read.csv(file = "sysvak.csv", sep = ";")
sysvak$KONSULTASJONSDATO <- as.Date(sysvak$KONSULTASJONSDATO)

# Subset to relevant time period and individuals (based on current mbrn subset)
df_sysvak <- subset(sysvak, subset = (KONSULTASJONSDATO >= min(df_mbrn$lmpdate) &
                                        KONSULTASJONSDATO <= max(df_mbrn$FDATO)))
df_sysvak <- subset(df_sysvak, subset = w20_1290_lopenr_person %in% df_mbrn$W20_1290_LOPENR_PERSON)

# Inspect vaccination codes using full SYSVAK extract, identify C19 vaccination codes
table(sysvak$VAKSINEKODE)
table(sysvak$PREPARATBESKRIVELSE, sysvak$VAKSINEKODE)
vcode <- c("ASZ03", "BNT03", "CBA01", "CBA45", "CBB15", "CSH03", "CUR03", "JAN03", "MOD03",
           "NUV03", "SBA01", "SBA45", "SIN03", "VAL03", "XCO03", "XCV03", "XXC03", "XXS03")

# Subset to C19 vaccination codes in data subset
df_sysvak <- subset(df_sysvak, subset = VAKSINEKODE %in% vcode)
vcode <- subset(vcode, subset = vcode %in% df_sysvak$VAKSINEKODE)

# Inspect details for included vaccination codes
for(code in vcode){
  print("------------------------------------------------------------------------------")
  print(code)
  print(unique(df_sysvak$PREPARATBESKRIVELSE[which(df_sysvak$VAKSINEKODE == code)]))
  print(summary(df_sysvak$KONSULTASJONSDATO[which(df_sysvak$VAKSINEKODE == code)]))
  print("------------------------------------------------------------------------------")
}

# Exclude individuals with an unrealistic vaccine date (before the first C19 vaccine in Norway) from both mbrn and sysvak subsets
ids <- unique(df_sysvak$w20_1290_lopenr_person[which(df_sysvak$KONSULTASJONSDATO < as.Date("2020-12-27"))])
df_mbrn <- subset(df_mbrn, subset = !(W20_1290_LOPENR_PERSON %in% ids))
df_sysvak <- subset(df_sysvak, subset = !(w20_1290_lopenr_person %in% ids))
sort(df_sysvak$KONSULTASJONSDATO)[1:10]

# Remove duplicates for ID and vaccination date - ignoring differing vaccine codes
rownames(df_sysvak) <- NULL
if(length(which(duplicated(df_sysvak[,c("w20_1290_lopenr_person", "KONSULTASJONSDATO")]))) > 0){
  df_sysvak <- df_sysvak[-which(duplicated(df_sysvak[,c("w20_1290_lopenr_person", "KONSULTASJONSDATO")])),]
}

# Subset to keep only relevant variables
df_sysvak <- subset(df_sysvak, select = c(w20_1290_lopenr_person, KONSULTASJONSDATO))

# Sort and number vaccine doses, and reshape to wide format
df_sysvak <- df_sysvak[order(df_sysvak$KONSULTASJONSDATO),]
df_sysvak <- df_sysvak[order(df_sysvak$w20_1290_lopenr_person),]
df_sysvak$dose <- ave(seq_along(df_sysvak$KONSULTASJONSDATO), df_sysvak$w20_1290_lopenr_person, FUN = seq_along)
df_sysvak <- reshape(df_sysvak, idvar = "w20_1290_lopenr_person", timevar = "dose", direction = "wide")



#### MERGE MBRN AND SYSVAK SUBSETS ####

# Merge mbrn and sysvak subsets
length(unique(df_mbrn$W20_1290_LOPENR_PERSON))
length(unique(df_sysvak$w20_1290_lopenr_person))
df <- left_join(df_mbrn, df_sysvak, by = c("W20_1290_LOPENR_PERSON" = "w20_1290_lopenr_person"))

# Remove vaccines before or after pregnancy
cinds <- which(startsWith(colnames(df), "KONSULTASJONSDATO"))
for(i in 1:length(cinds)){
  cind <- cinds[i]
  colnames(df)[cind] <- "vaccdate"
  df$vaccdate[which(df$vaccdate <= df$lmp | df$vaccdate >= df$FDATO)] <- NA
  colnames(df)[cind] <- paste0("vaccdate_", i)
}

# List and tabulate ordering of vaccine doses
df$dose_rec <- apply(df[,cinds], 1, function(x){
  nos <- NULL
  for(i in 1:length(cinds)){
    if(!is.na(x[i])){
      nos <- paste0(nos, i)
    }
  }
  if(is.null(nos)){
    nos <- "0"
  }
  return(nos)
})
tab <- table(df$dose_rec)

# Loop to fix ordering
while(length(names(tab)[!grepl("^[01]", names(tab))]) > 0){

  # Update each dose
  cinds <- which(startsWith(colnames(df), "vaccdate"))
  for(cind in cinds[1:(length(cinds)-1)]){
    rinds <- which(is.na(df[,cind]))
    df[rinds,cind] <- df[rinds,(cind+1)]
    df[rinds,(cind+1)] <- NA
  }
  
  # Update ordering
  df$dose_rec <- apply(df[,cinds], 1, function(x){
    nos <- NULL
    for(i in 1:length(cinds)){
      if(!is.na(x[i])){
        nos <- paste0(nos, i)
      }
    }
    if(is.null(nos)){
      nos <- "0"
    }
    return(nos)
  })
  tab <- table(df$dose_rec)
}

# Keep only first vaccination during pregnancy and calculate vaccination timing
df <- df[,1:which(colnames(df) == "vaccdate_1")]
df$vacctime <- as.numeric(df$vaccdate_1 - df$lmpdate)
df$vaccweek <- floor(df$vacctime/7)

# Number or vaccinated and unvaccinated
df$vacc_bin <- ifelse(!is.na(df$vaccdate_1), 1, 0)
table(df$vacc_bin)



#### DESCRIPTIVES ####

# Plot vaccination timing
pl_vacc <- ggplot(df, aes(x = vaccweek)) +
  geom_bar(position = position_dodge(width = 0.8), color = "black", fill = "#824F73", linewidth = 0.1) +
  scale_x_continuous(breaks = seq(from = 0, to = max(df$vaccweek, na.rm = TRUE), by = 5)) +
  labs(title = "Vaccination timing in pregnancy", x = "Gestational age (weeks)", y = NULL, fill = NULL, alpha = NULL) +
  theme_minimal(base_size = 13)

# Save plot as pdf
ggsave(file = "pl_vacc.pdf",
       plot = pl_vacc, width = 4.6, height = 2.1, dpi = 600)



#### ESTIMATE VACCINATION HAZARDs FOR THE COMPLETE PERIOD ####

# Survival analysis of GA of vaccination
df$stoptime <- ifelse(!is.na(df$vacctime), df$vacctime, df$ga)
km_ga_vacc <- survfit(Surv(stoptime, vacc_bin) ~ 1, data = df)

# Construct estimate dataframe, include frequencies
df_est <- data.frame(GA = seq(from = 0, to = max(km_ga_vacc$time)))
tmp <- as.data.frame(table(df$vacctime[which(df$vacc_bin == 1)]))
tmp$Var1 <- as.numeric(as.character(tmp$Var1))
df_est <- left_join(df_est, tmp, by = c("GA" = "Var1"))
df_est$Freq[which(is.na(df_est$Freq))] <- 0

# Get cumulative hazard estimates
df_est <- left_join(df_est, data.frame(GA = km_ga_vacc$time, cumhaz = km_ga_vacc$cumhaz),
                    by = c("GA" = "GA"))

# Approximate missing cumulative hazard estimates
df_est$cumhaz_appr <- df_est$cumhaz
df_est$cumhaz_appr[1] <- 0
df_est$cumhaz_appr <- approx(x = df_est$GA[which(!is.na(df_est$cumhaz_appr))],
                             y = df_est$cumhaz_appr[which(!is.na(df_est$cumhaz_appr))],
                             xout = df_est$GA)$y

# Calculate hazards
df_est$haz <- c(df_est$cumhaz_appr[1], diff(df_est$cumhaz_appr))

# Subset
haz.exposure.default <- df_est$haz[2:302]

# Save data as .rds and .rda
saveRDS(haz.exposure.default, file = "haz.exposure.default.rds")
save(haz.exposure.default, file = "haz.exposure.default.rda")



#### ESTIMATE VACCINATION HAZARDS SPLIT ON DIFFERENT PERIODS ####

# Define periods by LMP cohorts
lmp_cohorts <- as.Date(c("2020-10-01", "2021-01-01", "2021-03-01",
                         "2021-05-01", "2021-07-01", "2021-09-01",
                         "2021-11-01", "2022-01-01", "2023-01-01"))

# Back up full dataframe
df_full <- df

# Initiate dataframe
haz.exposure.lmpcohorts <- data.frame(GA = c(1:301))

# Loop through lmp cohorts
for(i in 1:8){
  
  # Subset full dataset
  df <- subset(df_full, subset = (lmpdate >= lmp_cohorts[i] & lmpdate < lmp_cohorts[i+1]))
  
  # Survival analysis of GA of vaccination
  df$stoptime <- ifelse(!is.na(df$vacctime), df$vacctime, df$ga)
  km_ga_vacc <- survfit(Surv(stoptime, vacc_bin) ~ 1, data = df)
  
  # Construct estimate dataframe, include frequencies
  df_est <- data.frame(GA = seq(from = 0, to = max(km_ga_vacc$time)))
  tmp <- as.data.frame(table(df$vacctime[which(df$vacc_bin == 1)]))
  tmp$Var1 <- as.numeric(as.character(tmp$Var1))
  df_est <- left_join(df_est, tmp, by = c("GA" = "Var1"))
  df_est$Freq[which(is.na(df_est$Freq))] <- 0
  
  # Get cumulative hazard estimates
  df_est <- left_join(df_est, data.frame(GA = km_ga_vacc$time, cumhaz = km_ga_vacc$cumhaz),
                      by = c("GA" = "GA"))
  
  # Approximate missing cumulative hazard estimates
  df_est$cumhaz_appr <- df_est$cumhaz
  df_est$cumhaz_appr[1] <- 0
  df_est$cumhaz_appr <- approx(x = df_est$GA[which(!is.na(df_est$cumhaz_appr))],
                               y = df_est$cumhaz_appr[which(!is.na(df_est$cumhaz_appr))],
                               xout = df_est$GA)$y
  
  # Calculate hazards
  df_est$haz <- c(df_est$cumhaz_appr[1], diff(df_est$cumhaz_appr))
  
  # Add to dataframe
  haz.exposure.lmpcohorts <- left_join(haz.exposure.lmpcohorts, subset(df_est,
                                                                       select = c(GA, haz),
                                                                       subset = GA <= 301),
                                       by = c("GA" = "GA"))
  colnames(haz.exposure.lmpcohorts)[ncol(haz.exposure.lmpcohorts)] <-
    paste0("lmp_", lmp_cohorts[i], "_", (lmp_cohorts[i+1]-1))
}

# Replace NAs
haz.exposure.lmpcohorts[is.na(haz.exposure.lmpcohorts)] <- 0

# Save data as .rds and .rda
saveRDS(haz.exposure.lmpcohorts, file = "haz.exposure.lmpcohorts.rds")
save(haz.exposure.lmpcohorts, file = "haz.exposure.lmpcohorts.rda")



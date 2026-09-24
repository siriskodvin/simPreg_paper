
library(simPreg)
library(survival)
library(ggplot2)



#### LOOP TO SIMULATE ####

# Initiate list to store results
res_list <- list()

t <- Sys.time()
set.seed(20)
for(hr in c(0.5, 1, 2)){
  
  # Input
  hr.spont.livebirth.base <- hr
  hr.spont.livebirth <- c(rep(hr.spont.livebirth.base, 258), rep(1, (301 - 258)))
  hr.nonspont.livebirth.base <- hr
  hr.nonspont.livebirth <- c(rep(hr.nonspont.livebirth.base, 258), rep(1, (301 - 258)))
  hr.late.miscarriage.stillbirth.base <- 1
  hr.late.miscarriage.stillbirth <- rep(hr.late.miscarriage.stillbirth.base, 301)

  # Function call simPreg
  df.prop <- simPregProp(hr.spont.livebirth = hr.spont.livebirth,
                         hr.nonspont.livebirth = hr.nonspont.livebirth,
                         hr.late.miscarriage.stillbirth = hr.late.miscarriage.stillbirth)
  sumProp <- sum(df.prop$Prop)
  
  # Vary sample size 1000 times
  for(i in 1:1000){
    n <- sample(c(1000:50000), size = 1, replace = FALSE)
    
    # Print status
    if((i/100) %in% c(1:10)){
      print(Sys.time())
      print(paste0("HR = ", hr))
      print(paste0("n = ", n))
      print(paste0("i = ", i))
      print("--------------------------------")
    }
    
    # Simulate
    df.samp <- simPregSamp(n = n,
                           df = df.prop,
                           expand = TRUE)
    
    # Add random IDs
    df <- df.samp
    rownames(df) <- NULL
    df$id <- paste0("SID_", seq(from = 1, to = nrow(df), by = 1))
    
    # Add delayed entry and subset to exclude short pregnancies and late miscarriage/stillbirth
    delayed_entry <- 20*7
    df$entry <- delayed_entry
    df <- subset(df, subset = (GA > delayed_entry & Outcome != "late_miscarriage_stillbirth"))
    
    # Define preterm variables
    df$vacc_bin_pre <- ifelse(!is.na(df$ExpGA) & df$ExpGA < 258, 1, 0)
    df$preterm <- ifelse(df$GA < 259, 1, 0)
    df$ga_cens <- pmin(df$GA, 258)
    
    # Re-structure and subset
    df <- df[,c("id", "GA", "ExpGA", "entry", "ga_cens", "vacc_bin_pre", "preterm")]
    
    # Find numbers
    n_vacc <- length(which(df$vacc_bin_pre == 1))
    n_pre <- length(which(df$preterm == 1))
    n_pre_vacc <- length(which(df$vacc_bin_pre == 1 & df$preterm == 1))
    
    # Cox, time-varying
    df_sub <- tmerge(data1 = df,
                     data2 = df,
                     id = id,
                     tstart = entry,
                     tstop = ga_cens,
                     event = event(ga_cens, preterm),
                     vacc = tdc(ExpGA))
    mod <- coxph(Surv(tstart, tstop, event) ~ vacc, data = df_sub)
    est_log <- summary(mod)$coefficients["vacc", "coef"]
    est_exp <- summary(mod)$coefficients["vacc", "exp(coef)"]
    se_log <- summary(mod)$coefficients["vacc", "se(coef)"]
    cil_log <- suppressMessages(confint(mod))["vacc", "2.5 %"]
    ciu_log <- suppressMessages(confint(mod))["vacc", "97.5 %"]

    # Save results in dataframe and add to list
    df_res_new <- data.frame(hr = hr,
                             samplesize = n,
                             sumProp = sumProp,
                             n = nrow(df),
                             n_vacc = n_vacc,
                             n_pre = n_pre,
                             n_pre_vacc = n_pre_vacc,
                             est_log = est_log,
                             est_exp = est_exp,
                             se_log = se_log,
                             cil_log = cil_log,
                             ciu_log = ciu_log)
    res_list[[length(res_list)+1]] <- df_res_new
    
    # Clean
    rm(mod)
  }
}
print(Sys.time() - t)

# Collapse result list into dataframe
df_res <- do.call(rbind, res_list)



#### SUMMARIZE/PLOT ####

# Adjust n (by thousands)
df_res$n_1000 <- df_res$n/1000

# Adjust samplesize (by thousands)
df_res$samplesize_1000 <- df_res$samplesize/1000

# Make HR factor variable
df_res$hr_fac <- factor(df_res$hr, levels = c("2", "1", "0.5"))

# Calculate exp CIs
df_res$cil_exp <- exp(df_res$cil_log)
df_res$ciu_exp <- exp(df_res$ciu_log)

# Plot estimates with CIs over sample size, for three different true HR values
pl_prec <- ggplot(data = df_res, aes(x = n_1000, y = est_exp, ymin = cil_exp, ymax = ciu_exp, color = hr_fac, fill = hr_fac, shape = hr_fac)) +
  geom_hline(aes(yintercept = hr, color = hr_fac), linetype = 6, linewidth = 0.5, show.legend = FALSE) +
  geom_point(size = 0.5, show.legend = FALSE) +
  geom_errorbar(linewidth = 0.3, show.legend = FALSE) +
  scale_color_manual(values = c("0.5" = "#6796A7",
                                "1" = "#EE756A",
                                "2" = "#824F73")) +
  scale_shape_manual(values = c("0.5" = 16,
                                "1" = 15,
                                "2" = 17)) +
  scale_y_continuous(trans = scales::log_trans(base = 2),
                     breaks = c(0.0625, 0.125, 0.25, 0.5, 1, 2, 4),
                     labels = c("0.0625", "0.125", "0.25", "0.5", "1", "2", "4")) +
  labs(title = NULL, x = "Sample size (thousands)", y = "Effect estimate", color = NULL) +
  theme_minimal(base_size = 13) +
  theme(panel.border =  element_rect(color = "black", fill = NA, linewidth = 0.5))

# Save as pdf
ggsave("pl_prec.pdf", plot = pl_prec, width = 5.63, height = 2.25)



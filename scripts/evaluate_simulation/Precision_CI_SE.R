
library(simPreg)
library(survival)
library(dplyr)
library(ggplot2)



#### LOOP TO SIMULATE ####

# Initiate list to store results
res_list <- list()

# Save time information
t <- Sys.time()

# Vary hazard ratios
set.seed(51)
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
  
  # Vary sample size
  for(samplesize in c(1200, 3000, 6000, 12000, 30000, 60000)){
    
    # Set repetitions
    for(i in 1:1000){
      
      # Print status
      if((i/100) %in% c(1:10)){
        print(Sys.time())
        print(paste0("HR = ", hr))
        print(paste0("n = ", samplesize))
        print(paste0("i = ", i))
        print("--------------------------------")
      }
      
      # Simulate
      df.samp <- simPregSamp(n = samplesize,
                             df = df.prop,
                             expand = TRUE)
      
      # Add random IDs
      df <- df.samp
      df$id <- paste0("SID_", seq(from = 1, to = nrow(df), by = 1))
      
      # Add delayed entry and subset to exclude short pregnancies and late miscarriage/stillbirth
      delayed_entry <- 20*7
      df$entry <- delayed_entry
      df <- subset(df, subset = (GA > delayed_entry & Outcome != "late_miscarriage_stillbirth"))
      rownames(df) <- NULL
      
      # Subset to get fixed sample size
      n <- samplesize*5/6
      df <- df[sample(1:nrow(df), size = n, replace = FALSE),]
      rownames(df) <- NULL
      
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
                               sumProp = sumProp,
                               n = n,
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
}
print(Sys.time() - t)

# Collapse result list into dataframe
df_res <- do.call(rbind, res_list)



#### SUMMARIZE/PLOT ####

# Check coverage of CIs
df_res$coverage <- ifelse(df_res$cil_log < log(df_res$hr) & df_res$ciu_log > log(df_res$hr), 1, 0)

# Aggregate results
df_res_tab <- aggregate(coverage ~ hr + n, data = df_res, FUN = mean)
tmp <- aggregate(est_log ~ hr + n, data = df_res, FUN = mean)
colnames(tmp)[which(colnames(tmp) == "est_log")] <- "est_log_mean"
df_res_tab <- left_join(df_res_tab, tmp, by = c("hr" = "hr", "n" = "n"))
tmp <- aggregate(se_log ~ hr + n, data = df_res, FUN = mean)
colnames(tmp)[which(colnames(tmp) == "se_log")] <- "se_log_mean"
df_res_tab <- left_join(df_res_tab, tmp, by = c("hr" = "hr", "n" = "n"))
tmp <- aggregate(est_log ~ hr + n, data = df_res, FUN = sd)
colnames(tmp)[which(colnames(tmp) == "est_log")] <- "se_log_sim"
df_res_tab <- left_join(df_res_tab, tmp, by = c("hr" = "hr", "n" = "n"))

# Coverage table
coverage_table <- reshape(
  df_res_tab[, c("hr", "n", "coverage")],
  idvar = "hr",
  timevar = "n",
  direction = "wide"
)

# Remove the HR = 0.5 and n = 1,000 scenario due to convergence issues
df_res_tab <- subset(df_res_tab, subset = (hr != 0.5 | n != 1000))

# Format for plotting
df_res_tab$hr <- factor(df_res_tab$hr, levels = c("0.5", "1", "2"))
df_res_tab$n <- factor(df_res_tab$n, levels = c("1000", "2500", "5000", "10000", "25000", "50000"))

# Plot mean and sim SEs
pl_se <- ggplot(data = df_res_tab, aes(x = se_log_mean, y = se_log_sim, color = hr, shape = hr, size = n)) +
  geom_abline(intercept = 0, slope = 1, color = "gray30", lty = "dashed", linewidth = 0.5) +
  geom_point() +
  scale_color_manual(values = c("0.5" = "#6796A7",
                                "1" = "#EE756A",
                                "2" = "#824F73"),
                     guide = guide_legend(order = 1)) +
  scale_shape_manual(values = c("0.5" = 16,
                                "1" = 15,
                                "2" = 17),
                     guide = guide_legend(order = 1)) +
  scale_size_manual(values = c("1000" = 1.75,
                               "2500" = 2.5,
                               "5000" = 3.25,
                               "10000" = 4,
                               "25000" = 4.75,
                               "50000" = 5.5),
                    guide = guide_legend(order = 2)) +
  labs(title = NULL, x = "Mean model SE", y = "Simulation SE", color = "HR", shape = "HR", size = "Sample size") +
  theme_minimal(base_size = 13) +
  theme(panel.border =  element_rect(color = "black", fill = NA, linewidth = 0.5),
        legend.title = element_text(hjust = 0.5),
        legend.position = c(1.35, 0.55),
        plot.margin = margin(t = -10, r = 200, b = 5, l = 5))

# Save as pdf
ggsave("pl_se.pdf", plot = pl_se, width = 4.5, height = 2.25)

# Save results
df_res_se <- df_res
saveRDS(df_res_se, file = "df_res_se.rds")
save(df_res_se, file = "df_res_se.rda")



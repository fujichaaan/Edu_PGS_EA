# Meta-analysis for Edu PGS
library(tidyverse); library(data.table); library(nephro); library(gtsummary); library(rqlm); library(ggsci); library(patchwork); library(meta)

output_dir <- "./output/"


# PRS edu (R2) ------------------------------------------------------------

read_csv("r2.csv") |>
  ggplot(aes(x = group, y = r2, fill = gwas)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(x = group, y = r2, label = round(r2, 4)),
            position = position_dodge(0.9), vjust = 2, size = 6, show.legend = FALSE, color = "white") +
  labs(y = "Nagelkerke's R2") +
  scale_y_continuous(limits = c(0, 0.025)) +
  scale_fill_npg() +
  scale_color_npg() +
  theme_classic() +
  theme(axis.title = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 14),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.position = "top")

ggsave(paste0(output_dir, "Figure S4.png"), width = 16, height = 8, dpi = 600)



# PGS edu vs. lifestyle ---------------------------------------------------

oee_pgs_life <- read_csv("output/oee_02_edupgs_lifestyle.csv")
jap_pgs_life <- read_csv("output/jap_02_edupgs_lifestyle.csv")
total_pgs_life <- rbind(oee_pgs_life, jap_pgs_life) |>
  mutate(Array = c(rep("OEE", 8), rep("JAP", 8)))

traits <- c("smoking", "drinking", "less_exercise", "inbalanced_food")
meta_results <- list()

for (tr in traits) {
  for (md in c("Model1", "Model2")) {
    
    dt_tmp <- total_pgs_life |>
      filter(Trait == tr & Model == md)
    
    m.gen <- metagen(
      TE = Estimate,
      seTE = SE,
      data = dt_tmp,
      sm = "SMD",
      common = TRUE,
      random = FALSE,
      method.tau = "REML",
      method.random.ci = TRUE,
      title = paste("Trait:", tr, "|", md)
    )
    
    meta_results[[paste(tr, md, sep = "_")]] <- data.frame(
      Trait        = tr,
      Model        = md,
      TE_common    = round(m.gen$TE.common, 4),
      seTE_common  = round(m.gen$seTE.common, 4),
      lower_common = round(m.gen$lower.common, 4),
      upper_common = round(m.gen$upper.common, 4),
      all_common   = paste0(
        sprintf("%2.2f", m.gen$TE.common * 100), " (",
        sprintf("%2.2f", m.gen$lower.common * 100), ", ",
        sprintf("%2.2f", m.gen$upper.common * 100), ")"
      ),
      pval_common  = m.gen$pval.common,
      Q            = m.gen$Q,
      pval_Q       = m.gen$pval.Q
    )
  }
}

final_meta_table <- do.call(rbind, meta_results)
print(final_meta_table)

write_csv(final_meta_table, paste0(output_dir,"meta_02_edupgs_lifestyle.csv"))



## visualize ----

final_meta_table |>
  mutate(Trait = factor(Trait, levels = c("inbalanced_food", "less_exercise", "drinking", "smoking"))) |>
  ggplot(aes(x = TE_common, y = Trait, color = Model)) +
  geom_vline(xintercept = 0, color = "grey", linetype = "dashed") +
  geom_pointrange(aes(xmin = lower_common, xmax = upper_common), position = position_dodge(0.6) , size = 1.5) +
  # geom_text(aes(x = TE_common, y = Trait, label = all_common, color = Model), position = position_dodge(0.8),
  #          vjust = -1.8, size = 5) +
  scale_x_continuous(limits = c(-0.026, 0.026), 
                     breaks = c(-0.025, 0, 0.025),
                     labels = scales::percent) +
  scale_y_discrete(labels = c("Unbalanced food", "Physical Inactivity", "Drinking", "Smoking")) +
  labs(x = "PD (95% CI) \n per 1-SD Edu PGS") +
  scale_color_manual(values = c("#008080", "#4b0082")) +
  theme_classic() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 16, hjust = 0),
        legend.position = "none")

ggsave(paste0(output_dir, "Figure 3A.png"), width = 5, height = 8, dpi = 600)




# PGS edu vs. phenotype ---------------------------------------------------

oee_pgs_pheno <- read_csv("output/oee_03_edupgs_phenotype.csv")
jap_pgs_pheno <- read_csv("output/jap_03_edupgs_phenotype.csv")
total_pgs_pheno <- rbind(oee_pgs_pheno, jap_pgs_pheno) |>
  mutate(Array = c(rep("OEE", 8), rep("JAP", 8)))

traits <- c("obesity", "ht", "dyslipid", "dm")
meta_results <- list()


for (tr in traits) {
  for (md in c("Model1", "Model2")) {
    
    dt_tmp <- total_pgs_pheno |>
      filter(Trait == tr & Model == md)
    
    m.gen <- metagen(
      TE = Estimate,
      seTE = SE,
      data = dt_tmp,
      sm = "SMD",
      common = TRUE,
      random = FALSE,
      method.tau = "REML",
      method.random.ci = TRUE,
      title = paste("Trait:", tr, "|", md)
    )
    
    meta_results[[paste(tr, md, sep = "_")]] <- data.frame(
      Trait        = tr,
      Model        = md,
      TE_common    = round(m.gen$TE.common, 4),
      seTE_common  = round(m.gen$seTE.common, 4),
      lower_common = round(m.gen$lower.common, 4),
      upper_common = round(m.gen$upper.common, 4),
      all_common   = paste0(
        sprintf("%2.2f", m.gen$TE.common * 100), " (",
        sprintf("%2.2f", m.gen$lower.common * 100), ", ",
        sprintf("%2.2f", m.gen$upper.common * 100), ")"
      ),
      pval_common  = m.gen$pval.common,
      Q            = m.gen$Q,
      pval_Q       = m.gen$pval.Q
    )
  }
}

final_meta_table <- do.call(rbind, meta_results)
print(final_meta_table)

write_csv(final_meta_table, paste0(output_dir,"meta_03_edupgs_phenotype.csv"))



## visualize ----

final_meta_table |>
  mutate(Trait = factor(Trait, levels = c("dm", "dyslipid", "ht", "obesity"))) |>
  ggplot(aes(x = TE_common, y = Trait, color = Model)) +
  geom_vline(xintercept = 0, color = "grey", linetype = "dashed") +
  geom_pointrange(aes(xmin = lower_common, xmax = upper_common), position = position_dodge(0.6), size = 1.5) +
  # geom_text(aes(x = TE_common, y = Trait, label = all_common, color = Model), 
  # position = position_dodge(0.8), vjust = -1.8, size = 5) +
  scale_x_continuous(limits = c(-0.026, 0.026), 
                     breaks = c(-0.025, 0, 0.025),
                     labels = scales::percent) +
  scale_y_discrete(labels = c("T2DM", "Dyslipidemia", "Hypertension", "Obesity")) +
  labs(x = "PD (95% CI) \n per 1-SD Edu PGS") +
  scale_color_manual(values = c("#008080", "#4b0082")) +
  theme_classic() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 16, hjust = 0),
        legend.position = "none")

ggsave(paste0(output_dir, "Figure 4A.png"), width = 5, height = 8, dpi = 600)



# PGS edu vs. mortality ---------------------------------------------------

oee_pgs_mortality <- read_csv("output/oee_05_edupgs_poisson.csv")
jap_pgs_mortality <- read_csv("output/jap_05_edupgs_poisson.csv")
total_pgs_mortality <- rbind(oee_pgs_mortality, jap_pgs_mortality) |>
  mutate(Array = c(rep("OEE", 6), rep("JAP", 6)))

traits <- c("death_10yr", "death_cvd_10yr", "death_cancer_10yr")
meta_results <- list()

for (tr in traits) {
  
  for(md in c("Model1", "Model2")) {
    
    dt_trait <- total_pgs_mortality |> 
      filter(Trait == tr & Model == md)
    
    m.gen <- metagen(
      TE = log(Estimate),
      seTE = SE,
      data = dt_trait,
      sm = "RR",
      common = TRUE,
      random = FALSE,
      method.tau = "REML",
      method.random.ci = TRUE,
      title = paste("Trait:", tr)
    )
    
    meta_results[[paste(tr, md, sep = "_")]] <- data.frame(
      Trait        = tr,
      Model        = md,
      TE_common    = round(exp(m.gen$TE.common), 4),
      seTE_common  = round(m.gen$seTE.common, 4),
      lower_common = round(exp(m.gen$lower.common), 4),
      upper_common = round(exp(m.gen$upper.common), 4),
      all_common   = paste0(sprintf("%1.2f", exp(m.gen$TE.common)), " (",
                            sprintf("%1.2f", exp(m.gen$lower.common)), ", ",
                            sprintf("%1.2f", exp(m.gen$upper.common)), ")"),
      pval_common  = m.gen$pval.common,
      Q            = m.gen$Q,
      pval_Q       = m.gen$pval.Q
    )
  }
}

final_meta_table <- do.call(rbind, meta_results)
print(final_meta_table)

write_csv(final_meta_table, paste0(output_dir,"meta_05_edupgs_mortality.csv"))



## visualize ----

final_meta_table |>
  mutate(Trait = factor(Trait, levels = c("death_cancer_10yr", "death_cvd_10yr", "death_10yr"))) |>
  ggplot(aes(x = TE_common, y = Trait, color = Model)) +
  geom_vline(xintercept = 1, color = "grey", linetype = "dashed") +
  geom_pointrange(aes(xmin = lower_common, xmax = upper_common), position = position_dodge(0.6), size = 1.5) +
  # geom_text(aes(x = TE_common, y = Trait, label = all_common, color = Model), 
  #           position = position_dodge(0.8), vjust = -1.8, size = 5) +
  scale_x_continuous(limits = c(0.75, 1.15), 
                     breaks = c(0.8, 0.9, 1),
                     trans = "log") +
  scale_y_discrete(labels = c("Cancer", "CVD", "All-cause")) +
  labs(x = "RR (95% CI) \n per 1-SD Edu PGS") +
  scale_color_manual(values = c("#008080", "#4b0082")) +
  theme_classic() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 16, hjust = 0),
        legend.position = "none")

ggsave(paste0(output_dir, "Figure 5A.png"), width = 5, height = 8, dpi = 600)



# PGS edu + Edu vs. lifestyle ---------------------------------------------

oee_pgs_edu_life <- read_csv("output/oee_06_edupgs_edu_lifestyle.csv")
jap_pgs_edu_life <- read_csv("output/jap_06_edupgs_edu_lifestyle.csv")
total_pgs_edu_life <- rbind(oee_pgs_edu_life, jap_pgs_edu_life) |>
  mutate(Array = c(rep("OEE", 24), rep("JAP", 24)))


traits <- c("smoking", "drinking", "less_exercise", "inbalanced_food")
categories <- c("2nd", "3rd", "4th", "5th", "6th")
meta_results <- list()

for (tr in traits) {
  
  meta_results[[paste(tr, "1st", sep = "_")]] <- data.frame(
    Trait        = tr,
    Category     = "1st",
    TE_common    = 0,
    seTE_common  = NA,
    lower_common = 0,
    upper_common = 0,
    all_common   = "0 (Reference)",
    pval_common  = NA,
    Q            = NA,
    pval_Q       = NA
  )
  
  for(ca in categories){
    
    dt_trait <- total_pgs_edu_life |> 
      filter(Traits == tr & Category == ca)
    
    m.gen <- metagen(
      TE = Estimate,
      seTE = SE,
      data = dt_trait,
      sm = "SMD",
      common = TRUE,
      random = FALSE,
      method.tau = "REML",
      method.random.ci = TRUE,
      title = paste("Trait:", tr)
    )
    
    meta_results[[paste(tr, ca, sep = "_")]] <- data.frame(
      Trait        = tr,
      Category     = ca,
      TE_common    = round(m.gen$TE.common, 4),
      seTE_common  = round(m.gen$seTE.common, 4),
      lower_common = round(m.gen$lower.common, 4),
      upper_common = round(m.gen$upper.common, 4),
      all_common   = paste0(
        sprintf("%2.2f", m.gen$TE.common * 100), " (",
        sprintf("%2.2f", m.gen$lower.common * 100), ", ",
        sprintf("%2.2f", m.gen$upper.common * 100), ")"
      ),
      pval_common  = m.gen$pval.common,
      Q            = m.gen$Q,
      pval_Q       = m.gen$pval.Q
    )
  }
}

final_meta_table <- do.call(rbind, meta_results) |>
  mutate(Edu = factor(rep(c("Short", "Long"), 12),
                      levels = c("Short", "Long")),
         PGS = factor(rep(c("Low", "Low", "Middle", "Middle", "High", "High"), 4),
                      levels = c("Low", "Middle", "High")))
print(final_meta_table)

write_csv(final_meta_table, paste0(output_dir,"meta_06_edupgs_edu_lifestyle.csv"))


## visualize ----

Fig3_smoking <- final_meta_table |>
  filter(Trait == "smoking") |>
  ggplot(aes(x = PGS, y = TE_common, color = Edu)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_pointrange(aes(ymin = lower_common, ymax = upper_common), 
                  position = position_dodge(width = 0.5), size = 1.5) +
  annotate("text", x = -Inf, y = Inf, label = "Smoking",
           hjust = -.2, vjust = 2, size = 5) +
  scale_y_continuous(limits = c(-0.18, 0.18), 
                     breaks = c(-0.15, -0.10, -0.05, 0, 0.05, 0.10, 0.15),
                     labels = scales::percent) +
  labs(x = "Edu PGS", y = "PD (95% CI)") +
  scale_color_nejm() +
  theme_classic() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        legend.position = "none")

Fig3_drinking <- final_meta_table |>
  filter(Trait == "drinking") |>
  ggplot(aes(x = PGS, y = TE_common, color = Edu)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_pointrange(aes(ymin = lower_common, ymax = upper_common), 
                  position = position_dodge(width = 0.5), size = 1.5) +
  annotate("text", x = -Inf, y = Inf, label = "Drinking",
           hjust = -.2, vjust = 2, size = 5) +
  scale_y_continuous(limits = c(-0.18, 0.18), 
                     breaks = c(-0.15, -0.10, -0.05, 0, 0.05, 0.10, 0.15),
                     labels = scales::percent) +
  labs(x = "Edu PGS", y = "PD (95% CI)") +
  scale_color_nejm() +
  theme_classic() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        legend.position = "none")

Fig3_exercise <- final_meta_table |>
  filter(Trait == "less_exercise") |>
  ggplot(aes(x = PGS, y = TE_common, color = Edu)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_pointrange(aes(ymin = lower_common, ymax = upper_common), 
                  position = position_dodge(width = 0.5), size = 1.5) +
  annotate("text", x = -Inf, y = Inf, label = "Physical inactivity",
           hjust = -.2, vjust = 2, size = 5) +
  scale_y_continuous(limits = c(-0.18, 0.18), 
                     breaks = c(-0.15, -0.10, -0.05, 0, 0.05, 0.10, 0.15),
                     labels = scales::percent) +
  labs(x = "Edu PGS", y = "PD (95% CI)") +
  scale_color_nejm() +
  theme_classic() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        legend.position = "none")

Fig3_food <- final_meta_table |>
  filter(Trait == "inbalanced_food") |>
  ggplot(aes(x = PGS, y = TE_common, color = Edu)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_pointrange(aes(ymin = lower_common, ymax = upper_common), 
                  position = position_dodge(width = 0.5), size = 1.5) +
  annotate("text", x = -Inf, y = Inf, label = "Unbalanced food",
           hjust = -.2, vjust = 2, size = 5) +
  scale_y_continuous(limits = c(-0.18, 0.18), 
                     breaks = c(-0.15, -0.10, -0.05, 0, 0.05, 0.10, 0.15),
                     labels = scales::percent) +
  labs(x = "Edu PGS", y = "PD (95% CI)") +
  scale_color_nejm() +
  theme_classic() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        legend.position = "none")

Fig3B <- (Fig3_smoking | Fig3_drinking) / (Fig3_exercise | Fig3_food)
Fig3B

ggsave(paste0(output_dir, "Figure 3B.png"), width = 10, height = 8, dpi = 600)




# PGS edu + Edu vs. phenotype ---------------------------------------------

oee_pgs_edu_pheno <- read_csv("output/oee_07_edupgs_edu_phenotype.csv")
jap_pgs_edu_pheno <- read_csv("output/jap_07_edupgs_edu_phenotype.csv")
total_pgs_edu_pheno <- rbind(oee_pgs_edu_pheno, jap_pgs_edu_pheno) |>
  mutate(Array = c(rep("OEE", 24), rep("JAP", 24)))


traits <- c("obesity", "ht", "dyslipid", "dm")
categories <- c("2nd", "3rd", "4th", "5th", "6th")
meta_results <- list()

for (tr in traits) {
  
  meta_results[[paste(tr, "1st", sep = "_")]] <- data.frame(
    Trait        = tr,
    Category     = "1st",
    TE_common    = 0,
    seTE_common  = NA,
    lower_common = 0,
    upper_common = 0,
    all_common   = "0 (Reference)",
    pval_common  = NA,
    Q            = NA,
    pval_Q       = NA
  )
  
  for(ca in categories){
    
    dt_trait <- total_pgs_edu_pheno |> 
      filter(Traits == tr & Category == ca)
    
    m.gen <- metagen(
      TE = Estimate,
      seTE = SE,
      data = dt_trait,
      sm = "SMD",
      common = TRUE,
      random = FALSE,
      method.tau = "REML",
      method.random.ci = TRUE,
      title = paste("Trait:", tr)
    )
    
    meta_results[[paste(tr, ca, sep = "_")]] <- data.frame(
      Trait        = tr,
      Category     = ca,
      TE_common    = round(m.gen$TE.common, 4),
      seTE_common  = round(m.gen$seTE.common, 4),
      lower_common = round(m.gen$lower.common, 4),
      upper_common = round(m.gen$upper.common, 4),
      all_common   = paste0(
        sprintf("%2.2f", m.gen$TE.common * 100), " (",
        sprintf("%2.2f", m.gen$lower.common * 100), ", ",
        sprintf("%2.2f", m.gen$upper.common * 100), ")"
      ),
      pval_common  = m.gen$pval.common,
      Q            = m.gen$Q,
      pval_Q       = m.gen$pval.Q
    )
  }
}

final_meta_table <- do.call(rbind, meta_results) |>
  mutate(Edu = factor(rep(c("Short", "Long"), 12),
                      levels = c("Short", "Long")),
         PGS = factor(rep(c("Low", "Low", "Middle", "Middle", "High", "High"), 4),
                      levels = c("Low", "Middle", "High")))
print(final_meta_table)

write_csv(final_meta_table, paste0(output_dir,"meta_07_edupgs_edu_phenotype.csv"))



## visualize ----

Fig4_obesity <- final_meta_table |>
  filter(Trait == "obesity") |>
  ggplot(aes(x = PGS, y = TE_common, color = Edu)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_pointrange(aes(ymin = lower_common, ymax = upper_common), 
                  position = position_dodge(width = 0.5), size = 1.5) +
  annotate("text", x = -Inf, y = Inf, label = "Obesity",
           hjust = -.2, vjust = 2, size = 5) +
  scale_y_continuous(limits = c(-0.065, 0.065), 
                     breaks = c(-0.05, 0, 0.05),
                     labels = scales::percent) +
  labs(x = "Edu PGS", y = "PD (95% CI)") +
  scale_color_nejm() +
  theme_classic() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        legend.position = "none")

Fig4_ht <- final_meta_table |>
  filter(Trait == "ht") |>
  ggplot(aes(x = PGS, y = TE_common, color = Edu)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_pointrange(aes(ymin = lower_common, ymax = upper_common), 
                  position = position_dodge(width = 0.5), size = 1.5) +
  annotate("text", x = -Inf, y = Inf, label = "Hypertension",
           hjust = -.2, vjust = 2, size = 5) +
  scale_y_continuous(limits = c(-0.065, 0.065), 
                     breaks = c(-0.05, 0, 0.05),
                     labels = scales::percent) +
  labs(x = "Edu PGS", y = "PD (95% CI)") +
  scale_color_nejm() +
  theme_classic() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        legend.position = "none")

Fig4_dyslipid <- final_meta_table |>
  filter(Trait == "dyslipid") |>
  ggplot(aes(x = PGS, y = TE_common, color = Edu)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_pointrange(aes(ymin = lower_common, ymax = upper_common), 
                  position = position_dodge(width = 0.5), size = 1.5) +
  annotate("text", x = -Inf, y = Inf, label = "Dyslipidemia",
           hjust = -.2, vjust = 2, size = 5) +
  scale_y_continuous(limits = c(-0.065, 0.065), 
                     breaks = c(-0.05, 0, 0.05),
                     labels = scales::percent) +
  labs(x = "Edu PGS", y = "PD (95% CI)") +
  scale_color_nejm() +
  theme_classic() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        legend.position = "none")

Fig4_dm <- final_meta_table |>
  filter(Trait == "dm") |>
  ggplot(aes(x = PGS, y = TE_common, color = Edu)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_pointrange(aes(ymin = lower_common, ymax = upper_common), 
                  position = position_dodge(width = 0.5), size = 1.5) +
  annotate("text", x = -Inf, y = Inf, label = "T2DM",
           hjust = -.2, vjust = 2, size = 5) +
  scale_y_continuous(limits = c(-0.065, 0.065), 
                     breaks = c(-0.05, 0, 0.05),
                     labels = scales::percent) +
  labs(x = "Edu PGS", y = "PD (95% CI)") +
  scale_color_nejm() +
  theme_classic() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        legend.position = "none")

Fig4B <- (Fig4_obesity | Fig4_ht) / (Fig4_dyslipid | Fig4_dm)
Fig4B

ggsave(paste0(output_dir, "Figure 4B.png"), width = 10, height = 8, dpi = 600)




# PGS edu + Edu vs. mortality ---------------------------------------------


oee_pgs_edu_mortality <- read_csv("output/oee_08_edupgs_edu_mortality.csv")
jap_pgs_edu_mortality <- read_csv("output/jap_08_edupgs_edu_mortality.csv")
total_pgs_edu_mortality <- rbind(oee_pgs_edu_mortality, jap_pgs_edu_mortality) |>
  mutate(Array = c(rep("OEE", 18), rep("JAP", 18)))


traits <- c("death_10yr", "death_cvd_10yr", "death_cancer_10yr")
categories <- c("2nd", "3rd", "4th", "5th", "6th")
meta_results <- list()

for (tr in traits) {
  
  meta_results[[paste(tr, "1st", sep = "_")]] <- data.frame(
    Trait        = tr,
    Category     = "1st",
    TE_common    = 1,
    seTE_common  = NA,
    lower_common = 1,
    upper_common = 1,
    all_common   = "1 (Reference)",
    pval_common  = NA,
    Q            = NA,
    pval_Q       = NA
  )
  
  for(ca in categories){
    
    dt_trait <- total_pgs_edu_mortality |> 
      filter(Traits == tr & Category == ca)
    
    m.gen <- metagen(
      TE = log(Estimate),
      seTE = SE,
      data = dt_trait,
      sm = "RR",
      common = TRUE,
      random = FALSE,
      method.tau = "REML",
      method.random.ci = TRUE,
      title = paste("Trait:", tr)
    )
    
    meta_results[[paste(tr, ca, sep = "_")]] <- data.frame(
      Trait        = tr,
      Category     = ca,
      TE_common    = round(exp(m.gen$TE.common), 4),
      seTE_common  = round(m.gen$seTE.common, 4),
      lower_common = round(exp(m.gen$lower.common), 4),
      upper_common = round(exp(m.gen$upper.common), 4),
      all_common   = paste0(
        sprintf("%1.2f", exp(m.gen$TE.common)), " (",
        sprintf("%1.2f", exp(m.gen$lower.common)), ", ",
        sprintf("%1.2f", exp(m.gen$upper.common)), ")"
      ),
      pval_common  = m.gen$pval.common,
      Q            = m.gen$Q,
      pval_Q       = m.gen$pval.Q
    )
  }
}

final_meta_table <- do.call(rbind, meta_results) |>
  mutate(Edu = factor(rep(c("Short", "Long"), 9),
                      levels = c("Short", "Long")),
         PGS = factor(rep(c("Low", "Low", "Middle", "Middle", "High", "High"), 3),
                      levels = c("Low", "Middle", "High")))
print(final_meta_table)

write_csv(final_meta_table, paste0(output_dir,"meta_08_edupgs_edu_mortality.csv"))


## visualize ----

Fig5_all <- final_meta_table |>
  filter(Trait == "death_10yr") |>
  ggplot(aes(x = PGS, y = TE_common, color = Edu)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey") +
  geom_pointrange(aes(ymin = lower_common, ymax = upper_common), 
                  position = position_dodge(width = 0.5), size = 1.5) +
  annotate("text", x = -Inf, y = Inf, label = "All-cause",
           hjust = -.2, vjust = 2, size = 5) +
  scale_y_continuous(limits = c(0.25, 2.6), 
                     breaks = c(0.25, 0.5, 1, 2),
                     trans = "log") +
  labs(x = "Edu PGS", y = "RR (95% CI)") +
  scale_color_nejm() +
  theme_classic() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        legend.position = "none")

Fig5_cvd <- final_meta_table |>
  filter(Trait == "death_cvd_10yr") |>
  ggplot(aes(x = PGS, y = TE_common, color = Edu)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey") +
  geom_pointrange(aes(ymin = lower_common, ymax = upper_common), 
                  position = position_dodge(width = 0.5), size = 1.5) +
  annotate("text", x = -Inf, y = Inf, label = "CVD",
           hjust = -.2, vjust = 2, size = 5) +
  scale_y_continuous(limits = c(0.25, 2.6), 
                     breaks = c(0.25, 0.5, 1, 2),
                     trans = "log") +
  labs(x = "Edu PGS", y = "RR (95% CI)") +
  scale_color_nejm() +
  theme_classic() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        legend.position = "none")

Fig5_cancer <- final_meta_table |>
  filter(Trait == "death_cancer_10yr") |>
  ggplot(aes(x = PGS, y = TE_common, color = Edu)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey") +
  geom_pointrange(aes(ymin = lower_common, ymax = upper_common), 
                  position = position_dodge(width = 0.5), size = 1.5) +
  annotate("text", x = -Inf, y = Inf, label = "Cancer",
           hjust = -.2, vjust = 2, size = 5) +
  scale_y_continuous(limits = c(0.25, 2.6), 
                     breaks = c(0.25, 0.5, 1, 2),
                     trans = "log") +
  labs(x = "Edu PGS", y = "RR (95% CI)") +
  scale_color_nejm() +
  theme_classic() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        legend.position = "none")

Fig5B <- (Fig5_all | Fig5_cvd) / (Fig5_cancer | Fig5_cancer)
Fig5B

ggsave(paste0(output_dir, "Figure 5B.png"), width = 10, height = 8, dpi = 600)




# Edu vs. lifestyle w Edu PGS ---------------------------------------------


oee_pgs_life <- read_csv("output/oee_09_edu_lifestyle.csv")
jap_pgs_life <- read_csv("output/jap_09_edu_lifestyle.csv")
total_pgs_life <- rbind(oee_pgs_life, jap_pgs_life) |>
  mutate(Array = c(rep("OEE", 8), rep("JAP", 8)))

traits <- c("smoking", "drinking", "less_exercise", "inbalanced_food")
meta_results <- list()

for (tr in traits) {
  for (va in c("univ_comp", "prs_edu_binary")) {
    
    dt_tmp <- total_pgs_life |>
      filter(Traits == tr & Variable == va)
    
    m.gen <- metagen(
      TE = Estimate,
      seTE = SE,
      data = dt_tmp,
      sm = "SMD",
      common = TRUE,
      random = FALSE,
      method.tau = "REML",
      method.random.ci = TRUE,
      title = paste("Trait:", tr, "|", va)
    )
    
    meta_results[[paste(tr, va, sep = "_")]] <- data.frame(
      Trait        = tr,
      Variable        = va,
      TE_common    = round(m.gen$TE.common, 4),
      seTE_common  = round(m.gen$seTE.common, 4),
      lower_common = round(m.gen$lower.common, 4),
      upper_common = round(m.gen$upper.common, 4),
      all_common   = paste0(
        sprintf("%2.2f", m.gen$TE.common * 100), " (",
        sprintf("%2.2f", m.gen$lower.common * 100), ", ",
        sprintf("%2.2f", m.gen$upper.common * 100), ")"
      ),
      pval_common  = m.gen$pval.common,
      Q            = m.gen$Q,
      pval_Q       = m.gen$pval.Q
    )
  }
}

final_meta_table <- do.call(rbind, meta_results)
print(final_meta_table)

write_csv(final_meta_table, paste0(output_dir,"meta_09_edu_lifestyle.csv"))



## visualize ----

final_meta_table |>
  mutate(Trait = factor(Trait, levels = c("inbalanced_food", "less_exercise", "drinking", "smoking"))) |>
  ggplot(aes(x = TE_common, y = Trait, color = Variable)) +
  geom_vline(xintercept = 0, color = "grey", linetype = "dashed") +
  geom_pointrange(aes(xmin = lower_common, xmax = upper_common), , position = position_dodge(0.6) , size = 1.2) +
  #  geom_text(aes(x = TE_common, y = Trait, label = all_common), 
  #            color = "black", vjust = -1.8, size = 5) +
  scale_color_manual(values = c("#e5a323", "#65ab31")) +
  scale_x_continuous(limits = c(-0.15, 0.15), 
                     breaks = c(-0.15, -0.10, -0.05, 0, 0.05, 0.10, 0.15),
                     labels = scales::percent) +
  scale_y_discrete(labels = c("Unbalanced food", "Physical Inactivity", "Drinking", "Smoking")) +
  labs(x = "PD (95% CI)") +
  theme_classic() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 16, hjust = 0),
        legend.position = "none")

ggsave(paste0(output_dir, "Figure S1.png"), width = 6, height = 8, dpi = 600)



# Edu vs. phenotype w Edu PGS ---------------------------------------------


oee_pgs_pheno <- read_csv("output/oee_10_edu_phenotype.csv")
jap_pgs_pheno <- read_csv("output/jap_10_edu_phenotype.csv")
total_pgs_pheno <- rbind(oee_pgs_pheno, jap_pgs_pheno) |>
  mutate(Array = c(rep("OEE", 8), rep("JAP", 8)))

traits <- c("obesity", "ht", "dyslipid", "dm")
meta_results <- list()

for (tr in traits) {
  for (va in c("univ_comp", "prs_edu_binary")) {
    
    dt_tmp <- total_pgs_pheno |>
      filter(Traits == tr & Variable == va)
    
    m.gen <- metagen(
      TE = Estimate,
      seTE = SE,
      data = dt_tmp,
      sm = "SMD",
      common = TRUE,
      random = FALSE,
      method.tau = "REML",
      method.random.ci = TRUE,
      title = paste("Trait:", tr, "|", va)
    )
    
    meta_results[[paste(tr, va, sep = "_")]] <- data.frame(
      Trait        = tr,
      Variable        = va,
      TE_common    = round(m.gen$TE.common, 4),
      seTE_common  = round(m.gen$seTE.common, 4),
      lower_common = round(m.gen$lower.common, 4),
      upper_common = round(m.gen$upper.common, 4),
      all_common   = paste0(
        sprintf("%2.2f", m.gen$TE.common * 100), " (",
        sprintf("%2.2f", m.gen$lower.common * 100), ", ",
        sprintf("%2.2f", m.gen$upper.common * 100), ")"
      ),
      pval_common  = m.gen$pval.common,
      Q            = m.gen$Q,
      pval_Q       = m.gen$pval.Q
    )
  }
}

final_meta_table <- do.call(rbind, meta_results)
print(final_meta_table)

write_csv(final_meta_table, paste0(output_dir,"meta_10_edu_pheno.csv"))



## visualize ----

final_meta_table |>
  mutate(Trait = factor(Trait, levels = c("dm", "dyslipid", "ht", "obesity"))) |>
  ggplot(aes(x = TE_common, y = Trait, color = Variable)) +
  geom_vline(xintercept = 0, color = "grey", linetype = "dashed") +
  geom_pointrange(aes(xmin = lower_common, xmax = upper_common), , position = position_dodge(0.6) , size = 1.2) +
  #  geom_text(aes(x = TE_common, y = Trait, label = all_common), 
  #            color = "black", vjust = -1.8, size = 5) +
  scale_color_manual(values = c("#e5a323", "#65ab31")) +
  scale_x_continuous(limits = c(-0.05, 0.05), 
                     breaks = c(-0.05, -0.025, 0, 0.025, 0.05),
                     labels = scales::percent) +
  scale_y_discrete(labels = c("T2DM", "Dyslipidemia", "Hypertension", "Obesity")) +
  labs(x = "PD (95% CI)") +
  theme_classic() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 16, hjust = 0),
        legend.position = "none")

ggsave(paste0(output_dir, "Figure S2.png"), width = 6, height = 8, dpi = 600)



# Edu vs. mortality w Edu PGS ---------------------------------------------


oee_pgs_mortality <- read_csv("output/oee_11_edu_poisson.csv")
jap_pgs_mortality <- read_csv("output/jap_11_edu_poisson.csv")
total_pgs_mortality <- rbind(oee_pgs_mortality, jap_pgs_mortality) |>
  mutate(Array = c(rep("OEE", 6), rep("JAP", 6)))

traits <- c("death_10yr", "death_cvd_10yr", "death_cancer_10yr")
meta_results <- list()

for (tr in traits) {
  for (va in c("univ_comp", "prs_edu_binary")) {
    
    dt_trait <- total_pgs_mortality |> 
      filter(Traits == tr & Variable == va)
    
    m.gen <- metagen(
      TE = log(Estimate),
      seTE = SE,
      data = dt_trait,
      sm = "RR",
      common = TRUE,
      random = FALSE,
      method.tau = "REML",
      method.random.ci = TRUE,
      title = paste("Trait:", tr, "|", va)
    )
    
    meta_results[[paste(tr, va, sep = "_")]] <- data.frame(
      Trait        = tr,
      Variable        = va,
      TE_common    = round(exp(m.gen$TE.common), 4),
      seTE_common  = round(m.gen$seTE.common, 4),
      lower_common = round(exp(m.gen$lower.common), 4),
      upper_common = round(exp(m.gen$upper.common), 4),
      all_common   = paste0(sprintf("%1.2f", exp(m.gen$TE.common)), " (",
                            sprintf("%1.2f", exp(m.gen$lower.common)), ", ",
                            sprintf("%1.2f", exp(m.gen$upper.common)), ")"),
      pval_common  = m.gen$pval.common,
      Q            = m.gen$Q,
      pval_Q       = m.gen$pval.Q
    )
  }
}

final_meta_table <- do.call(rbind, meta_results)
print(final_meta_table)

write_csv(final_meta_table, paste0(output_dir,"meta_11_edu_poisson.csv"))



## visualize ----

final_meta_table |>
  mutate(Trait = factor(Trait, levels = c("death_cancer_10yr", "death_cvd_10yr", "death_10yr"))) |>
  ggplot(aes(x = TE_common, y = Trait, color = Variable)) +
  geom_vline(xintercept = 1, color = "grey", linetype = "dashed") +
  geom_pointrange(aes(xmin = lower_common, xmax = upper_common), position = position_dodge(0.6), size = 1.2) +
  #  geom_text(aes(x = TE_common, y = Trait, label = all_common), 
  #            color = "black", vjust = -1.8, size = 5) +
  scale_color_manual(values = c("#e5a323", "#65ab31")) +
  scale_x_continuous(limits = c(0.4, 1.5), 
                     breaks = c(0.5, 0.75, 1, 1.5),
                     trans = "log") +
  scale_y_discrete(labels = c("Cancer", "CVD", "All-cause")) +
  labs(x = "RR (95% CI)") +
  theme_classic() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 16, hjust = 0),
        legend.position = "none")

ggsave(paste0(output_dir, "Figure S3.png"), width = 6, height = 8, dpi = 600)




# Interaction (PGS Edu + Edu) for lifestlye ----

oee_pgs_lifestyle_int <- read_csv("output/oee_12_edupgs_interaction.csv")
jap_pgs_lifestyle_int <- read_csv("output/jap_12_edupgs_interaction.csv")
total_pgs_lifestyle_int <- rbind(oee_pgs_lifestyle_int, jap_pgs_lifestyle_int) |>
  mutate(Array = c(rep("OEE", 4), rep("JAP", 4)))

traits <- c("smoking", "drinking", "less_exercise", "inbalanced_food")
meta_results <- list()

for (tr in traits) {
  
  dt_tmp <- total_pgs_lifestyle_int |>
    filter(Trait == tr)
  
  m.gen <- metagen(
    TE = Estimate,
    seTE = SE,
    data = dt_tmp,
    sm = "MD",
    common = TRUE,
    random = FALSE,
    method.tau = "REML",
    method.random.ci = TRUE,
    title = paste("Interaction | Trait:", tr)
  )
  
  meta_results[[tr]] <- data.frame(
    Trait        = tr,
    TE_common    = round(m.gen$TE.common, 4),
    seTE_common  = round(m.gen$seTE.common, 4),
    lower_common = round(m.gen$lower.common, 4),
    upper_common = round(m.gen$upper.common, 4),
    all_common   = paste0(
      sprintf("%2.2f", m.gen$TE.common * 100), " (",
      sprintf("%2.2f", m.gen$lower.common * 100), ", ",
      sprintf("%2.2f", m.gen$upper.common * 100), ")"
    ),
    pval_common  = m.gen$pval.common,
    Q            = m.gen$Q,
    pval_Q       = m.gen$pval.Q
  )
}

final_meta_interaction <- do.call(rbind, meta_results)
print(final_meta_interaction)

write_csv(
  final_meta_interaction,
  paste0(output_dir, "meta_12_edupgs_edu_lifestyle_int.csv")
)



# Interaction (PGS Edu + Edu) for phenotype ----

oee_pgs_pheno_int <- read_csv("output/oee_13_edupgs_interaction.csv")
jap_pgs_pheno_int <- read_csv("output/jap_13_edupgs_interaction.csv")
total_pgs_pheno_int <- rbind(oee_pgs_pheno_int, jap_pgs_pheno_int) |>
  mutate(Array = c(rep("OEE", 4), rep("JAP", 4)))

traits <- c("obesity", "ht", "dyslipid", "dm")
meta_results <- list()

for (tr in traits) {
  
  dt_tmp <- total_pgs_pheno_int |>
    filter(Trait == tr)
  
  m.gen <- metagen(
    TE = Estimate,
    seTE = SE,
    data = dt_tmp,
    sm = "MD",
    common = TRUE,
    random = FALSE,
    method.tau = "REML",
    method.random.ci = TRUE,
    title = paste("Interaction | Trait:", tr)
  )
  
  meta_results[[tr]] <- data.frame(
    Trait        = tr,
    TE_common    = round(m.gen$TE.common, 4),
    seTE_common  = round(m.gen$seTE.common, 4),
    lower_common = round(m.gen$lower.common, 4),
    upper_common = round(m.gen$upper.common, 4),
    all_common   = paste0(
      sprintf("%2.2f", m.gen$TE.common * 100), " (",
      sprintf("%2.2f", m.gen$lower.common * 100), ", ",
      sprintf("%2.2f", m.gen$upper.common * 100), ")"
    ),
    pval_common  = m.gen$pval.common,
    Q            = m.gen$Q,
    pval_Q       = m.gen$pval.Q
  )
}

final_meta_interaction <- do.call(rbind, meta_results)
print(final_meta_interaction)

write_csv(
  final_meta_interaction,
  paste0(output_dir, "meta_13_edupgs_edu_phenotype_int.csv")
)




# Interaction (PGS Edu + Edu) for mortality ----

oee_pgs_mortality_int <- read_csv("output/oee_14_edupgs_interaction.csv")
jap_pgs_mortality_int <- read_csv("output/jap_14_edupgs_interaction.csv")
total_pgs_mortality_int <- rbind(oee_pgs_mortality_int, jap_pgs_mortality_int) |>
  mutate(Array = c(rep("OEE", 3), rep("JAP", 3)))

traits <- c("death_10yr", "death_cvd_10yr", "death_cancer_10yr")
meta_results <- list()

for (tr in traits) {
  
  dt_trait <- total_pgs_mortality_int |> 
    filter(Trait == tr)
  
  m.gen <- metagen(
    TE = log(Estimate),
    seTE = SE,
    data = dt_trait,
    sm = "RR",
    common = TRUE,
    random = FALSE,
    method.tau = "REML",
    method.random.ci = TRUE,
    title = paste("Trait:", tr)
  )
  
  meta_results[[paste(tr, md, sep = "_")]] <- data.frame(
    Trait        = tr,
    Model        = md,
    TE_common    = round(exp(m.gen$TE.common), 4),
    seTE_common  = round(m.gen$seTE.common, 4),
    lower_common = round(exp(m.gen$lower.common), 4),
    upper_common = round(exp(m.gen$upper.common), 4),
    all_common   = paste0(sprintf("%1.2f", exp(m.gen$TE.common)), " (",
                          sprintf("%1.2f", exp(m.gen$lower.common)), ", ",
                          sprintf("%1.2f", exp(m.gen$upper.common)), ")"),
    pval_common  = m.gen$pval.common,
    Q            = m.gen$Q,
    pval_Q       = m.gen$pval.Q
  )
}

final_meta_table <- do.call(rbind, meta_results)
print(final_meta_table)

write_csv(final_meta_table, paste0(output_dir,"meta_14_edupgs_edu_mortality_int.csv"))
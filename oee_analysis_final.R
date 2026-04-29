# Edu PGS analysis in J-MICC (OEE)
library(tidyverse); library(data.table); library(nephro); library(gtsummary); library(rqlm); library(ggsci); library(patchwork); library(performance)

output_dir <- "./output/"
load("dataset_oee_imp.RData")


# Table 1 -----------------------------------------------------------------

dt |>
  group_by(edu) |>
  summarise(count = n(),
            per = count / nrow(dt) * 100)

Table1 <- dt |>
  select(prs_edu3, adi2010_2, site, age, male, sbp, dbp, bmi, ldl, tg, fbs, ht, dyslipid, dm, smoking, drinking, univ_comp, less_exercise, inbalanced_food) |>
  tbl_summary(by = prs_edu3,
              label = list(male         ~ "Men",
                           age          ~ "Age",
                           sbp          ~ "SBP",
                           dbp          ~ "DBP",
                           bmi          ~ "BMI",
                           ldl          ~ "LDL-C",
                           tg           ~ "TG",
                           fbs         ~ "Fasting BS",
                           ht   ~ "Hypertension",
                           dyslipid      ~ "Dyslipidemia",
                           dm      ~ "Type 2 Diabetes",
                           smoking   ~ "Current/Ever smoker",
                           drinking  ~ "Curret/Ever drinker",
                           less_exercise  ~ "Non Habitual exercise",
                           univ_comp       ~ "University Completion"),
              digits = list(age      ~ c(0, 1),
                            sbp      ~ c(1, 1),
                            dbp      ~ c(1, 1)),
              statistic = list(age   ~ c("{mean} ({sd})"),
                               sbp   ~ c("{mean} ({sd})"),
                               dbp   ~ c("{mean} ({sd})"),
                               bmi   ~ c("{mean} ({sd})")))

gtsummary::as_hux_xlsx(Table1, file = paste0(output_dir, "oee_00_table1.xlsx"))


all <- dt |>
  group_by(prs_edu10) |>
  summarise(univ_comp = mean(univ_comp, na.rm = TRUE))

sex <- dt |>
  group_by(prs_edu10, sex) |>
  summarise(univ_comp = mean(univ_comp, na.rm = TRUE)) |>
  pivot_wider(names_from = sex, values_from = univ_comp) |>
  rename(Men = `1`,
         Women = `2`)

all |>
  cbind(sex) |>
  write_csv(paste0(output_dir,"oee_01_edupgs_univ.csv"))


## Education -----

dt |>
  group_by(prs_edu10) |>
  summarise(univ_comp = mean(univ_comp, na.rm = TRUE)) |>
  ggplot(aes(x = prs_edu10, y = univ_comp)) +
  geom_point(color = "#BC3C29FF", size = 8) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Edu PGS decile", y = "Mean % of university completion") +
  theme_classic() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 16))


ggplot(dt, aes(y = fct_rev(prs_edu10), fill = fct_rev(edu_int))) +
  geom_bar(position = "fill") +
  scale_fill_manual(breaks = c("ES & JHS", "High School", "Junior College & Tech School", "Univ & Grad School"),
                    values = c("#BC3C29FF", "#0072B5FF", "#E18727FF", "#7876B1FF"),
                    labels = c("ES & JHS", "High School", "Junior College & Tech School", "College")) +
  scale_y_discrete(labels = c("10th", "9th", "8th", "7th", "6th",
                              "5th", "4th", "3rd", "2nd", "1st")) +
  scale_x_continuous(labels = scales::percent) +
  labs(x = "Percentage", y = "Edu PGS Decile", fill = "Educational attainment") +
  theme_minimal() +
  theme(legend.position = "top",
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12))

ggsave(paste0(output_dir, "oee_fig1A.png"), width = 10, height = 8, dpi = 600)


## Other data viz ----

library(ggdist)
ggplot(dt, aes(x = prs_edu, y = edu_int, color = edu_int, fill = edu_int)) +
  stat_dotsinterval(alpha = 0.8) +
  scale_fill_manual(breaks = c("ES & JHS", "High School", "Junior College & Tech School", "Univ & Grad School"),
                    values = c("#BC3C29FF", "#0072B5FF", "#E18727FF", "#7876B1FF"),
                    labels = c("ES & JHS", "High School", "Junior College & Tech School", "College")) +
  scale_color_manual(breaks = c("ES & JHS", "High School", "Junior College & Tech School", "Univ & Grad School"),
                     values = c("#BC3C29FF", "#0072B5FF", "#E18727FF", "#7876B1FF"),
                     labels = c("ES & JHS", "High School", "Junior College & Tech School", "College")) +
  scale_y_discrete(labels = c("ES & JHS", "High School", "Junior College & Tech School", "College")) +
  labs(x = "Edu PGS", y = "Educational attainement") +
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4)) +
  theme_minimal() +
  theme(axis.title = element_text(size = 18),
        axis.text.y = element_text(size = 14, hjust = 0),
        axis.text.x = element_text(size = 14, hjust = 0),
        legend.position = "none")

ggsave(paste0(output_dir, "oee_fig1B.png"), width = 9, height = 6, dpi = 600)


model <- glm(univ_comp ~ prs_edu, data = dt, family = binomial)
r2_nagelkerke(model)



# PGS edu vs. lifestyles -----------------------------------------

# Lifestyles
traits <- c("smoking", "drinking", "less_exercise", "inbalanced_food")
results <- list()

for (trait in traits) {
  
  dt_tmp <- dt[!is.na(dt[[trait]]), ]
  
  # ----- Model 1: prs only -----
  formula1 <- as.formula(
    paste(trait, "~ prs_edu + sex + age + I(age^2) + adi2010_2 +
          factor(birth_year) + V3 + V4 + V5 + V6 + V7")
  )
  model1 <- rqlm(formula1, data = dt_tmp, family = "gaussian", digits = 30)
  
  # ----- Model 2: prs + univ -----
  formula2 <- as.formula(
    paste(trait, "~ prs_edu + univ_comp + sex + age + I(age^2) + adi2010_2 +
          factor(birth_year) + V3 + V4 + V5 + V6 + V7")
  )
  model2 <- rqlm(formula2, data = dt_tmp, family = "gaussian", digits = 30)
  
  N <- nrow(dt_tmp)
  N.case <- sum(dt_tmp[[trait]] == 1, na.rm = TRUE)
  
  tmp_res <- rbind(
    
    data.frame(
      Trait = trait,
      Model = "Model1",
      Estimate = model1["prs_edu", "coef"],
      SE       = model1["prs_edu", "SE"],
      Lower    = model1["prs_edu", "CL"],
      Upper    = model1["prs_edu", "CU"],
      P.value  = model1["prs_edu", "P-value"],
      N_used   = N,
      N_case   = N.case,
      all      = paste0(
        sprintf("%2.2f", model1["prs_edu", "coef"] * 100), " (",
        sprintf("%2.2f", model1["prs_edu", "CL"] * 100), ", ",
        sprintf("%2.2f", model1["prs_edu", "CU"] * 100), ")"
      )
    ),
    
    data.frame(
      Trait = trait,
      Model = "Model2",
      Estimate = model2["prs_edu", "coef"],
      SE       = model2["prs_edu", "SE"],
      Lower    = model2["prs_edu", "CL"],
      Upper    = model2["prs_edu", "CU"],
      P.value  = model2["prs_edu", "P-value"],
      N_used   = N,
      N_case   = N.case,
      all      = paste0(
        sprintf("%2.2f", model2["prs_edu", "coef"] * 100), " (",
        sprintf("%2.2f", model2["prs_edu", "CL"] * 100), ", ",
        sprintf("%2.2f", model2["prs_edu", "CU"] * 100), ")"
      )
    )
  )
  
  results[[trait]] <- tmp_res
}

final_results <- do.call(rbind, results)
print(final_results)

write_csv(final_results, paste0(output_dir,"oee_02_edupgs_lifestyle.csv"))




# ----- Interaction model -----

interaction_results <- list()

for (trait in traits) {
  
  dt_tmp <- dt[!is.na(dt[[trait]]), ]
  
  formula_int <- as.formula(
    paste(trait, "~ prs_edu * univ_comp + sex + age + I(age^2) + adi2010_2 +
          factor(birth_year) + V3 + V4 + V5 + V6 + V7")
  )
  
  model_int <- rqlm(formula_int, data = dt_tmp, family = "gaussian", digits = 30)
  
  term_name <- "prs_edu:univ_comp"
  
  N <- nrow(dt_tmp)
  N.case <- sum(dt_tmp[[trait]] == 1, na.rm = TRUE)
  
  tmp_int <- data.frame(
    Trait = trait,
    Term  = term_name,
    Estimate = model_int[term_name, "coef"],
    SE       = model_int[term_name, "SE"],
    Lower    = model_int[term_name, "CL"],
    Upper    = model_int[term_name, "CU"],
    P.value  = model_int[term_name, "P-value"],
    N_used   = N,
    N_case   = N.case,
    all      = paste0(
      sprintf("%2.2f", model_int[term_name, "coef"] * 100), " (",
      sprintf("%2.2f", model_int[term_name, "CL"] * 100), ", ",
      sprintf("%2.2f", model_int[term_name, "CU"] * 100), ")"
    )
  )
  
  interaction_results[[trait]] <- tmp_int
}

interaction_results_df <- do.call(rbind, interaction_results)

print(interaction_results_df)

write_csv(
  interaction_results_df,
  paste0(output_dir, "oee_12_edupgs_interaction.csv")
)



# PGS edu vs. phenotype ---------------------------------------------------


traits <- c("obesity", "ht", "dyslipid", "dm")
results <- list()

for (trait in traits) {
  
  dt_tmp <- dt[!is.na(dt[[trait]]), ]
  
  # ----- Model 1: prs only -----
  formula1 <- as.formula(
    paste(trait, "~ prs_edu + sex + age + I(age^2) + adi2010_2 +
          factor(birth_year) + V3 + V4 + V5 + V6 + V7")
  )
  model1 <- rqlm(formula1, data = dt_tmp, family = "gaussian")
  
  # ----- Model 2: prs + univ -----
  formula2 <- as.formula(
    paste(trait, "~ prs_edu + univ_comp + sex + age + I(age^2) + adi2010_2 +
          factor(birth_year) + V3 + V4 + V5 + V6 + V7")
  )
  model2 <- rqlm(formula2, data = dt_tmp, family = "gaussian")
  
  N <- nrow(dt_tmp)
  N.case <- sum(dt_tmp[[trait]] == 1, na.rm = TRUE)
  
  tmp_res <- rbind(
    
    data.frame(
      Trait = trait,
      Model = "Model1",
      Estimate = model1["prs_edu", "coef"],
      SE       = model1["prs_edu", "SE"], 
      Lower    = model1["prs_edu", "CL"],
      Upper    = model1["prs_edu", "CU"],
      P.value  = model1["prs_edu", "P-value"],
      N_used   = N,
      N_case   = N.case,
      all      = paste0(
        sprintf("%2.2f", model1["prs_edu", "coef"] * 100), " (",
        sprintf("%2.2f", model1["prs_edu", "CL"] * 100), ", ",
        sprintf("%2.2f", model1["prs_edu", "CU"] * 100), ")"
      )
    ),
    
    data.frame(
      Trait = trait,
      Model = "Model2",
      Estimate = model2["prs_edu", "coef"],
      SE       = model2["prs_edu", "SE"], 
      Lower    = model2["prs_edu", "CL"],
      Upper    = model2["prs_edu", "CU"],
      P.value  = model2["prs_edu", "P-value"],
      N_used   = N,
      N_case   = N.case,
      all      = paste0(
        sprintf("%2.2f", model2["prs_edu", "coef"] * 100), " (",
        sprintf("%2.2f", model2["prs_edu", "CL"] * 100), ", ",
        sprintf("%2.2f", model2["prs_edu", "CU"] * 100), ")"
      )
    )
  )
  
  results[[trait]] <- tmp_res
}

final_results <- do.call(rbind, results)
print(final_results)

write_csv(final_results, paste0(output_dir,"oee_03_edupgs_phenotype.csv"))




# ----- Interaction model -----

interaction_results <- list()

for (trait in traits) {
  
  dt_tmp <- dt[!is.na(dt[[trait]]), ]
  
  formula_int <- as.formula(
    paste(trait, "~ prs_edu * univ_comp + sex + age + I(age^2) + adi2010_2 +
          factor(birth_year) + V3 + V4 + V5 + V6 + V7")
  )
  
  model_int <- rqlm(formula_int, data = dt_tmp, family = "gaussian", digits = 30)
  
  term_name <- "prs_edu:univ_comp"
  
  N <- nrow(dt_tmp)
  N.case <- sum(dt_tmp[[trait]] == 1, na.rm = TRUE)
  
  tmp_int <- data.frame(
    Trait = trait,
    Term  = term_name,
    Estimate = model_int[term_name, "coef"],
    SE       = model_int[term_name, "SE"],
    Lower    = model_int[term_name, "CL"],
    Upper    = model_int[term_name, "CU"],
    P.value  = model_int[term_name, "P-value"],
    N_used   = N,
    N_case   = N.case,
    all      = paste0(
      sprintf("%2.2f", model_int[term_name, "coef"] * 100), " (",
      sprintf("%2.2f", model_int[term_name, "CL"] * 100), ", ",
      sprintf("%2.2f", model_int[term_name, "CU"] * 100), ")"
    )
  )
  
  interaction_results[[trait]] <- tmp_int
}

interaction_results_df <- do.call(rbind, interaction_results)

print(interaction_results_df)

write_csv(
  interaction_results_df,
  paste0(output_dir, "oee_14_edupgs_interaction.csv")
)




# PGS edu vs. Mortality ---------------------------------------------------

table <- dt |>
  group_by(prs_edu3) |>
  summarise(case1 = sum(death_10yr),
            case2 = sum(death_cvd_10yr),
            case3 = sum(death_cancer_10yr),
            total = n())
table

write.csv(table, paste0(output_dir, "oee_04_edupgs_mortality.csv"))


# Lifestyles
traits <- c("death_10yr", "death_cvd_10yr", "death_cancer_10yr")
results <- list()

for (trait in traits) {
  
  dt_tmp <- dt[!is.na(dt[[trait]]), ]
  
  # ----- Model 1: prs only -----
  formula1 <- as.formula(
    paste(trait, "~ prs_edu + sex + age + I(age^2) + adi2010_2 +
          factor(birth_year) + V3 + V4 + V5 + V6 + V7")
  )
  model1 <- rqlm(formula1, data = dt_tmp, family = "poisson", eform = TRUE)
  
  # ----- Model 2: prs + univ -----
  formula2 <- as.formula(
    paste(trait, "~ prs_edu + univ_comp + sex + age + I(age^2) + adi2010_2 +
          factor(birth_year) + V3 + V4 + V5 + V6 + V7")
  )
  model2 <- rqlm(formula2, data = dt_tmp, family = "poisson", eform = TRUE)
  
  N <- nrow(dt_tmp)
  N.case <- sum(dt_tmp[[trait]] == 1, na.rm = TRUE)
  
  tmp_res <- rbind(
    
    data.frame(
      Trait = trait,
      Model = "Model1",
      Estimate = model1["prs_edu", "exp(coef)"],
      SE       = model1["prs_edu", "SE"], 
      Lower    = model1["prs_edu", "CL"],
      Upper    = model1["prs_edu", "CU"],
      P.value  = model1["prs_edu", "P-value"],
      N_used   = N,
      N_case   = N.case,
      all      = paste0(
        sprintf("%2.2f", model1["prs_edu", "exp(coef)"]), " (",
        sprintf("%2.2f", model1["prs_edu", "CL"]), ", ",
        sprintf("%2.2f", model1["prs_edu", "CU"]), ")"
      )
    ),
    
    data.frame(
      Trait = trait,
      Model = "Model2",
      Estimate = model2["prs_edu", "exp(coef)"],
      SE       = model2["prs_edu", "SE"], 
      Lower    = model2["prs_edu", "CL"],
      Upper    = model2["prs_edu", "CU"],
      P.value  = model2["prs_edu", "P-value"],
      N_used   = N,
      N_case   = N.case,
      all      = paste0(
        sprintf("%2.2f", model2["prs_edu", "exp(coef)"]), " (",
        sprintf("%2.2f", model2["prs_edu", "CL"]), ", ",
        sprintf("%2.2f", model2["prs_edu", "CU"]), ")"
      )
    )
  )
  
  results[[trait]] <- tmp_res
}

final_results <- do.call(rbind, results)
print(final_results)

write.csv(final_results, paste0(output_dir, "oee_05_edupgs_poisson.csv"))




# Interaction model
interaction_results <- list()

for (trait in traits) {
  
  dt_tmp <- dt[!is.na(dt[[trait]]), ]
  
  # ----- Interaction model -----
  formula_int <- as.formula(
    paste(trait, "~ prs_edu * univ_comp + sex + age + I(age^2) + adi2010_2 +
          factor(birth_year) + V3 + V4 + V5 + V6 + V7")
  )
  
  model_int <- rqlm(formula_int, data = dt_tmp, family = "poisson", eform = TRUE)
  
  term_name <- "prs_edu:univ_comp"
  
  if (!(term_name %in% rownames(model_int))) {
    term_name <- grep("prs_edu:univ_comp", rownames(model_int), value = TRUE)[1]
  }
  
  N <- nrow(dt_tmp)
  N.case <- sum(dt_tmp[[trait]] == 1, na.rm = TRUE)
  
  tmp_int <- data.frame(
    Trait = trait,
    Term  = term_name,
    Estimate = model_int[term_name, "exp(coef)"],
    SE       = model_int[term_name, "SE"],
    Lower    = model_int[term_name, "CL"],
    Upper    = model_int[term_name, "CU"],
    P.value  = model_int[term_name, "P-value"],
    N_used   = N,
    N_case   = N.case,
    all      = paste0(
      sprintf("%2.2f", model_int[term_name, "exp(coef)"]), " (",
      sprintf("%2.2f", model_int[term_name, "CL"]), ", ",
      sprintf("%2.2f", model_int[term_name, "CU"]), ")"
    )
  )
  
  interaction_results[[trait]] <- tmp_int
}

interaction_results_df <- do.call(rbind, interaction_results)

print(interaction_results_df)

write.csv(
  interaction_results_df,
  paste0(output_dir, "oee_15_edupgs_interaction.csv"),
  row.names = FALSE
)



# PGS edu + Edu vs. lifestyle ---------------------------------------------

traits <- c("smoking", "drinking", "less_exercise", "inbalanced_food")
results <- list()

for (trait in traits) {
  
  dt_tmp <- dt[!is.na(dt[[trait]]), ]
  
  formula <- as.formula(paste(trait, "~ prs_edu_edu + sex + age + I(age^2) + adi2010_2 + factor(birth_year) + V3 + V4 + V5 + V6 + V7"))
  
  model <- rqlm(formula, data = dt_tmp, family = "gaussian", digits = 30)
  
  coef_est <- c(0.00, model[2:6, "coef"])
  p_value <- c(NA, model[2:6, "P-value"])
  se <- c(0.00, model[2:6, "SE"])
  CL <- c(0.00, model[2:6, "CL"])
  CU <- c(0.00, model[2:6, "CU"])
  
  N_cat <- table(dt_tmp$prs_edu_edu)
  N_cat_vec <- as.numeric(N_cat)
  
  N_case_cat <- tapply(
    dt_tmp[[trait]] == 1,
    dt_tmp$prs_edu_edu,
    sum,
    na.rm = TRUE
  )
  N_case_vec <- as.numeric(N_case_cat)
  
  results[[trait]] <- data.frame(
    Traits = trait,
    PGS = rep(c("Low", "Middle", "High"), each = 2),
    Edu = rep(c("Short", "Long"), 3),
    Category = c("1st", "2nd", "3rd", "4th", "5th", "6th"),
    Estimate = coef_est,
    SE = se,
    CI.lower = CL,
    CI.upper = CU,
    N.use    = N_cat_vec,
    N.case   = N_case_vec,
    N = paste0(N_case_vec, "/", N_cat_vec),
    all      = paste0(
      sprintf("%2.2f", coef_est * 100), " (",
      sprintf("%2.2f", CL * 100), ", ",
      sprintf("%2.2f", CU * 100), ")"
    ),
    P.value = p_value
  )
}

final_results <- do.call(rbind, results)
print(final_results)

write_csv(final_results, paste0(output_dir,"oee_06_edupgs_edu_lifestyle.csv"))




# PGS edu + Edu vs. phenotype ---------------------------------------------------

traits <- c("obesity", "ht", "dyslipid", "dm")
results <- list()

for (trait in traits) {
  
  dt_tmp <- dt[!is.na(dt[[trait]]), ]
  
  formula <- as.formula(paste(trait, "~ prs_edu_edu + sex + age + I(age^2) + adi2010_2 + factor(birth_year) + V3 + V4 + V5 + V6 + V7"))
  
  model <- rqlm(formula, data = dt_tmp, family = "gaussian")
  
  coef_est <- c(0.00, model[2:6, "coef"])　
  p_value <- c(NA, model[2:6, "P-value"])
  se <- c(0.00, model[2:6, "SE"])
  CL <- c(0.00, model[2:6, "CL"])
  CU <- c(0.00, model[2:6, "CU"])
  
  N_cat <- table(dt_tmp$prs_edu_edu)
  N_cat_vec <- as.numeric(N_cat)
  
  N_case_cat <- tapply(
    dt_tmp[[trait]] == 1,
    dt_tmp$prs_edu_edu,
    sum,
    na.rm = TRUE
  )
  N_case_vec <- as.numeric(N_case_cat)
  
  results[[trait]] <- data.frame(
    Traits = trait,
    PGS = rep(c("Low", "Middle", "High"), each = 2),
    Edu = rep(c("Short", "Long"), 3),
    Category = c("1st", "2nd", "3rd", "4th", "5th", "6th"),
    Estimate = coef_est,
    SE = se,
    CI.lower = CL,
    CI.upper = CU,
    N.use    = N_cat_vec,
    N.case   = N_case_vec,
    N = paste0(N_case_vec, "/", N_cat_vec),
    all      = paste0(
      sprintf("%2.2f", coef_est * 100), " (",
      sprintf("%2.2f", CL * 100), ", ",
      sprintf("%2.2f", CU * 100), ")"
    ),
    P.value = p_value
  )
}

final_results <- do.call(rbind, results)
print(final_results)

write_csv(final_results, paste0(output_dir,"oee_07_edupgs_edu_phenotype.csv"))



# PGS edu + Edu vs. Mortality ---------------------------------------------

traits <- c("death_10yr", "death_cvd_10yr", "death_cancer_10yr")
results <- list()

for (trait in traits) {
  
  formula <- as.formula(paste(trait, "~ prs_edu_edu + sex + age + I(age^2) + adi2010_2 + factor(birth_year) + V3 + V4 + V5 + V6 + V7"))
  
  model <- rqlm(formula, data = dt, family = "poisson", eform = TRUE)
  
  coef_est <- c(1.00, model[2:6, "exp(coef)"])
  p_value <- c(NA, model[2:6, "P-value"])
  se <- c(0.00, model[2:6, "SE"])
  CL <- c(1.00, model[2:6, "CL"])
  CU <- c(1.00, model[2:6, "CU"])
  
  N_cat <- table(dt_tmp$prs_edu_edu)
  N_cat_vec <- as.numeric(N_cat)
  
  N_case_cat <- tapply(
    dt_tmp[[trait]] == 1,
    dt_tmp$prs_edu_edu,
    sum,
    na.rm = TRUE
  )
  N_case_vec <- as.numeric(N_case_cat)
  
  results[[trait]] <- data.frame(
    Traits = trait,
    PGS = rep(c("Low", "Middle", "High"), each = 2),
    Edu = rep(c("Short", "Long"), 3),
    Category = c("1st", "2nd", "3rd", "4th", "5th", "6th"),
    Estimate = coef_est,
    SE = se,
    CI.lower = CL,
    CI.upper = CU,
    N.use    = N_cat_vec,
    N.case   = N_case_vec,
    N = paste0(N_case_vec, "/", N_cat_vec),
    all      = paste0(
      sprintf("%2.2f", coef_est), " (",
      sprintf("%2.2f", CL), ", ",
      sprintf("%2.2f", CU), ")"
    ),
    P.value = p_value
  )
}

final_results <- do.call(rbind, results)
print(final_results)

write_csv(final_results, paste0(output_dir,"oee_08_edupgs_edu_mortality.csv"))



# Edu vs. lifestyle Adjusted for Edu PGS -----------------------

# Univ Comp: 26.0%
dt <- dt |>
  mutate(prs_edu_binary = cut(prs_edu, breaks = c(-Inf, quantile(prs_edu, 0.74), Inf), labels = c(0, 1)))


dt |>
  group_by(prs_edu_binary) |>
  summarise(count = n(),
            per = count / nrow(dt) * 100)


traits <- c("smoking", "drinking", "less_exercise", "inbalanced_food")
results <- list()

for (trait in traits) {
  
  dt_tmp <- dt[!is.na(dt[[trait]]), ]
  
  formula1 <- as.formula(paste(trait, "~ univ_comp + sex + age + I(age^2) + adi2010_2 + factor(birth_year) + V3 + V4 + V5 + V6 + V7"))
  formula2 <- as.formula(paste(trait, "~ prs_edu_binary + sex + age + I(age^2) + adi2010_2 + factor(birth_year) + V3 + V4 + V5 + V6 + V7"))
  
  model1 <- rqlm(formula1, data = dt_tmp, family = "gaussian", digits = 100)
  model2 <- rqlm(formula2, data = dt_tmp, family = "gaussian", digits = 100)
  
  res1 <- data.frame(
    Traits = trait,
    Variable = "univ_comp",
    Estimate = model1["univ_comp", "coef"],
    SE = model1["univ_comp", "SE"],
    CI.lower = model1["univ_comp", "CL"],
    CI.upper = model1["univ_comp", "CU"],
    P.value = model1["univ_comp", "P-value"]
  )
  
  res2 <- data.frame(
    Traits = trait,
    Variable = "prs_edu_binary",
    Estimate = model2["prs_edu_binary", "coef"],
    SE = model2["prs_edu_binary", "SE"],
    CI.lower = model2["prs_edu_binary", "CL"],
    CI.upper = model2["prs_edu_binary", "CU"],
    P.value = model2["prs_edu_binary", "P-value"]
  )
  
  res_tmp <- rbind(res1, res2)
  
  res_tmp$N.used <- nrow(dt_tmp)
  res_tmp$N.case <- sum(dt_tmp[[trait]] == 1, na.rm = TRUE)
  
  res_tmp$all <- paste0(
    sprintf("%2.2f", res_tmp$Estimate * 100), " (",
    sprintf("%2.2f", res_tmp$CI.lower * 100), ", ",
    sprintf("%2.2f", res_tmp$CI.upper * 100), ")"
  )
  
  results[[trait]] <- res_tmp
}

final_results <- do.call(rbind, results)
print(final_results)

write_csv(final_results, paste0(output_dir,"oee_09_edu_lifestyle.csv"))


# Edu vs. Phenotype Adjusted for Edu PGS ---------------------------------------------------


# Lifestyles
traits <- c("obesity", "ht", "dyslipid", "dm")
results <- list()

for (trait in traits) {
  
  dt_tmp <- dt[!is.na(dt[[trait]]), ]
  
  formula1 <- as.formula(paste(trait, "~ univ_comp + sex + age + I(age^2) + adi2010_2 + factor(birth_year) + V3 + V4 + V5 + V6 + V7"))
  formula2 <- as.formula(paste(trait, "~ prs_edu_binary + sex + age + I(age^2) + adi2010_2 + factor(birth_year) + V3 + V4 + V5 + V6 + V7"))
  
  model1 <- rqlm(formula1, data = dt_tmp, family = "gaussian", digits = 100)
  model2 <- rqlm(formula2, data = dt_tmp, family = "gaussian", digits = 100)
  
  res1 <- data.frame(
    Traits = trait,
    Variable = "univ_comp",
    Estimate = model1["univ_comp", "coef"],
    SE = model1["univ_comp", "SE"],
    CI.lower = model1["univ_comp", "CL"],
    CI.upper = model1["univ_comp", "CU"],
    P.value = model1["univ_comp", "P-value"]
  )
  
  res2 <- data.frame(
    Traits = trait,
    Variable = "prs_edu_binary",
    Estimate = model2["prs_edu_binary", "coef"],
    SE = model2["prs_edu_binary", "SE"],
    CI.lower = model2["prs_edu_binary", "CL"],
    CI.upper = model2["prs_edu_binary", "CU"],
    P.value = model2["prs_edu_binary", "P-value"]
  )
  
  res_tmp <- rbind(res1, res2)
  
  res_tmp$N.used <- nrow(dt_tmp)
  res_tmp$N.case <- sum(dt_tmp[[trait]] == 1, na.rm = TRUE)
  
  res_tmp$all <- paste0(
    sprintf("%2.2f", res_tmp$Estimate * 100), " (",
    sprintf("%2.2f", res_tmp$CI.lower * 100), ", ",
    sprintf("%2.2f", res_tmp$CI.upper * 100), ")"
  )
  
  results[[trait]] <- res_tmp
}

final_results <- do.call(rbind, results)
print(final_results)

write_csv(final_results, paste0(output_dir,"oee_10_edu_phenotype.csv"))



# Edu vs. Mortality Adjusted for Edu PGS ---------------------------------------------------

traits <- c("death_10yr", "death_cvd_10yr", "death_cancer_10yr")
results <- list()

for (trait in traits) {
  
  dt_tmp <- dt[!is.na(dt[[trait]]), ]
  
  formula1 <- as.formula(paste(trait, "~ univ_comp + sex + age + I(age^2) + adi2010_2 + factor(birth_year) + V3 + V4 + V5 + V6 + V7"))
  
  formula2 <- as.formula(paste(trait, "~ prs_edu_binary + sex + age + I(age^2) + adi2010_2 + factor(birth_year) + V3 + V4 + V5 + V6 + V7"))
  
  model1 <- rqlm(formula1, data = dt_tmp, family = "poisson", eform = TRUE)
  model2 <- rqlm(formula2, data = dt_tmp, family = "poisson", eform = TRUE)
  
  res1 <- data.frame(
    Traits = trait,
    Variable = "univ_comp",
    Estimate = model1["univ_comp", "exp(coef)"],
    SE = model1["univ_comp", "SE"],
    CI.lower = model1["univ_comp", "CL"],
    CI.upper = model1["univ_comp", "CU"],
    P.value = model1["univ_comp", "P-value"]
  )
  
  res2 <- data.frame(
    Traits = trait,
    Variable = "prs_edu_binary",
    Estimate = model2["prs_edu_binary", "exp(coef)"],
    SE = model2["prs_edu_binary", "SE"],
    CI.lower = model2["prs_edu_binary", "CL"],
    CI.upper = model2["prs_edu_binary", "CU"],
    P.value = model2["prs_edu_binary", "P-value"]
  )
  
  res_tmp <- rbind(res1, res2)
  
  res_tmp$N.used <- nrow(dt_tmp)
  res_tmp$N.case <- sum(dt_tmp[[trait]] == 1, na.rm = TRUE)
  
  res_tmp$all <- paste0(
    sprintf("%2.2f", res_tmp$Estimate), " (",
    sprintf("%2.2f", res_tmp$CI.lower), ", ",
    sprintf("%2.2f", res_tmp$CI.upper), ")"
  )
  
  results[[trait]] <- res_tmp
}


final_results <- do.call(rbind, results)
print(final_results)

write.csv(final_results, paste0(output_dir, "oee_11_edu_poisson.csv"))
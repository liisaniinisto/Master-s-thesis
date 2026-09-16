######################
## setup
#######################

library(data.table)
library(dplyr)
library(tidyr)
library(bigrquery)
library(survey)
library(ggplot2)
library(forcats)
library(logistf)
library(colorspace)
library(stringr)

#######################---------------------------------------------------------
# Functions
#######################---------------------------------------------------------

as_rate = function(p, base){
  paste0(format(p*base, scientific=F, big.mark = " ", digits=1L), " in ", 
         format(base, scientific=F, big.mark = " "))
}

# Clopper-Pearson conf.int. based on the equation presented in thesis Ch. 4.1
clopper_pearson = function(x, n, conf.level){
  alpha = 1 - conf.level
  if(x==0){
    c(0, (1+(n-x)/((x+1)*qf(1-alpha/2, 2*(x+1), 2*(n-x))))^-1)
  }
  else if(x==n){
    c((1+(n-x+1)/(x*qf(alpha/2, 2*x, 2*(n-x+1))))^-1, 1)
  }
  else {
  c(
    (1+(n-x+1)/(x*qf(alpha/2, 2*x, 2*(n-x+1))))^-1,
    (1+(n-x)/((x+1)*qf(1-alpha/2, 2*(x+1), 2*(n-x))))^-1
  )
  }
}
#######################---------------------------------------------------------
# Data
#######################---------------------------------------------------------

# Read genotype information matrix (confidential)
genotypes <- fread("/genotype_matrix.txt")

# Reading information table on the variants (Source: ClinVar)
variant_info <- fread("~/RAS_variantit.txt")

# Read all endpoints (individual level, confidential)
endpoints <- tbl(sandbox, "endpoint_cohorts_r13_v3")

# Read endpoint definition file
endpoint_definitions <- fread("/finngen_R13_endpoint_definitions_3.0.txt")

# Read covariate information (individual level, confidential)
phenotype_data <- fread("/finngen_R13_minimum_extended_2.0.txt.gz")

# collecting summary stats of Finland's age and sex dist (source: Tilastokeskus 2020)
pop_counts_age_sex <- fread("/vaestorakenne_suomi_2020.csv") %>%
  gather(key = "age_group", value = "n", -Sukupuoli) %>%
  mutate(SEX = fct_recode(
    Sukupuoli,
    "male"   = "Miehet",
    "female" = "Naiset"
  )) %>%
  select(-Sukupuoli) %>%
  mutate(
    osuus = n / sum(n),
    osuus = ifelse(SEX == "male", -osuus, osuus),
    source = "Finnish population"
  )



################################################################################
# Preparing data tables for the analysis
################################################################################

# cases and controls for all endpoints (wide format) 
endpoints.wide <- endpoints %>%
  filter(CONTROL_CASE_EXCL <= 1) %>% # removing individuals who fulfill the exclusion criteria
  select(FINNGENID, ENDPOINT, CONTROL_CASE_EXCL) %>%
  pivot_wider(names_from = ENDPOINT, values_from = CONTROL_CASE_EXCL)


# create an individual level table where  ras_variant_status = 1, if >0 RAS-variants are observed
#                                         ras variant_status = 0, if 0 RAS-variants are observed
sample_RAS_status <- genotypes %>%
  mutate(ras_variant_status = ifelse(rowSums(across(-1), na.rm = TRUE) > 0, 1, 0))%>%
  select(FINNGENID, ras_variant_status)

# Dividing age into bins of 5 and 10 years and extracting covariate information
age_groups <- unique(pop_counts_age_sex$age_group)
age_groups2 <- c("0 - 9", "10 - 19", "20 - 29", "30 - 39", "40 - 49", "50 - 59", "60 - 69", "70 - 79", "80 -")


sample_data <- sample_RAS_status %>%
  left_join(endpoints.wide %>%                    # extract Q17_NOONAN endpoint for all FinnGen individuals
              select(FINNGENID, Q17_NOONAN)%>%
              collect(), by="FINNGENID")%>%
  left_join(phenotype_data %>%                    # extract sex and sample collection age for all FinnGen individuals
              select(FINNGENID, BL_AGE, SEX)%>%
              mutate(
                age_group = cut(
                  BL_AGE,
                  breaks = c(seq(0, 85, by = 5), Inf),
                  labels = age_groups,
                  right = TRUE
                ),
                age_group2 = cut(
                  BL_AGE,
                  breaks = c(seq(0, 80, by = 10), Inf),
                  labels = age_groups2,
                  right = TRUE)
              ), by="FINNGENID")

colSums(is.na(sample_data))


################################################################################
# Summary tables by variant and gene (Thesis Ch. 6.4)
################################################################################

# Finding individuals with more than one effective RAS-variant (or 2 alleles of one variant):
RAS_variants_multi <- genotypes %>%
  filter(rowSums(across(-1), na.rm = TRUE) > 1)%>%
  select(FINNGENID, where(~ is.numeric(.x) && sum(.x, na.rm = TRUE) > 0))

nrow(RAS_variants_multi)

# Recoding the genotype data so that 
# 0 = no effective alleels and 
# 1 = one or more effective alleels
lookup <- c('0'=0L, '1'=1L, '2'=1L)

genotypes_binary <-
  genotypes %>%
  mutate(across(-1, ~ lookup[as.character(.)]))

# Frequencies of RAS-variants seen in genotype file (Thesis Appendix A)
variant_freq <- genotypes_binary %>%
  summarize(across(-1, ~sum(.x, na.rm = TRUE))) %>%
  gather(key = "variant", value = "FinnGen_freq.") %>%
  filter(FinnGen_freq. > 0)%>%
  left_join(variant_info, by = "variant")%>%
  arrange(desc(FinnGen_freq.))

# Frequencies of genes (Thesis Table 5)
gene_freq <-
  variant_freq %>%
  group_by(variant_freq$`Gene(s)`)%>%
  summarise('Varianttien lkm.' = n_distinct(variant, na.rm = TRUE), 'Yksilöiden lkm.' = sum(FinnGen_freq., na.rm = TRUE))%>%
  mutate('Selitysosuus (%)'=round(`Yksilöiden lkm.`/sum(`Yksilöiden lkm.`)*100, 1)) %>%
  rename(Geeni = 'variant_freq$`Gene(s)`')%>%
  arrange(desc(`Yksilöiden lkm.`))


################################################################################
# FINNGEN DESCRIPTIVE STATISTICS (Thesis Ch. 6.1)
################################################################################
#-------------------------------------------------------------------------------
# Sample age and sex distribution
#-------------------------------------------------------------------------------

summary_counts_sex = sample_data %>%
  count(SEX)%>%
  mutate('%-osuus' = format(n / sum(n)*100, digit=2L))

(median_age = median(sample_data$BL_AGE, na.rm=T))

# Counts into proportions and manipulating the values of male
# in the negative axis to aid in visualization
summary_counts_age_sex <- sample_data %>%
  group_by(age_group)%>%
  count(age_group, SEX) %>%
  drop_na() %>% 
  ungroup() %>%
  mutate(
    osuus = n / sum(n),
    osuus = ifelse(SEX == "male", -osuus, osuus),
    source = "FinnGen"
  )

# (Thesis Figure 3)
p1 <- ggplot()+
  geom_col(data = summary_counts_age_sex,
           aes(x = osuus, y = age_group, fill = SEX),
           width = 0.5) +
  scale_fill_manual(
    values = c("male" = "#1f5cff", "female" = "#ff6b4a"),
    labels = c("male" = "Miehet (FinnGen)", "female" = "Naiset (FinnGen)")
  ) +
  scale_x_continuous(
    labels = function(x) format(abs(x), big.mark = " ", scientific = FALSE)
  ) +
  geom_col(data = pop_counts_age_sex,
           aes(x = osuus, y = age_group, fill = SEX),
           width = 1,
           color="grey",
           alpha = 0.3)+
  labs(
    x = "Suhteellinen osuus FinnGenissä/väestössä",
    y = "Ikäryhmä",
    fill = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    panel.grid.major.y = element_blank()
  )


################################################################################
# RASOPATHY PREVALENCE (Thesis Ch. 6.3)
################################################################################
#-------------------------------------------------------------------------------
# Prevalence based on Noonan diagnosis (ICD-10 Q87.14)
#-------------------------------------------------------------------------------

# Prevalence based on Q17_NOONAN endpoint (Noonan-diagnosis)
sum(sample_data$Q17_NOONAN, na.rm = T)
(prev_diag <- mean(sample_data$Q17_NOONAN, na.rm = T))
as_rate(prev_diag, 100000)

clopper_pearson(sum(sample_data$Q17_NOONAN, na.rm = T), nrow(sample_data), 0.95)

#-----------------------------------------------------------------------------
# Prevalence based on RAS-variant frequency
#-----------------------------------------------------------------------------

sum(sample_data$ras_variant_status, na.rm = T)
(prev_RAS <- mean(sample_data$ras_variant_status))
as_rate(prev_RAS, 100000)

clopper_pearson(sum(sample_data$ras_variant_status, na.rm = T), nrow(sample_data), 0.95)

# Cross-table of Noonan-diagnosis and having ras-variant
xtabs(~ ras_variant_status + Q17_NOONAN, data = sample_data)

#-----------------------------------------------------------------------------
# Comparing prevalences based on diagnosis and RAS-variant frequency grouped by age and sex
#-----------------------------------------------------------------------------

# checking if there is a statistical difference between sexes for Noonan diagnosis or RAS-status
chisq.test(sample_data$Q17_NOONAN, sample_data$SEX)
chisq.test(sample_data$ras_variant_status, sample_data$SEX)

# prevalence of Noonan diagnosis and RAS-variants by sex
summary_prev_sex <- sample_data %>%
  group_by(SEX) %>%
  summarise(
    n_total = n(), # total number of individuals 
    n_diag = sum(Q17_NOONAN == 1, na.rm = TRUE), # Number of Noonan-diagnosed ind.
    n_ras = sum(ras_variant_status == 1, na.rm = TRUE), # Number of carriers of RAS-variants
    prev_diag = n_diag / n_total*100,
    L95_dia = clopper_pearson(n_diag, n_total, 0.95)[1]*100,
    U95_dia = clopper_pearson(n_diag, n_total, 0.95)[2]*100,
    prev_ras = n_ras / n_total*100,
    L95_ras = clopper_pearson(n_ras, n_total, 0.95)[1]*100,
    U95_ras = clopper_pearson(n_ras, n_total, 0.95)[2]*100
  )%>%
  drop_na()

# prevalence of Noonan diagnosis and RAS-variants by age
summary_prev_age <- sample_data %>%
  group_by(age_group2) %>%
  summarise(
    n_total = n(),
    n_dia = sum(Q17_NOONAN == 1, na.rm = TRUE),
    n_ras = sum(ras_variant_status == 1, na.rm = TRUE),
    prev_dia = n_dia / n_total*100,
    L95_dia = clopper_pearson(n_dia, n_total, 0.95)[1]*100,
    U95_dia = clopper_pearson(n_dia, n_total, 0.95)[2]*100,
    prev_ras = n_ras / n_total*100,
    L95_ras = clopper_pearson(n_ras, n_total, 0.95)[1]*100,
    U95_ras = clopper_pearson(n_ras, n_total, 0.95)[2]*100
  )

summary_prev_age_long<-summary_prev_age %>%
  select(-c(n_total, n_dia, n_ras))%>%
  pivot_longer(cols=-1, names_pattern = "(.*)_(...)$", names_to = c("limit", "name")) %>% 
  mutate(limit=ifelse(limit=="", "value", limit)) %>%
  pivot_wider(id_cols = c(age_group2, name), names_from = limit, values_from = value, names_repair = "check_unique")%>%
  drop_na()

# (Thesis Figure 5)
p2 = summary_prev_age_long %>% ggplot(aes(x=age_group2, y=prev, fill=name))+
  geom_col( 
           position = position_dodge(), width = 0.9,
           alpha = 0.8)+
  labs(x="Ikäryhmä FinnGen näytteenotolle", 
       y = "Vallitsevuus (%)",
       fill = NULL)+
  scale_fill_manual(
    values = c("#1f5cff", "#ff6b4a"),
    labels = c("Noonan-diagnoosi", "Rasopatiavariantti")
  )+
  geom_errorbar(aes(ymin = L95, ymax = U95),
                position = position_dodge(), width = 0.9, colour = '#676765')+
  theme_minimal(base_size = 12)+
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 12))+
  geom_text(
    data = summary_prev_age%>%drop_na()%>%select(n_total, age_group2),
    aes(x = age_group2, y = 0.47, label = paste0("N = ", n_total)),
    inherit.aes = FALSE
  )

#-----------------------------------------------------------------------------
# # Weighted estimate for Noonan diagnosis prevalence
#-----------------------------------------------------------------------------

# by age

# calculating population proportions by age groups
pop_totals_age <- pop_counts_age_sex %>%
  group_by(age_group)%>%
  summarise(Freq = sum(n))%>%
  ungroup()%>%
  mutate(w = Freq/sum(Freq))

# calculating sample prevalence by age group
summary_prev_age <- sample_data %>%
  group_by(age_group) %>%
  drop_na()%>%
  summarise(
    prev_sample = mean(Q17_NOONAN == 1, na.rm = TRUE),
    .groups ="drop"
  )%>%
  left_join(pop_totals_age, by = "age_group")

summary_prev_age$age_group = factor(summary_prev_age$age_group)

# weighted prevalence (manually)
(prev_diag_w <- sum(summary_prev_age$w*summary_prev_age$prev_sample))
as_rate(prev_diag_w, 100000)

# weighted prevalence (using survey package)

design_age <- svydesign(
  ids = ~1,
  data = subset(sample_data, !is.na(age_group) & !is.na(Q17_NOONAN))
)

design_post_age <- postStratify(
  design_age,
  ~age_group,
  summary_prev_age[, c("age_group", "Freq")]
)

svymean(~Q17_NOONAN, design_post_age)
confint(svymean(~Q17_NOONAN, design_post_age))


#-----------------------------------------------------------------------------
# # Weighted estimate for RAS-variant prevalence
#-----------------------------------------------------------------------------

# by endpoint E4_ENDOGLAND


summary_prev_E4 <- sample_data %>%
  left_join(endpoints.wide %>%
              select(E4_ENDOGLAND, FINNGENID)%>%
              collect(), by = "FINNGENID")%>%
  mutate(E4_ENDOGLAND = as.factor(E4_ENDOGLAND))%>%
  group_by(E4_ENDOGLAND) %>%
  summarise(
    prev_ras = mean(ras_variant_status == 1, na.rm = TRUE),
    Freq = sum(ras_variant_status),
    .groups ="drop"
  )%>%
  drop_na()

# Prevalence of E4_ENDOGLAND in population (source: risteys.finngen.fi)
prev_E4_pop <- 0.0299
pop_totals_E4 <- data.frame(
  E4_ENDOGLAND = factor(c(0, 1)),
  Freq = c((1-prev_E4_pop)*10000, prev_E4_pop*10000)
)

# Weighted estimate for RAS-variant prevalence
(prev_RAS_w <- summary_prev_E4$prev_ras[summary_prev_E4$E4_ENDOGLAND == 0] * (1 - prev_E4_pop) +
    summary_prev_E4$prev_ras[summary_prev_E4$E4_ENDOGLAND == 1] * prev_E4_pop)


# weighted prevalence (using survey package)

design_E4 <- svydesign(
  ids = ~1,
  data = subset(sample_data %>%
                  left_join(endpoints.wide %>%
                              select(E4_ENDOGLAND, FINNGENID)%>%
                              collect()%>%
                              mutate(E4_ENDOGLAND = as.factor(E4_ENDOGLAND)), by = "FINNGENID"), !is.na(E4_ENDOGLAND) & !is.na(ras_variant_status)
                              )
)

design_post_E4 <- postStratify(
  design_E4,
  ~E4_ENDOGLAND,
  pop_totals_E4
)

svymean(~ras_variant_status, design_post_E4)
confint(svymean(~ras_variant_status, design_post_E4))


################################################################################
# PheWAS (thesis Ch. 6.2)
################################################################################

#-----------------------------------------------------------------------------
# # Preparing data and functions needed for pheWAS
#-----------------------------------------------------------------------------

# calculate the counts for each endpoint
endpoint_counts <- endpoints.wide %>%
  summarise(across(-FINNGENID, ~sum(.x, na.rm =T))) %>%
  collect()

# Endpoints which have more than 20 observations in FinnGen are used in PheWAS
ep.names_selected <- endpoint_counts%>%
  gather()%>%
  filter(value > 20)%>%
  pull(key)

# Filter endpoint and age data for PheWAS-analysis
endpoints.long <- endpoints %>%
  filter(CONTROL_CASE_EXCL <= 1) %>%
  select(FINNGENID, ENDPOINT, CONTROL_CASE_EXCL, AGE) 

# Table with covariates sex and RAS-status
ras_status_sex <- sample_data %>%
  select(FINNGENID, ras_variant_status, SEX)

# Dividing endpoints to chunks to improve computational speed
chunk_size <- 25 
endpoint_chunks <- split(
  ep.names_selected, 
  ceiling(seq_along(ep.names_selected) / chunk_size) 
)

# Function for extracting logistf output values
extract_logistf <- function(fit, endpoint_name = NA) {
  
  k <- length(fit$coefficients)
  
  data.frame(
    endpoint    = rep(endpoint_name, k),
    term        = names(fit$coefficients),
    beta        = fit$coefficients,
    se          = sqrt(diag(fit$var)),
    OR          = exp(fit$coefficients),
    p_penalized = fit$prob,              
    ci_lower_95 = fit$ci.lower,          
    ci_upper_95 = fit$ci.upper,
    N           = rep(fit$n, k),
    events      = rep(sum(fit$y), k),          
    conv_LL     = fit$conv[1],
    conv_score  = fit$conv[2],
    conv_beta   = fit$conv[3],
    loglik_null = rep(fit$loglik[1], k),
    loglik_full = rep(fit$loglik[2], k),
    row.names   = NULL
  )
}

#-----------------------------------------------------------------------------
# # PheWAS analysis
#-----------------------------------------------------------------------------

start.time <- Sys.time()
all_results <- list()
failed_chunks <- character()

for(i in seq_along(endpoint_chunks)){
  
  chunk = endpoint_chunks[[i]]
  chunk_results <- list()
  
  message("processing chunk ", i)
  
  # collecting endpoint and age data from BigQuery
  dat_chunk <- tryCatch({
    endpoints.long %>%
    filter(ENDPOINT %in% chunk) %>%
    compute()%>%
    collect()
  }, error = function(e) {
    message("Collection failed for chunk ", i)
    failed_chunks <- c(failed_chunks, paste(chunk, collapse = ","))
    return(NULL)
  })
  
  if(is.null(dat_chunk)) next
  
  # Forming the design matrix
  dat_chunk <- ras_status_sex %>%
    left_join(dat_chunk, by="FINNGENID")
  
  # Fitting Firth's penalized regression for each endpoint separately (by chunk)
  for(j in chunk){
    
    dat_endpoint <- dat_chunk %>%
      filter(ENDPOINT == j)
    
    # Fitted with all covariates
    fit = try(logistf(CONTROL_CASE_EXCL ~ ras_variant_status + SEX + AGE, family = binomial, data = dat_endpoint), silent=T)
    if (!inherits(fit, "try-error")) {
      results <- extract_logistf(fit, endpoint_name = j)
      
    }else{
      # Fitted without SEX due to separation issues
      fit=try(logistf(CONTROL_CASE_EXCL ~ ras_variant_status + AGE, family = binomial, data = dat_endpoint), silent=T)
      if (!inherits(fit, "try-error")) {
        results <- extract_logistf(fit, endpoint_name = j)%>%
          mutate(sex_specific = 1)
      }
      
    }
    chunk_results[[j]] <- results
    
  }
  all_results[[as.character(i)]] <- bind_rows(chunk_results)
  message(round(Sys.time() - start.time,2))
}
round(Sys.time() - start.time,2)

result.file <- bind_rows(all_results)


#-----------------------------------------------------------------------------
# # PheWAS result file generation
#-----------------------------------------------------------------------------

# adding the limits for convergence
result.file <- result.file %>%
  mutate(converged =
           abs(conv_LL) < 1e-5 &
           abs(conv_score) < 1e-5 &
           abs(conv_beta) < 1e-5)

# Summary table including only association of interest: ras-status
summary_results <- result.file %>%
  filter(term == "ras_variant_status") %>%
  mutate(
    OR       = round(exp(beta), 2),
    CI_lower_95 = round(exp(ci_lower_95), 2),
    CI_upper_95 = round(exp(ci_upper_95), 2)
  ) %>%
  mutate(p_adj = p.adjust(p_penalized, method = "BH")) %>%  #FDR
  select(
    endpoint,
    events,
    OR,
    CI_lower_95,
    CI_upper_95,
    p_penalized,
    p_adj, 
    sex_specific,
    converged, 
    run_without_age
  ) %>%
  arrange(p_penalized)

# check if there are endpoints which did not converge
ep.names_nonconv <- summary_results%>%
  filter(!converged)%>%
  pull(endpoint)

#-----------------------------------------------------------------------------
# # PheWAS rerun for P16 endpoints
#-----------------------------------------------------------------------------

# (non-converged endpoints and) endpoints starting with P16 (fetus and newborn related endpoints)
# were run again without age as covariate

ep.names_P16 <- summary_results%>%
  filter(str_starts(endpoint, "P16"))%>%
  pull(endpoint)

# PheWAS-analysis

start.time <- Sys.time()
all_results <- list()
failed_endpoints <- character()
  
dat_P16_endpoints <- tryCatch({
  endpoints.wide %>%
    select(all_of(ep.names_P16), FINNGENID) %>%
    compute()%>%
    collect()
}, error = function(e) {
  message("Collection failed for endpoint ", i)
  failed_endpoints <- c(failed_endpoints, paste(ep.names_P16, collapse = ","))
  return(NULL)
})
  
dat_P16_endpoints <- ras_status_sex %>%
   left_join(dat_P16_endpoints, by="FINNGENID")
  
for(i in seq_along(ep.names_P16)){
  fit = try(logistf(dat_P16_endpoints[[ep.names_P16[i]]] ~ ras_variant_status + SEX, family = binomial, data = dat_P16_endpoints), silent=T)
  if (!inherits(fit, "try-error")) {
      results <- extract_logistf(fit, endpoint_name = ep.names_P16[i])%>%
        mutate(run_without_age = TRUE)
  }
    
  all_results[[ep.names_P16[i]]] <- bind_rows(results)
}
round(Sys.time() - start.time,2)

result.file_P16 <- bind_rows(all_results)

result.file_P16 <- result.file_P16 %>%
  mutate(converged =
           abs(conv_LL) < 1e-5 &
           abs(conv_score) < 1e-5 &
           abs(conv_beta) < 1e-5)


summary_results_P16 <- result.file_P16 %>%
  filter(term == "ras_variant_status") %>%
  mutate(
    OR       = round(exp(beta), 2),
    CI_lower_95 = round(exp(ci_lower_95), 2),
    CI_upper_95 = round(exp(ci_upper_95), 2)
  ) %>%
  mutate(p_adj = p.adjust(p_penalized, method = "BH")) %>%
  select(
    endpoint,
    events,
    OR,
    CI_lower_95,
    CI_upper_95,
    p_penalized,
    p_adj, 
    converged, 
  ) %>%
  arrange(p_penalized)

# check if there are endpoints which did not converge
ep.names_nonconv_P16 <- summary_results_P16%>%
  filter(converged=="FALSE")%>%
  pull(endpoint)

#-----------------------------------------------------------------------------
# # PheWAS result analysis
#-----------------------------------------------------------------------------

# Replace original P16 results with rerun results
summary_results_updated <- summary_results %>%
  filter(!endpoint %in% ep.names_P16) %>%
  bind_rows(summary_results_P16)

# Significant endpoints (using FDR-correction) from PheWAS (Thesis Table 2)
significant_endpoints = summary_results_updated%>%
  filter(p_adj<=0.1)

ep.names_significant <- significant_endpoints%>% pull(endpoint)

# cross-tabulating significant endpoints with ras-variant status (Thesis Table 3)
sig_endpoints_data <- endpoints.wide %>%
  select(FINNGENID, all_of(ep.names_significant)) %>%
  compute()%>%
  collect()

sig_endpoints_data <- ras_status_sex %>%
  left_join(sig_endpoints_data, by="FINNGENID")

all_tabs <- list()
for(i in seq_along(ep.names_significant)){
  formula <- as.formula(
    paste("~ ras_variant_status +", ep.names_significant[i])
  )
  tab <- as.data.frame(xtabs(formula, data = sig_endpoints_data))
  print(tab)
  all_tabs[[i]] <- bind_rows(tab)
}
all_tabs <- bind_rows(all_tabs)


# Manhattan plot (Thesis Figure 4)

# List of ep.names often seen with Noonan-patients
ep.names_interest <- c("I9_NONRHEVALV", "IQ17_SEPTA_DEFEC", "Q17_ASD", "Q17_AVSD",
                       "Q17_RVOTO","Q17_CONGENITAL_STENOSIS_OF_PULMONARY_VALVE", "Q17_CONOTR_DEFEC", "I9_CARDMYOHYP",
                       "I9_HYPERTROCARDMYOP", "E4_SHORT", "E4_SHORT_STATURE_STARTING_BEFORE_BIRTH",
                       "Q17_UNDESCENDED_TESTIS_UNILATERAL", "Q17_BILATERAL_CRYPTORCHIDISM", 
                       "H8_HL_CON_NAS", "H8_HL_CON_NAS", "H8_HL_MIX_NAS", "Q17_CONGENITAL_PTOSIS",
                       "Q17_CONGENITAL_PTOSIS_WIDE", "Q17_NOONAN")


summary_results_updated  <- summary_results_updated  %>%
  left_join(endpoint_definitions%>%select(NAME, LONGNAME, TAGS), by = c("endpoint" = "NAME" ))%>%
  arrange(endpoint)%>%
  mutate(
    logp = -log10(p_penalized),
    endpoint_index = row_number(),
    endpoint_of_interest = endpoint %in% ep.names_interest,
    Kategoria = as.factor(sub("^#", "", sub(",.*", "", TAGS))))%>%
  arrange(p_penalized)
  

p_fdr_cutoff <- max(summary_results_updated$p_penalized[summary_results_updated$p_adj <0.1], na.rm=T)


n_cat <- length(unique(summary_results_updated$Kategoria))
cols <- qualitative_hcl(n_cat, palette = "Dark 3")


p3 <- ggplot(summary_results_updated, aes(x = Kategoria, y = logp, color = Kategoria)) +
  geom_jitter(width = 0.5, height = 0, size =1)+
  geom_hline(yintercept = -log10(p_fdr_cutoff),
             linetype = "dashed",
             colour = "black") +
  labs(x = "Fenotyyppikategoriat",
       y = expression(-log[10](p))) +
  geom_text(
    data = subset(summary_results_updated, p_adj < 0.25),
    aes(label = endpoint),
    size =3,
    vjust = -0.5
  )+
  scale_color_manual(values = cols)+
  theme_bw()+
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "none"
  )

                       
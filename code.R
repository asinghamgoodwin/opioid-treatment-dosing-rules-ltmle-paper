library(tidyverse)
library(mice)
library(mitools)
library(SuperLearner)
library(gam)
library(xtable)
library(ltmle)


# Multiple Imputation -----------------------------------------------------

#set random seed for reproduceability
set.seed(123)

#Note: `patients_with_outcomes_02` is a data frame including all individual-level data, with one entry per patient
#We will be imputing values for all of the variables in `comorbidities`.
#  There is no missing data among the `demog` (demographic) variables
#  We will not be imputing outcomes, but using the relapse_date to help our imputations (the other outcomes are redundant)
#  how many missing values are there? --> no columns are missing for more than 8% of patients

patients_to_impute = patients_with_outcomes_02 %>%
  select(all_of(demog), all_of(treatment_info), all_of(comorbidities), relapse_date)

#we're just running a "fake" imputation to get the matrices set up
setup = mice(patients_to_impute, maxit = 0) 
my_methods = setup$method
my_predMatrix = setup$predictorMatrix

#alter the prediction matrix, to exclude some variables from being predictors during imputation
my_predMatrix[, c("who")] = 0
my_predMatrix[, c("project")] = 0 #exclude b/c completely determined by site
my_predMatrix[, c("medicine")] = 0 #exclude b/c completely determined by trt
my_predMatrix[, c("end_of_detox")] = 0 #exclude b/c completely determined by rand_dt

ALT_patients_to_impute = ALT_patients_with_outcomes_02 %>%
  select(all_of(demog), all_of(treatment_info), all_of(comorbidities), relapse_date)

ALT_patients_imputed_03 = mice(ALT_patients_to_impute, 
                               method = my_methods, 
                               predictorMatrix = my_predMatrix, 
                               nnet.MaxNWts = 10000,
                               m = 5)



# Preparing data for LTMLE analysis format --------------------------------

transform_data_for_ltmle = function(original_weekly_data) {
  weeks_with_outcomes_02 = original_weekly_data
  
  #Remove patients who can't be analyzed
  #only keep columns we need for LTMLE. we'll add in demographics later using imputed datasets
  ltmle_prep1 = weeks_with_outcomes_02 %>% 
    dplyr::select(who, switched_meds, never_initiated, rand_dt, relapse_date, medicine, project,
                  week_of_intervention, relapse_this_week, use_this_week, dose_this_week) %>% 
    filter(!switched_meds & !never_initiated) %>% 
    filter(week_of_intervention <= 24) %>% 
    dplyr::select(-switched_meds, -never_initiated, -rand_dt, -relapse_date)

  #Change data format for LTMLE
  ltmle_prep2 = ltmle_prep1 %>% 
    group_by(who) %>% 
    # if it's the last week we have recorded for this patient AND it's earlier than week 24,
    # mark them as relapsed this week. otherwise leave outcome as-is
    mutate(relapse_this_week = ifelse(week_of_intervention == max(week_of_intervention) &
                                        week_of_intervention < 24,
                                      TRUE,
                                      relapse_this_week),
           #if they've ever had a previous relapse, mark them as relapsed for every following week
           # (because our outcome must be monotonic for the ltmle function)
           relapse_this_week = ifelse(lag(relapse_this_week, default = FALSE) |
                                        lag(relapse_this_week, default = FALSE, n = 2) |
                                        lag(relapse_this_week, default = FALSE, n = 3) |
                                        lag(relapse_this_week, default = FALSE, n = 4) |
                                        lag(relapse_this_week, default = FALSE, n = 5) |
                                        lag(relapse_this_week, default = FALSE, n = 6) |
                                        lag(relapse_this_week, default = FALSE, n = 7) |
                                        lag(relapse_this_week, default = FALSE, n = 8) |
                                        lag(relapse_this_week, default = FALSE, n = 9) |
                                        lag(relapse_this_week, default = FALSE, n = 10) |
                                        lag(relapse_this_week, default = FALSE, n = 11) |
                                        lag(relapse_this_week, default = FALSE, n = 12) |
                                        lag(relapse_this_week, default = FALSE, n = 13) |
                                        lag(relapse_this_week, default = FALSE, n = 14) |
                                        lag(relapse_this_week, default = FALSE, n = 15) |
                                        lag(relapse_this_week, default = FALSE, n = 16) |
                                        lag(relapse_this_week, default = FALSE, n = 17) |
                                        lag(relapse_this_week, default = FALSE, n = 18) |
                                        lag(relapse_this_week, default = FALSE, n = 19) |
                                        lag(relapse_this_week, default = FALSE, n = 20) |
                                        lag(relapse_this_week, default = FALSE, n = 21) |
                                        lag(relapse_this_week, default = FALSE, n = 22) |
                                        lag(relapse_this_week, default = FALSE, n = 23) |
                                        lag(relapse_this_week, default = FALSE, n = 24),
                                      TRUE,
                                      relapse_this_week)) %>%
    ungroup() %>% 
    # Create a default "treatment node" - whether or not **this week's dose is higher than last week's**
    # Note: if our treatment definition is more complex (ex. over a dose threshold), this'll get updated later on.
    mutate(dose_increase_this_week = dose_this_week > lag(dose_this_week, default = 0),
           # BUT if they've already relapsed, count it as no dose increase (can't have treatment after outcome)
           dose_increase_this_week = ifelse(relapse_this_week, FALSE, dose_increase_this_week)
    ) %>% 
    # change TRUE/FALSE to 1/0
    mutate(relapse_this_week = as.numeric(relapse_this_week),
           use_this_week = as.numeric(use_this_week),
           dose_increase_this_week = as.numeric(dose_increase_this_week))
  
  
  #Pivot to wide format
  ltmle_prep3 = ltmle_prep2 %>% 
    # pivot to wide format
    pivot_wider(names_from = week_of_intervention,
                names_glue = "wk{week_of_intervention}.{.value}",
                values_from = c(dose_this_week,
                                use_this_week,
                                relapse_this_week,
                                dose_increase_this_week),
                values_fill = NA) %>% 
    # LTMLE doesn't play nice with NAs (which are generated by pivot_wider for any weeks that weren't documented)
    # so here I'm replacing any NAs with a default value (1 for relapse, 0 for the others).
    # These won't get used in LTMLE computations anyways, since they all occur after a relapse
    mutate_at(vars(contains("relapse_this_week")), ~replace_na(., 1)) %>%
    mutate_at(vars(contains("dose_this_week")), ~replace_na(., 0)) %>%
    mutate_at(vars(contains("use_this_week")), ~replace_na(., 0)) %>% 
    mutate_at(vars(contains("dose_increase_this_week")), ~replace_na(., 0))
  
  # Wide dataset to return:
  ltmle_prep3
}

#Create a set of different cases to analyze

case_attributes = c("medicine", "projects", "dose_threshold", "name", "max_weeks", "responsive")
# where responsive can be TRUE, FALSE (for only looking at any dose increase), 
# ...or "comp" for a comparisson between responsive and any dose increase, or "hybrid" for comparing hybrid vs. static dose

caseA = list("bup", c(27, 30, 51), 32, "Buprenorphine (all 3 projects), responsive dose increase", 24, TRUE)
names(caseA) = case_attributes

caseB = list("bup", c(27, 30, 51), 16, "Buprenorphine (all 3 projects), any dose increase if <16mg", 24, FALSE)
names(caseB) = case_attributes

caseC = list("bup", c(27, 30, 51), 16, "Buprenorphine (all 3 projects), comp: increase if <16mg vs. responsive dose", 24, "comp")
names(caseC) = case_attributes

caseD = list("met", c(27), 150, "Methodone (p27), responsive dose increase", 24, TRUE)
names(caseD) = case_attributes

caseE = list("met", c(27), 100, "Methodone (p27), any dose increase if <100mg", 24, FALSE)
names(caseE) = case_attributes

caseF = list("met", c(27), 100, "Methodone (p27), comp: increase if <100mg vs. responsive dose", 24, "comp")
names(caseF) = case_attributes

caseG = list("bup", c(27, 30, 51), 16, "Buprenorphine (all 3 projects), hybrid rule vs. increase if <16mg", 24, "hybrid")
names(caseG) = case_attributes

caseH = list("bup", c(27, 30, 51), 32, "Buprenorphine (all 3 projects), hybrid rule vs. responsive dose", 24, "hybrid")
names(caseH) = case_attributes

caseI = list("met", c(27), 100, "Methodone (p27), hybrid rule vs. increase if <100mg", 24, "hybrid")
names(caseI) = case_attributes

caseJ = list("met", c(27), 150, "Methodone (p27), hybrid rule vs. responsive dose", 24, "hybrid")
names(caseJ) = case_attributes

caseK = list("bup", c(27, 30, 51), 0, "Buprenorphine (all 3 projects), hybrid rule vs. static dose", 24, "hybrid")
names(caseK) = case_attributes

caseL = list("met", c(27), 0, "Methodone (p27), hybrid rule vs. static", 24, "hybrid")
names(caseL) = case_attributes

cases_NEW = list(caseA, caseB, caseC, caseD, caseE, caseF, caseG, caseH, caseI, caseJ)
cases_HYBRID_v_Static = list(caseK, caseL)


#Function: take in a case, output a dataset and parameters ready for LTMLE
ltmle_case_prep = function(data, case) {
  
  #Step 1: filter the data down to only the medicine, project, and weeks specified by the case
  
  filtered_data = data %>% 
    filter(project %in% case[["projects"]] & medicine == case[["medicine"]]) %>% 
    dplyr::select(-starts_with(paste0("wk", (case[["max_weeks"]]+1):25)))
  
  #assign an absolute max for the medicine (useful in comparisson situations)
  med_max = ifelse(case[["medicine"]] == "bup", 32, 150)
  
  #Step 2: check whether there are enough events in each week to include that week in the analysis
  
  #Are there any weeks (of 4-24) where no-one's outcome changed (no new relapses)?
  #...If so, exclude those weeks from the analysis alltogether.
  outcome_weeks = c()
  for (w in 4:case[["max_weeks"]]) {
    this_week_relapse_col = paste0("wk", w, ".relapse_this_week")
    last_week_relapse_col = paste0("wk", w-1, ".relapse_this_week")
    if (!all(filtered_data[[sym(this_week_relapse_col)]] == filtered_data[[sym(last_week_relapse_col)]], na.rm = TRUE)) {
      outcome_weeks = c(outcome_weeks, w)
    }
  }
  
  #Are there any weeks with 0 people who fit the treatment rule?
  #...If so, exclude those weeks from the Anodes, Lnodes, and Abar.
  #...also, exclude any weeks past week 12.
  weeks_with_treatments = c()
  # only go up to whichever is first, 12 weeks or max_weeks
  for (w in 4:min(case[["max_weeks"]], 12)) {
    increases_this_week = filtered_data[[sym(paste0("wk", w, ".dose_increase_this_week"))]]
    # under_threshold_last_week is always TRUE if our threshold is the max
    under_threshold_last_week = filtered_data[[sym(paste0("wk", w - 1, ".dose_this_week"))]] < case[["dose_threshold"]]
    under_max_last_week = filtered_data[[sym(paste0("wk", w - 1, ".dose_this_week"))]] < med_max
    # use_last_week is always TRUE if we don't actually care about being responsive
    use_last_week = ifelse(rep(case[["responsive"]] == TRUE, nrow(filtered_data)), 
                           filtered_data[[sym(paste0("wk", w - 1, ".use_this_week"))]],
                           rep(TRUE, nrow(filtered_data)))
    
    treatA = sum(increases_this_week & under_threshold_last_week)
    treatB = sum(increases_this_week & use_last_week & under_max_last_week)
    
    if (case[["responsive"]] == "comp") {
      treatment = min(treatA, treatB) #not sure why this isn't actually choosing the min... but it's ok because they're always above 0
    } else if (case[["responsive"]] == "hybrid") {
      treatment = sum(treatA, treatB)
    } else {
      # this is just a hybrid that'll be correct regardless of which treatment we care about
      treatment = sum(increases_this_week & under_threshold_last_week & use_last_week)
    }
    
    #cat(paste0("Treatments this week: ", sum(treatment), "\n"))
    if (treatment > 0) {
      weeks_with_treatments = c(weeks_with_treatments, w)
    }
  }
  
  #Step 4: create the lists of Anodes, Lnodes, Ynodes, and abar
  
  #Anodes (treatment nodes) - each includes the treatment node for this and all previous weeks...
  #...restricted to weeks where a relapse outcome and a treatment were both possible (so, 4-24)
  Anodes = vector("list", length(outcome_weeks))
  for (i in 1:length(outcome_weeks)) {
    Anodes[[i]] <- #c("wk3.dose_increase_this_week", 
      c(paste0("wk", outcome_weeks[1:i][outcome_weeks[1:i] %in% weeks_with_treatments], ".dose_increase_this_week"))
  }
  
  #Lnodes (time-varying covariates) - each includes `dose_this_week` lagged by 1 (3-24)...
  #... and `use_this_week` lagged by 1 for this and all previous weeks (2-24).
  #  (for use this week, include a week if the next week was a week with a treatment)
  Lnodes = vector("list", length(outcome_weeks))
  for (i in 1:length(outcome_weeks)) {
    Lnodes[[i]] <- c(paste0("wk", c(3, outcome_weeks[1:i-1][outcome_weeks[1:i-1] %in% c(weeks_with_treatments-1)]), ".dose_this_week"),
                     paste0("wk", c(2, 3, outcome_weeks[1:i-1][outcome_weeks[1:i-1] %in% weeks_with_treatments]), ".use_this_week"))
  }
  
  #Ynodes (outcome nodes) - `relapse_this_week` for each week...
  #...restricted to weeks where a relapse outcome was possible (so, 4-24)
  Ynodes = paste0("wk", outcome_weeks, ".relapse_this_week")
  
  #abar (treatment rule)
  abar1 = matrix()
  abar0 = matrix()
  ## First part of the rule: did they have use the week before? (this is necessary for RESPONSIVE treatment rules)
  #...for each observation at each time point, was there any opioid use the week before?
  abarA <- filtered_data[, paste0("wk", c(outcome_weeks - 1)[outcome_weeks %in% weeks_with_treatments], ".use_this_week")] == 1
  
  ## Second part of the rule: were they able to be increased? (i.e. under the max dose last week)
  #...for each observation at each time point, were they under the max possible dose the week before?
  abarB <- filtered_data[, paste0("wk", c(outcome_weeks - 1)[outcome_weeks %in% weeks_with_treatments], ".dose_this_week")] < med_max
  
  ## Third part of the rule: were they below the threshold we're interested in? (if no threshold, max do is passed in)
  #...for each observation at each time point, were they under the threshold dose the week before?
  #... (meaning: if we're interested in dose increases for anyone who is under 16mg, 
  #     we "count" them if last week they were under 16, and this week they got a dose increase)
  abarC <- filtered_data[, paste0("wk", c(outcome_weeks - 1)[outcome_weeks %in% weeks_with_treatments], ".dose_this_week")] < case[["dose_threshold"]]
  
  if (case[["responsive"]] == "comp") {
    #"treatment" = dose increase up till threshold = abar1
    abar1 <- I(abarB == TRUE & abarC == TRUE)
    
    #"no treatment" = comp = dynamic dose increase = abar0
    abar0 <- I(abarA == TRUE & abarB == TRUE)
    
  } else if (case[["responsive"]] == "hybrid") {
    
    #B: under the max dose AND [ (C: below our threshold) OR (above our threshold & A: had use last week) ]
    # note - we can leave out the "above our threshold" bit because it's implied by the logic
    abar1 <- I(abarB == TRUE & ((abarC == TRUE) | (abarA)))
    
    # the "no treatment" / comparisson group is either the dynamic dose or the dose increase till threshold
    # and we can tell which by referring to the dose threshold case attribute
    # if the dose threshold case attribute is 0, that's code for "compare against static dose"
    
    if (case[["dose_threshold"]] == med_max) {
      #cases H and J -- dynamic dose
      abar0 <- I(abarA == TRUE & abarB == TRUE)
    } else if (case[["dose_threshold"]] == 0) {
      #cases K and L -- static dose
      abar0 <- matrix(rep(FALSE, length(weeks_with_treatments)*nrow(filtered_data)), ncol = length(weeks_with_treatments))
    } else {
      #cases G and I -- increase up to threshold
      abar0 <- I(abarB == TRUE & abarC == TRUE)
    }
  } else if (case[["responsive"]] == TRUE) {
    abar1 <- I(abarA == TRUE & abarB == TRUE)
    
    #...a matrix of what the treatment would be for every observation at every time point,
    #...in a hypothetical population with no dynamic dose increase (so, all FALSE)
    abar0 <- matrix(rep(FALSE, length(weeks_with_treatments)*nrow(filtered_data)), ncol = length(weeks_with_treatments))
    
  } else {
    abar1 <- I(abarB == TRUE & abarC == TRUE)
    
    #...a matrix of what the treatment would be for every observation at every time point,
    #...in a hypothetical population with no dynamic dose increase (so, all FALSE)
    abar0 <- matrix(rep(FALSE, length(weeks_with_treatments)*nrow(filtered_data)), ncol = length(weeks_with_treatments))
  }
  
  
  #Step 5: create a table showing the counts of different types of patient events at each week
  # for the comparisson tables, treatment_this_week means "dose increase to above THRESHOLD"
  weekly_counts = tibble(week = 1:case[["max_weeks"]]) %>% 
    group_by(week) %>% 
    mutate(patients = sum(filtered_data[, paste0("wk", week, ".relapse_this_week")] == 0),
           dose_increase_this_week = sum(filtered_data[, paste0("wk", week, ".dose_increase_this_week")] == 1),
           use_last_week = sum(filtered_data[, paste0("wk", max(week - 1, 1), ".use_this_week")] == 1),
           # in comparisson situations, this first vbl represents how many got the dose increase to a threshold
           treatment_this_week = sum(filtered_data[, paste0("wk", week, ".dose_increase_this_week")] == 1 &
                                       ifelse(rep(case[["responsive"]] == TRUE, nrow(filtered_data)), 
                                              filtered_data[, paste0("wk", max(week - 1, 1), ".use_this_week")] == 1, 
                                              rep(TRUE, nrow(filtered_data))) &
                                       filtered_data[, paste0("wk", max(week - 1, 1), ".dose_this_week")] < case[["dose_threshold"]]),
           # in comparisson situations, this second vbl represents how many got the dynamic dose treatment (increase, use last week, under max)
           DELETE_or_alt_trt_for_comparisson = sum(filtered_data[, paste0("wk", week, ".dose_increase_this_week")] == 1 &
                                                     filtered_data[, paste0("wk", max(week - 1, 1), ".use_this_week")] == 1 &
                                                     filtered_data[, paste0("wk", max(week - 1, 1), ".dose_this_week")] < med_max),
           #we want to create summary counts, but only want patients who haven't yet relapsed (the [[1,2]])
           dose_mean = aggregate(filtered_data[, paste0("wk", week, ".dose_this_week")],
                                 by = filtered_data[, paste0("wk", week, ".relapse_this_week")],
                                 FUN = mean)[[1,2]],
           dose_median = aggregate(filtered_data[, paste0("wk", week, ".dose_this_week")],
                                   by = filtered_data[, paste0("wk", week, ".relapse_this_week")],
                                   FUN = median)[[1,2]],
           dose_min = aggregate(filtered_data[, paste0("wk", week, ".dose_this_week")],
                                by = filtered_data[, paste0("wk", week, ".relapse_this_week")],
                                FUN = min)[[1,2]],
           dose_max = aggregate(filtered_data[, paste0("wk", week, ".dose_this_week")],
                                by = filtered_data[, paste0("wk", week, ".relapse_this_week")],
                                FUN = max)[[1,2]],
           dose_iqr = aggregate(filtered_data[, paste0("wk", week, ".dose_this_week")],
                                by = filtered_data[, paste0("wk", week, ".relapse_this_week")],
                                FUN = IQR)[[1,2]]
    ) %>% 
    ungroup() %>% 
    mutate(new_relapse_this_week = lag(patients, default = 0) - patients)
  
  #Step 6: return all inputs for the LTMLE function as a list with named elements
  input_names = c("dataset_missing_baseline_covariates",
                  "Anodes",
                  "Lnodes",
                  "Ynodes",
                  "abar1",
                  "abar0",
                  "outcome_weeks",
                  "weeks_with_treatments",
                  "weekly_counts",
                  "name",
                  "max_weeks")
  inputs = list(filtered_data,
                Anodes,
                Lnodes,
                Ynodes,
                abar1,
                abar0,
                outcome_weeks,
                weeks_with_treatments,
                weekly_counts,
                case[["name"]],
                case[["max_weeks"]])
  names(inputs) = input_names
  
  inputs
}

ALT_weekly_data_for_ltmle_04_NEW = list()
for (new_case in cases_NEW) {
  ALT_weekly_data_for_ltmle_04_NEW = c(ALT_weekly_data_for_ltmle_04_NEW,
                                       list(ltmle_case_prep(transform_data_for_ltmle(ALT_weeks_with_outcomes_02), new_case)))
}

ALT_weekly_data_for_ltmle_04_HYBRID_v_Static = list()
for (new_case in cases_HYBRID_v_Static) {
  ALT_weekly_data_for_ltmle_04_HYBRID_v_Static = c(ALT_weekly_data_for_ltmle_04_HYBRID_v_Static,
                                                   list(ltmle_case_prep(transform_data_for_ltmle(ALT_weeks_with_outcomes_02), new_case)))
}



# LTMLE analysis ----------------------------------------------------------

run_ltmle_case = function(prepared_info, mi_baseline_covariates, imputation_num, week_cutoff) {
  #Step 1: general prep
  #...1a: "unpack" the prepared info
  dataset_missing_baseline_covariates = prepared_info$dataset_missing_baseline_covariates %>% dplyr::select(., -medicine, -project)
  Anodes = prepared_info$Anodes
  Lnodes = prepared_info$Lnodes
  Ynodes = prepared_info$Ynodes
  abar1 = prepared_info$abar1
  abar0 = prepared_info$abar0
  outcome_weeks = prepared_info$outcome_weeks
  outcome_weeks = outcome_weeks[outcome_weeks <= week_cutoff]
  name = prepared_info$name
  
  #...1b: "unpack" the imputed baseline data
  baseline_covariates = complete(mi_baseline_covariates, action = imputation_num) %>%
    dplyr::select(all_of(c(demog, comorbidities)), site)
  
  #...1c: add baseline covariates into the dataset
  all_data = inner_join(baseline_covariates,
                        dataset_missing_baseline_covariates,
                        by = "who")
  
  #...1d: create helper variables
  all_cols = colnames(all_data)
  baseline_cols = colnames(baseline_covariates)
  
  # #...1e: create an empty vector where we'll store the output from each week's LTMLE estimation
  txest = vector("list", length(outcome_weeks))
  txvar = vector("list", length(outcome_weeks))
  cntest = vector("list", length(outcome_weeks))
  cntvar = vector("list", length(outcome_weeks))
  est = vector("list", length(outcome_weeks))
  var = vector("list", length(outcome_weeks))
  rrest = vector("list", length(outcome_weeks))
  rrvar = vector("list", length(outcome_weeks))
  
  
  #Step 2: iterate through all outcome weeks to get estimations
  for (i in 1:length(outcome_weeks)) {
    #...2a: prepare the data for use in the LTMLE function:
    #remove any columns that weren't used as A-, L-, or Y- nodes
    unused = setdiff(all_cols, c(baseline_cols, Anodes[[i]], Lnodes[[i]], Ynodes[1:i]))
    
    ltmle_data = all_data %>%
      dplyr::select(-all_of(unused)) %>%
      relocate(baseline_cols, Lnodes[[i]], Anodes[[i]], Ynodes[1:i]) %>% 
      dplyr::select(-who)
    
    ltmle_data = as.data.frame(ltmle_data)
    
    abar_this_time = list(abar1[,1:length(Anodes[[i]]), drop = FALSE], abar0[,1:length(Anodes[[i]]), drop = FALSE])
    
    #...2b: run the ltmle function
    output <- ltmle(
      data = ltmle_data,
      Anodes = Anodes[[i]],
      Lnodes = Lnodes[[i]],
      Ynodes = Ynodes[1:i],
      survivalOutcome = TRUE, #TRUE means the outcome is an event that can only occur once.
      abar = abar_this_time,
      estimate.time = FALSE, #if TRUE, this would just compute a rough estimate of runtime based on first 50 obs
      SL.library = c("SL.glm", "SL.mean", "SL.earth", "SL.xgboost"),
      variance.method = "ic"
    )
    
    #...2c: keep the output
    sumout = summary(output)
    
    txest[i] <- sumout$effect.measures$treatment$estimate
    txvar[i] <- sumout$effect.measures$treatment$std.dev ^ 2
    
    cntest[i] <- sumout$effect.measures$control$estimate
    cntvar[i] <- sumout$effect.measures$control$std.dev ^ 2
    
    est[i] <- sumout$effect.measures$ATE$estimate
    var[i] <- sumout$effect.measures$ATE$std.dev ^ 2
    
    rrest[i] <- sumout$effect.measures$RR$estimate
    rrvar[i] <- sumout$effect.measures$RR$std.dev ^ 2
  }
  
  #Step 3: package output together to return
  estimates = tibble(
    Outcome_week = outcome_weeks,
    Treatment_estimate = txest, 
    Treatment_variance = txvar,
    Control_estimate = cntest,
    Control_variance = cntvar,
    ATE_estimate = est,
    ATE_variance = var,
    RR_estimate = rrest,
    RR_variance = rrvar
  )
  
  estimates
}


run1 = function(case_num, imp_num) {
  result = run_ltmle_case(prepared_info = ALT_weekly_data_for_ltmle_04_NEW[[case_num]],
                          mi_baseline_covariates = ALT_patients_imputed_03,
                          imputation_num = imp_num,
                          week_cutoff = 5)
  write.csv(as.matrix(result), paste0("Case", case_num, "Imp", imp_num, ".csv"))
}

for (case in 1:length(cases_NEW)) {
  for (i in 1:5) {
    run1(case, i)
  }
}



# Pooling estimates across imputations ------------------------------------

combine_results_from_all_imputations = function(case_num) {
  #Convert CSVs to tibbles
  #Merge data from tibbles
  all_estimates = vector("list", 5)
  for (m in 1:5) {
    all_estimates[[m]] = read_csv(file = paste0("Case", case_num,"Imp", m, ".csv")) %>% 
      dplyr::select(-X1)
  }
  
  merged_all = full_join(all_estimates[[1]],
                         all_estimates[[2]],
                         by = "Outcome_week",
                         suffix = c("1", "2"))
  
  for (m in 3:5) {
    merged_all = full_join(
      all_estimates[[m]],
      merged_all,
      by = "Outcome_week",
      suffix = c(as.character(m), "3")
    ) #for some reason this is what works for labeling!
  }
  
  merged_all = merged_all %>%
    dplyr::select("Outcome_week", sort(colnames(.)))
  
  merged_all
}

#Use MIcombine to pool all results
pool_outcomes_from_imputations = function(combined_results) {
  merged_all = combined_results
  
  pooled_results = tibble(
    Outcome_week = numeric(),
    Treatment_estimate = numeric(),
    Treatment_variance = numeric(),
    Treatment_lowerCI = numeric(),
    Treatment_upperCI = numeric(),
    Control_estimate = numeric(),
    Control_variance = numeric(),
    Control_lowerCI = numeric(),
    Control_upperCI = numeric(),
    ATE_estimate = numeric(),
    ATE_variance = numeric(),
    ATE_lowerCI = numeric(),
    ATE_upperCI = numeric(),
    RR_estimate = numeric(),
    RR_variance = numeric(),
    RR_lowerCI = numeric(),
    RR_upperCI = numeric(),
  )
  
  
  for (week in merged_all$Outcome_week) {
    this_week_only = merged_all %>% filter(Outcome_week == week)
    
    treatment_pool = MIcombine(results = as.list(unlist(
      dplyr::select(this_week_only, starts_with("Treatment_estimate")),
      use.names = FALSE
    )),
    variances = as.list(unlist(
      dplyr::select(this_week_only, starts_with("Treatment_variance")),
      use.names = FALSE
    )))
    
    control_pool = MIcombine(results = as.list(unlist(
      dplyr::select(this_week_only, starts_with("Control_estimate")),
      use.names = FALSE
    )),
    variances = as.list(unlist(
      dplyr::select(this_week_only, starts_with("Control_variance")),
      use.names = FALSE
    )))
    
    ATE_pool = MIcombine(results = as.list(unlist(
      dplyr::select(this_week_only, starts_with("ATE_estimate")),
      use.names = FALSE
    )),
    variances = as.list(unlist(
      dplyr::select(this_week_only, starts_with("ATE_variance")),
      use.names = FALSE
    )))
    
    RR_pool = MIcombine(results = as.list(unlist(
      dplyr::select(this_week_only, starts_with("RR_estimate")),
      use.names = FALSE
    )),
    variances = as.list(unlist(
      dplyr::select(this_week_only, starts_with("RR_variance")),
      use.names = FALSE
    )))
    
    pooled_results = pooled_results %>% 
      add_row(
        Outcome_week = week,
        Treatment_estimate = treatment_pool$coefficients,
        Treatment_variance = treatment_pool$variance,
        Treatment_lowerCI = treatment_pool$coefficients - 1.96 * (treatment_pool$variance ^ .5),
        Treatment_upperCI = treatment_pool$coefficients + 1.96 * (treatment_pool$variance ^ .5),
        Control_estimate = control_pool$coefficients,
        Control_variance = control_pool$variance,
        Control_lowerCI = control_pool$coefficients - 1.96 * (control_pool$variance ^ .5),
        Control_upperCI = control_pool$coefficients + 1.96 * (control_pool$variance ^ .5),
        ATE_estimate = ATE_pool$coefficients,
        ATE_variance = ATE_pool$variance,
        ATE_lowerCI = ATE_pool$coefficients - 1.96 * (ATE_pool$variance ^ .5),
        ATE_upperCI = ATE_pool$coefficients + 1.96 * (ATE_pool$variance ^ .5),
        RR_estimate = RR_pool$coefficients,
        RR_variance = RR_pool$variance,
        RR_lowerCI = RR_pool$coefficients - 1.96 * (RR_pool$variance ^ .5),
        RR_upperCI = RR_pool$coefficients + 1.96 * (RR_pool$variance ^ .5),
      )
  }
  
  pooled_results
}


ltmle_results = vector(mode = "list", length = length(cases_NEW))

for (i in 1:length(cases_NEW)) {
  prepped_data = ALT_weekly_data_for_ltmle_04_NEW[[i]]
  case = cases_NEW[[i]]
  estimatesAll5 = combine_results_from_all_imputations(i)
  pooled = pool_outcomes_from_imputations(estimatesAll5)
  
  ltmle_results[[i]] <- pooled
}
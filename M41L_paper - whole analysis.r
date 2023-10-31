library(data.table)
library(tidyverse)
library(glue)
library(lubridate)
library(tableone)

google_bucket = Sys.getenv('WORKSPACE_BUCKET')

tt = function(x) table(table(x))
lu = function(x) length(unique(x))
g = pillar::glimpse

save_to_PD_and_GCS = function(x, dir = 'M41L_paper', PD = TRUE, GCS = TRUE) {
    stopifnot(is_tibble(x))
    fn = paste0(substitute(x), '.csv')
    google_bucket = Sys.getenv('WORKSPACE_BUCKET')
    if (PD) { write_csv(x = x, file = glue('./{dir}/{fn}')) }
    if (GCS) {
        system(glue('gsutil -m cp ./{dir}/{fn} {google_bucket}/{dir}/'), intern = TRUE)
        system(glue('gsutil ls {google_bucket}/{dir}/{fn}'), intern = TRUE)
    }
}

'%nin%' = function(a, b) (!(a %in% b))

demog_wgs = read_csv('aim1/demog_wgs_pts.csv', show_col_types = FALSE)

read_tsv('aim1/ancestry_preds.tsv', show_col_types = FALSE) %>% 
select(person_id = research_id, ancestry_pred) ->
ancestry_pred

read_tsv('aim1/genomic_metrics.tsv', show_col_types = FALSE) %>%
select(person_id = research_id, dragen_sex_ploidy, biosample_collection_date) ->
genomic_metrics

demog_wgs %>%
left_join(ancestry_pred) %>%
left_join(genomic_metrics) %>%
filter(!is.na(date_of_birth), !is.na(biosample_collection_date)) %>%
mutate(
    age_at_biosample = round(time_length(interval(date_of_birth, biosample_collection_date), unit = 'years')),
    decade_at_biosample = round(age_at_biosample, -1)
) %>%
# glimpse() ->
identity() ->
demog_plus

mutect_results_v6_files = list.files(path = 'mutect2_results_v6/', full.names = TRUE)
n = length(mutect_results_v6_files) %>% print()
groupsize = 1000
pts_by_group = split(seq(n), ceiling(seq(n)/groupsize))
length(pts_by_group)

collector = list()
last_group = length(pts_by_group)

for (group_idx in names(pts_by_group)[1:last_group]) {
    
    message(Sys.time(), ' reading group ', group_idx, ' of ', last_group)
    
    mutect_results_v6_files %>%
    `[`(pts_by_group[[group_idx]]) %>%
    map_df(fread, colClasses = list(character = 'AF')) %>%
    bind_rows() %>%
    filter(Gene.refGene == 'UBA1') ->
    collector[[group_idx]]
}

collector = bind_rows(collector)

write_csv(x = collector, file = 'M41L_paper/SM_in_UBA1_v6_unprocessed.csv')

# mutect_results_v7_files = list.files(path = 'mutect2_results_v7/mutect2_20230423_results/', full.names = TRUE)
# n = length(mutect_results_v7_files) %>% print()
# groupsize = 100
# pts_by_group = split(seq(n), ceiling(seq(n)/groupsize))
# length(pts_by_group)

# read_tsv(file = mutect_results_v7_files[1], comment = '##')

# collector = list()
# last_group = length(pts_by_group)

# for (group_idx in names(pts_by_group)[1:last_group]) {
    
#     message(Sys.time(), ' reading group ', group_idx, ' of ', last_group)
    
#     mutect_results_v7_files %>%
#     `[`(pts_by_group[[group_idx]]) %>%
#     map_df(fread, colClasses = list(character = 'AF')) %>%
# #     read_csv(show_col_types = FALSE) %>%
#     bind_rows() %>%
# #     filter(`Gene.refGene` %in% c(aim1_genes, 'MYD88')) %>%
#     identity() ->
#     collector[[group_idx]]
# }

# collector = bind_rows(collector)

# collector %>%
# # glimpse()
# # group_by(`Func.refGene`) %>% count() %>% arrange(-n)
# # # TODO: could be reasonable to include upstream and downstream
# filter(`Func.refGene` %in% c('exonic', 'splicing')) %>%
# # filter(`Gene.refGene` %in% c(aim1_genes, 'MYD88')) %>%
# group_by(`Gene.refGene`) %>% count() %>% arrange(-n)
# # group_by(`Gene.refGene`) %>% 
# # group_walk(~ write_csv(x = .x, file = paste0('aim1/', .y$`Gene.refGene`, "_v7.csv")))

allvars = fread(file = 'allvars_genes_of_Rob_Ashwin_Caitlyn_lastrow=146019.csv')

allvars %>%
filter(Gene.refGene == 'UBA1') ->
uba1_v7

write_csv(x = uba1_v7, file = 'M41L_paper/SM_in_UBA1_v7_unprocessed.csv')

uba1_col_types = paste0('ciicc', 'ccccc', 'ccccc', 'ccccc', 'ccccc', 'cccic', 'cciii', 'ddiii', 'iii_')

uba1_v6 = read_csv('M41L_paper/SM_in_UBA1_v6_unprocessed.csv', col_types = uba1_col_types)
uba1_v7 = read_csv('M41L_paper/SM_in_UBA1_v7_unprocessed.csv', col_types = uba1_col_types)

uba1 = bind_rows(uba1_v6, uba1_v7)

get_binom_p_val = function(a, b) { binom.test(x = a, n = a + b)$p.value }

uba1 %>%
select(-starts_with('Other')) %>%
filter(DP >= 20, minAD >= 3, F1R2_1 > 0, F1R2_2 > 0, F2R1_1 > 0, F2R1_2 > 0) %>%
# there are no multiallelic SNPs
# filter(!is.na(minAD_alt_multiallelic)) %>%
mutate(
    AF = as.numeric(AF), # can do this safely bc no multiallelic
    binom_p_val = map2(.x = maxAD, .y = minAD, .f = get_binom_p_val),
    person_id = parse_number(Sample)
) %>%
select(
    person_id,
    chr = Chr, pos = Start, ref = Ref, alt = Alt, 
    gene_name = Gene.refGene, func = Func.refGene, aa_change = AAChange.refGene, cosmic70,
    DP, AD, F1R2, F2R1, SB, AF, binom_p_val
) %>%
inner_join(demog_plus, by = 'person_id') %>%
filter(dragen_sex_ploidy == 'XY' | binom_p_val < 0.01) %>%
# glimpse() %>%
identity() ->
sm_in_uba1_genetics_passing

save_to_PD_and_GCS(sm_in_uba1_genetics_passing)

sm_in_uba1_genetics_passing %>%
separate_wider_delim(
    too_few = 'align_start',
    cols = aa_change, 
    delim = regex(':|,'),
    names = paste0(rep(c('gene_name', 'transcript', 'exon', 'dna', 'protein'), times = 2), rep(1:2, each = 5))) %>%
# only rows where these new columns can't be extracted are the splicing variants
# filter(!aa_change_ok, func != 'splicing') %>%
# protein1 always == protein2
# filter(protein1 != protein2) %>%
filter(str_detect(string = protein1, pattern = 'M41[A-Z]')) %>%
# glimpse() %>%
identity() ->
sm_in_uba1_at_m41

save_to_PD_and_GCS(sm_in_uba1_at_m41)

sm_in_uba1_at_m41 %>%
group_by(gene_name, dragen_sex_ploidy, sex_at_birth, dna1, protein1) %>% count() %>% arrange(-n)

sm_in_ubs1_at_m41 = read_csv('M41L_paper/sm_in_uba1_at_m41.csv')

# pull PID of participant w highest VAF
sm_in_ubs1_at_m41 %>% filter(AF > 0.3) %>% glimpse()
pid_w_highest_vaf = 2421575

# check whether any UBA1 M41L pts have CHIP
# system('gsutil cp gs://fc-secure-f1f8ca64-f52b-4690-a616-7bd09b94eb96/dsub/results/chip_annotation_and_filtering/caitlyn/batch2/CHIP_calls_batch2_nomyeloidcancers_06162023.txt .')
read_delim(file = 'CHIP_calls_batch2_nomyeloidcancers_06162023.txt') %>%
filter(person_id %in% sm_in_ubs1_at_m41$person_id) %>%
glimpse()
# this person has DNMT3A CHIP with VAF 20%
pid_w_chip = 1713303

cases_M41L = read_csv('M41L_paper/sm_in_uba1_at_m41.csv', show_col_types = FALSE)

cases_M41L %>%
group_by(age_at_biosample, dragen_sex_ploidy, ancestry_pred) %>%
summarise(
    num_cases = n(),
    num_controls_needed = 10*num_cases,
    case_pids = list(person_id)
) ->
control_demand_df

all_uba1_sm = read_csv('M41L_paper/sm_in_uba1_genetics_passing.csv', show_col_types = FALSE)

demog_plus %>%
filter(person_id %nin% all_uba1_sm$person_id) %>%
group_by(age_at_biosample, dragen_sex_ploidy, ancestry_pred) %>%
summarise(
    num_controls_avail = n(),
    avail_control_pids = list(person_id)
) ->
control_supply_df

set.seed(1)

options(width = 100)
control_demand_df %>%
left_join(control_supply_df) %>%
# mutate(deficit = num_controls_needed - num_controls_avail) %>% pull(deficit) %>% stem()
mutate(control_pids = map2(.x = avail_control_pids, .y = num_controls_needed, .f = sample, replace = FALSE)) %>%
pull(control_pids) %>% unlist() ->
control_pids

demog_plus %>%
filter(person_id %in% control_pids) %>%
glimpse() ->
controls_M41L

save_to_PD_and_GCS(controls_M41L)

cases_M41L = read_csv('M41L_paper/sm_in_uba1_at_m41.csv', show_col_types = FALSE)

demog_plus %>%
mutate(has_m41l = person_id %in% cases_M41L$person_id) %>%
glm(formula = has_m41l ~ age_at_biosample + sex_at_birth + race + ethnicity, family = 'binomial') %>%
summary()

demog_plus %>%
mutate(has_m41l = person_id %in% cases_M41L$person_id) %>%
glm(formula = has_m41l ~ age_at_biosample + dragen_sex_ploidy + ancestry_pred, family = 'binomial') %>%
summary()

demog_plus$age_at_biosample %>% range()

demog_plus %>%
CreateTableOne(
    data = .,
    vars = c('age_at_biosample', 'dragen_sex_ploidy', 'sex_at_birth', 'gender', 'ancestry_pred', 'race', 'ethnicity')
)

cases_M41L = read_csv('M41L_paper/sm_in_uba1_at_m41.csv', show_col_types = FALSE)
controls_M41L = read_csv('M41L_paper/controls_M41L.csv', show_col_types = FALSE)

demog_plus %>%
filter(person_id %in% c(cases_M41L$person_id, controls_M41L$person_id)) %>%
mutate(
    case_or_ctrl = case_when(
        (person_id %in% cases_M41L$person_id) & dragen_sex_ploidy == 'XY' ~ 'case_male',
        (person_id %in% cases_M41L$person_id) & dragen_sex_ploidy == 'XX' ~ 'case_female',
        (person_id %in% controls_M41L$person_id) & dragen_sex_ploidy == 'XY' ~ 'ctrl_male',
        (person_id %in% controls_M41L$person_id) & dragen_sex_ploidy == 'XX' ~ 'ctrl_female',
        .default = 'nomatch'
    )
) ->
demog_case_control

save_to_PD_and_GCS(demog_case_control)

demog_case_control %>%
mutate(wgs_case_ctrl = factor(case_or_ctrl, levels = c('ctrl_female', 'ctrl_male', 'case_female', 'case_male'))) %>%
# group_by(wgs_case_ctrl) %>% count()
CreateTableOne(
    data = .,
    strata = 'case_or_ctrl',
    vars = c('age_at_biosample', 'dragen_sex_ploidy', 'sex_at_birth', 'gender', 'ancestry_pred', 'race', 'ethnicity')
)

demog_case_control %>%
group_by(case_or_ctrl) %>%
summarise(youngest = min(age_at_biosample), oldest = max(age_at_biosample))

demog_plus %>%
filter(person_id == pid_w_highest_vaf)

demog_plus %>%
filter(person_id == pid_w_chip)

labs = read_csv('M41L_paper/labs_wgs_pts.csv', show_col_types = FALSE)

labs_case_control = labs_cases = filter(labs, person_id %in% c(cases_M41L$person_id, controls_M41L$person_id))

save_to_PD_and_GCS(labs_case_control)

demog_case_control = read_csv('M41L_paper/demog_case_control.csv')

labs_case_control %>%
mutate(
    lab_name = case_when(
        str_detect(standard_concept_name, regex('^leuko', ignore_case = TRUE)) ~ 'wbc',
        str_detect(standard_concept_name, regex('lympho', ignore_case = TRUE)) ~ 'lympho',
        str_detect(standard_concept_name, regex('platelet', ignore_case = TRUE)) ~ 'plt',
        str_detect(standard_concept_name, regex('hemoglobin', ignore_case = TRUE)) ~ 'hgb',
        str_detect(standard_concept_name, regex('sedimentation', ignore_case = TRUE)) ~ 'esr',
        str_detect(standard_concept_name, regex('reactive', ignore_case = TRUE)) ~ 'crp',
        str_detect(standard_concept_name, regex('MCV', ignore_case = TRUE)) ~ 'mcv',
        TRUE ~ 'nomatch'
    )
) %>%
# filter(lab_name == 'nomatch')
filter(value_as_number %nin% c(0, 1e7, 1e8)) %>%
filter(!is.na(value_as_number)) %>%
select(person_id, lab_name, lab_val = value_as_number,
       standard_concept_name, units = unit_concept_name, measurement_datetime) %>%
glimpse() ->
labs_w_names

labs_w_names %>%
filter(lab_name == 'wbc') %>%
mutate(units = ifelse(is.na(units), 'nounit', units)) %>%
mutate(lab_val = ifelse(units == 'per cubic millimeter', lab_val / 1000, lab_val)) %>%
# group_by(units) %>% summarise(n = n(), m = mean(lab_val, na.rm = TRUE))
filter(between(lab_val, 0.1, 50)) %>%
# ggplot(aes(x = lab_val, fill = units)) +
# geom_density()
identity() -> wbc

wbc %>% group_by(person_id) %>% summarise(min_wbc = min(lab_val)) %>%
left_join(demog_case_control) %>%
group_by(case_or_ctrl) %>%
summarise(
    num_pts = n(), 
    mean_min_wbc = round(mean(min_wbc), 1),
    sd_min_wbc = round(sd(min_wbc), 1)
)

labs_w_names %>%
filter(lab_name == 'lympho') %>%
mutate(units = ifelse(is.na(units), 'nounit', units)) %>%
filter(units %in% c('percent', 'No matching concept', 'no value', 'Percent', 'percent of white blood cells', 'nounits', 'Percentage unit')) %>%
# group_by(units) %>% summarise(n = n(), m = mean(lab_val, na.rm = TRUE)) %>% arrange(-n)
filter(between(lab_val, 0, 100)) %>%
# ggplot(aes(x = lab_val, fill = units)) +
# geom_density()
identity() -> lympho

lympho %>% group_by(person_id) %>% summarise(min_lympho = min(lab_val)) %>%
left_join(demog_case_control) %>%
group_by(case_or_ctrl) %>%
summarise(
    num_pts = n(), 
    mean_min_lympho = round(mean(min_lympho), 1),
    sd_min_lympho = round(sd(min_lympho), 1)
)

labs_w_names %>%
filter(lab_name == 'hgb') %>%
mutate(units = ifelse(is.na(units), 'nounit', units)) %>%
mutate(lab_val = ifelse(units == 'gram per liter', lab_val/10, lab_val)) %>%
# group_by(units) %>% summarise(n = n(), m = mean(lab_val, na.rm = TRUE)) %>% arrange(-n)
# ggplot(aes(x = lab_val, fill = units)) +
# geom_density()
identity() -> hgb

hgb %>% group_by(person_id) %>% summarise(min_hbg = min(lab_val)) %>%
left_join(demog_case_control) %>%
group_by(case_or_ctrl) %>%
summarise(
    num_pts = n(), 
    mean_min_hbg = round(mean(min_hbg), 1),
    sd_min_hbg = round(sd(min_hbg), 1)
)

labs_w_names %>%
filter(lab_name == 'mcv') %>%
mutate(units = ifelse(is.na(units), 'nounit', units)) %>%
# group_by(units) %>% summarise(n = n(), m = mean(lab_val, na.rm = TRUE)) %>% arrange(-n)
# ggplot(aes(x = lab_val, fill = units)) +
# geom_density()
identity() -> mcv

mcv %>% group_by(person_id) %>% summarise(max_mcv = max(lab_val)) %>%
left_join(demog_plus) %>%
mutate(
    wgs_case_ctrl = case_when(
        (person_id %in% cases_M41L$person_id) & dragen_sex_ploidy == 'XY' ~ 'case_male',
        (person_id %in% cases_M41L$person_id) & dragen_sex_ploidy == 'XX' ~ 'case_female',
        (person_id %in% controls_M41L$person_id) & dragen_sex_ploidy == 'XY' ~ 'ctrl_male',
        (person_id %in% controls_M41L$person_id) & dragen_sex_ploidy == 'XX' ~ 'ctrl_female',
        .default = 'nomatch'
    )
) %>%
group_by(wgs_case_ctrl) %>%
summarise(
    num_pts = n(), 
    mean_max_mcv = round(mean(max_mcv)),
    sd_max_mcv = round(sd(max_mcv))
)

labs_w_names %>%
filter(lab_name == 'esr') %>%
mutate(units = ifelse(is.na(units), 'nounit', units)) %>%
filter(units != 'mmol/h') %>%
# group_by(units) %>% summarise(n = n(), m = mean(lab_val, na.rm = TRUE)) %>% arrange(-n)
# ggplot(aes(x = lab_val, fill = units)) +
# geom_density()
identity() -> esr

esr %>% group_by(person_id) %>% summarise(max_esr = max(lab_val)) %>%
left_join(demog_case_control) %>%
group_by(case_or_ctrl) %>%
summarise(
    num_pts = n(), 
    mean_max_esr = round(mean(max_esr)),
    sd_max_esr = round(sd(max_esr))
)

labs_w_names %>%
filter(lab_name == 'crp') %>%
mutate(units = ifelse(is.na(units), 'nounit', units)) %>%
mutate(lab_val = ifelse(units == 'milligram per deciliter', lab_val * 10, lab_val)) %>%
filter(units != 'mg/L') %>%
# group_by(units) %>% summarise(n = n(), m = max(lab_val, na.rm = TRUE)) %>% arrange(-n)
# ggplot(aes(x = lab_val, fill = units)) +
# geom_density()
identity() -> crp

crp %>% group_by(person_id) %>% summarise(max_crp = max(lab_val)) %>%
left_join(demog_case_control) %>%
group_by(case_or_ctrl) %>%
summarise(
    num_pts = n(), 
    mean_max_crp = round(mean(max_crp)),
    sd_max_crp = round(sd(max_crp))
)

bind_rows(wbc, lympho, hgb, mcv, esr, crp) %>% filter(person_id == pid_w_highest_vaf)

bind_rows(wbc, lympho, hgb, mcv, esr, crp) %>% filter(person_id == pid_w_chip) %>%
ggplot(mapping = aes(x = as.Date(measurement_datetime), y = lab_val)) +
geom_point() +
facet_wrap(facets = vars(lab_name), scales = 'free')

dzs = read_csv('M41L_paper/dzs_wgs_pts.csv', show_col_types = FALSE)

glimpse(dzs)

paste(length(unique(dzs$person_id)), '(', round(100*length(unique(dzs$person_id))/nrow(demog_plus)), '%)')

dzs %>% mutate(y = year(as.Date(condition_start_datetime))) %>% 
group_by(person_id) %>% summarise(m = min(y)) %>% pull(m) -> year_of_first_snomed

paste(round(mean(2023 - year_of_first_snomed)), '(', round(sd(year_of_first_snomed)), ')')

dzs %>% select(person_id, standard_concept_name) %>% distinct() %>%
group_by(person_id) %>% count() %>% pull(n) -> num_snomed_codes

paste(round(mean(num_snomed_codes)), '(', round(sd(num_snomed_codes)), ')')
quantile(x = num_snomed_codes, probs = c(0.1, 0.5, 0.9))

males = bind_rows(cases_M41L, controls_M41L) %>% filter(dragen_sex_ploidy == 'XY')

dzs %>%
filter(person_id %in% c(cases_M41L$person_id, controls_M41L$person_id)) %>%
# filter(person_id %in% males$person_id) %>%
mutate(
    case_or_control = ifelse(person_id %in% cases_M41L$person_id, 'case', 'control'),
    snomed_name = standard_concept_name, 
    snomed_code = standard_concept_code
) %>%
select(-starts_with('standard')) %>%
glimpse() ->
dzs_case_and_control

save_to_PD_and_GCS(dzs_case_and_control)

demog_case_control = read_csv('M41L_paper/demog_case_control.csv')

dzs_case_and_control %>%
mutate(y = year(as.Date(condition_start_datetime))) %>%
group_by(person_id) %>% summarise(m = min(y), num_codes = length(unique(snomed_name))) %>% 
left_join(demog_case_control) %>%
group_by(case_or_ctrl) %>%
summarise(
    num_pts = n(), 
    emr_dur = round(mean(2023 - m)), emr_dur_sd = round(sd(m)),
    mean_num_codes = round(mean(num_codes)), sd_num_codes = round(sd(num_codes)),
    num_codes_q10 = quantile(num_codes, probs = 0.1),
    num_codes_q50 = quantile(num_codes, probs = 0.5),
    num_codes_q90 = quantile(num_codes, probs = 0.9),
)

dzs_case_and_control %>%
filter(person_id %in% pid_w_highest_vaf)

# interesting that she has langerhans cell histiocytosis, asthma,
# enthesopathy, tietze, pulm fibrosis, skin fibrosis
dzs_case_and_control %>%
filter(person_id %in% pid_w_chip) %>%
pull(snomed_name) %>% unique()

dzs_case_and_control %>%
select(
    person_id, case_or_control,
    snomed_name, snomed_code
) %>% 
distinct() %>% 
group_by(snomed_name) %>%
summarise(
    case_w = sum(case_or_control == 'case'),
    control_w = sum(case_or_control == 'control')) %>%
mutate(
    case_wo = nrow(cases_M41L) - case_w,
    control_wo = nrow(controls_M41L) - control_w,
) ->
snomed_scan

ms = list()
ors = list()
for (i in 1:nrow(snomed_scan)) {
    ms[[i]] = matrix(
        c(snomed_scan$control_wo[i], snomed_scan$control_w[i],
          snomed_scan$case_wo[i], snomed_scan$case_w[i]),
        nrow = 2, ncol = 2)
    ors[[i]] = fisher.test(ms[[i]])
}

snomed_scan$or = unlist(map(ors, 'estimate'))
snomed_scan$or_low = unlist(map(map(ors, 'conf.int'), function(x) { x[1] } ))
snomed_scan$or_hi = unlist(map(map(ors, 'conf.int'), function(x) { x[2] } ))
snomed_scan$p_val = unlist(map(ors, 'p.value'))

snomed_scan %>% arrange(p_val) %>% head(n = 30)

dzs_case_and_control %>%
select(
    person_id, case_or_control,
    snomed_name, snomed_code,
    starts_with('source_')
) %>% 
distinct() %>%
glimpse() ->
dzs_of_cases_and_controls_for_icd_scan

dzs_of_cases_and_controls_for_icd_scan %>%
group_by(source_vocabulary) %>% count() %>% arrange(-n)

dzs_of_cases_and_controls_for_icd_scan %>%
filter(source_vocabulary == 'ICD10CM') ->
dzs_icd10

read_csv('M41L_paper/icd10cmtoicd9gem.csv', show_col_types = FALSE) %>%
select(icd10_code = icd10cm, icd9_code = icd9cm) %>%
distinct(icd9_code, .keep_all = TRUE) %>%
# glimpse() %>%
identity() ->
crosswalk_icd9_icd10

dzs_of_cases_and_controls_for_icd_scan %>%
filter(source_vocabulary == 'ICD10CM') %>%
select(icd10_code = source_concept_code, icd10_name = source_concept_name) %>%
mutate(icd10_code = str_remove(string = icd10_code, pattern = '\\.')) %>%
distinct() %>%
glimpse() ->
icd10_names_and_codes

dzs_of_cases_and_controls_for_icd_scan %>%
filter(source_vocabulary == 'ICD9CM') %>%
mutate(icd9_dot_removed = str_remove(string = source_concept_code, pattern = '\\.')) %>%
inner_join(crosswalk_icd9_icd10, by = c('icd9_dot_removed' = 'icd9_code')) %>%
inner_join(icd10_names_and_codes) %>%
# glimpse() %>%
identity() ->
dzs_icd9_mapped_to_icd10

read_tsv('M41L_paper/tls_Icd10cmHumanReadableMap_US1000124_20230301.tsv', show_col_types = FALSE) %>% 
select(
    snomed_code = referencedComponentId, 
    snomed_name = referencedComponentName,
    icd_code = mapTarget, 
    icd_name = mapTargetName,
    mapGroup, mapPriority, mapRule, mapAdvice
) %>%
filter(mapGroup == 1, mapPriority == 1) %>%
select(starts_with('snomed'), starts_with('icd')) %>%
# glimpse() %>%
identity() ->
crosswalk_snomed_icd10

dzs_of_cases_and_controls_for_icd_scan %>%
filter(source_vocabulary == 'SNOMED') %>%
mutate(source_concept_code = as.numeric(source_concept_code)) %>%
inner_join(crosswalk_snomed_icd10) %>%
# glimpse() %>%
identity() ->
dzs_snomed_mapped_to_icd10

bind_rows(
    select(dzs_icd10, person_id, icd10_code = source_concept_code, icd10_name = source_concept_name),
    select(dzs_icd9_mapped_to_icd10, person_id, icd10_code, icd10_name),
    select(dzs_snomed_mapped_to_icd10, icd10_code = icd_code, icd10_name = icd_name)
) %>%
glimpse() %>%
mutate(case_or_control = ifelse(person_id %in% cases_M41L$person_id, 'case', 'control')) ->
icd10s_of_cases_and_controls

icd10s_of_cases_and_controls %>%
group_by(icd10_name, icd10_code) %>%
# group_by(icd10_name) %>% count() %>% arrange(-n)
summarise(
    case_w = sum(case_or_control == 'case'),
    control_w = sum(case_or_control == 'control')) %>%
mutate(
    case_wo = 74L - case_w,
    control_wo = 740L - control_w
) ->
icd_scan

ms = list()
ors = list()
for (i in 1:nrow(icd_scan)) {
    ms[[i]] = matrix(
        c(icd_scan$control_wo[i], icd_scan$control_w[i],
          icd_scan$case_wo[i], icd_scan$case_w[i]),
        nrow = 2, ncol = 2)
    ors[[i]] = fisher.test(ms[[i]])
}

icd_scan$or = unlist(map(ors, 'estimate'))
icd_scan$or_low = unlist(map(map(ors, 'conf.int'), function(x) { x[1] } ))
icd_scan$or_hi = unlist(map(map(ors, 'conf.int'), function(x) { x[2] } ))
icd_scan$p_val = unlist(map(ors, 'p.value'))

icd_scan %>%
arrange(p_val) %>%
head(n = 20)

icd10s_of_cases_and_controls %>%
mutate(icd10_base = word(string = icd10_code, sep = '\\.')) %>%
group_by(icd10_base) %>%
filter(!is.na(icd10_name)) %>%
summarise(
    commonest_icd_name = names(which.max(table(icd10_name))),
    case_w = sum(case_or_control == 'case'),
    control_w = sum(case_or_control == 'control')) %>%
mutate(
    case_wo = 74L - case_w,
    control_wo = 740L - control_w
) ->
icd_scan

ms = list()
ors = list()
for (i in 1:nrow(icd_scan)) {
    ms[[i]] = matrix(
        c(icd_scan$control_wo[i], icd_scan$control_w[i],
          icd_scan$case_wo[i], icd_scan$case_w[i]),
        nrow = 2, ncol = 2)
    ors[[i]] = fisher.test(ms[[i]])
}

icd_scan$or = unlist(map(ors, 'estimate'))
icd_scan$or_low = unlist(map(map(ors, 'conf.int'), function(x) { x[1] } ))
icd_scan$or_hi = unlist(map(map(ors, 'conf.int'), function(x) { x[2] } ))
icd_scan$p_val = unlist(map(ors, 'p.value'))

icd_scan %>%
arrange(p_val) %>%
head(n = 20)

crosswalk_icd10_phecode = read_csv('M41L_paper/Phecode_map_v1_2_icd10cm_beta.csv') %>% glimpse()

icd10s_of_cases_and_controls %>% 
inner_join(crosswalk_icd10_phecode, by = c('icd10_code' = 'icd10cm'), relationship = 'many-to-many') %>%
group_by(phecode_str) %>%
distinct(person_id, phecode, phecode_str, case_or_control) %>%
summarise(
    case_w = sum(case_or_control == 'case'),
    control_w = sum(case_or_control == 'control')) %>%
mutate(
    case_wo = 74L - case_w ,
    control_wo = 740L - control_w
) ->
phecode_scan

ms = list()
ors = list()
for (i in 1:nrow(phecode_scan)) {
    ms[[i]] = matrix(
        c(phecode_scan$control_wo[i], phecode_scan$control_w[i],
          phecode_scan$case_wo[i], phecode_scan$case_w[i]),
        nrow = 2, ncol = 2)
    ors[[i]] = fisher.test(ms[[i]])
}

phecode_scan$or = unlist(map(ors, 'estimate'))
phecode_scan$or_low = unlist(map(map(ors, 'conf.int'), function(x) { x[1] } ))
phecode_scan$or_hi = unlist(map(map(ors, 'conf.int'), function(x) { x[2] } ))
phecode_scan$p_val = unlist(map(ors, 'p.value'))

phecode_scan %>%
arrange(p_val) %>%
head(n = 20)

regexci = function(pattern) regex(pattern = pattern, ignore_case = TRUE)

dzs_case_and_control %>%
distinct(person_id, case_or_control, snomed_name, snomed_code) %>%
mutate(vexas_man = case_when(
    str_detect(snomed_name,
        regexci('iridocyclitis|(inflammatory eye)|(chondritis)|(sensorineural hearing)|cryopyrin|fever|(periorbital edema)')
    ) ~ 'head and neck', 
    str_detect(snomed_name,
        regexci('alveolitis|(pleural effusion)|(myocarditis)|(pulmonary infiltrates)')
    ) ~ 'thorax',
    str_detect(snomed_name,
        regexci('myelodysplastic|myeloma|anemia|leukopenia|thrombocytopenia')
    ) ~ 'marrow',
    str_detect(snomed_name,
        regexci('splenomegaly|colitis')
    ) ~ 'abdomen',
    str_detect(snomed_name,
        regexci('(rheumatoid arthritis)|lupus|sjogren|scleroderma|myositis|(psoriatic arthritis)|ankylosing|behcet')
    ) ~ 'inflammatory\narthritis',
    str_detect(snomed_name,
        regexci('epididymitis|orchitis')
    ) ~ 'pelvis',
    str_detect(snomed_name,
        regexci('arteritis|angiitis|vasculitis|aortitis|takayasu|(giant cell arteritis)|kawasaki')
    ) ~ 'vasculitis',
    str_detect(snomed_name,
        regexci('dermatosis|urticaria|psoriasis|sweet')
    ) ~ 'cutaneous',
    str_detect(snomed_name,
        regexci('thrombosis|(pulmonary embolism)|stroke'),
    ) ~ 'thrombosis',    
    TRUE ~ 'nomatch')
) %>%
# filter(vexas_man == 'cutaneous') %>%
# group_by(snomed_name) %>% count() %>% arrange(-n)
filter(vexas_man != 'nomatch') %>%
glimpse() ->
dzs_of_cases_and_controls_w_vexas_man

dzs_of_cases_and_controls_w_vexas_man %>%
group_by(vexas_man, snomed_name, snomed_code) %>%
count() %>% arrange(vexas_man, -n) %>%
write_csv('M41L_paper/snomed_names_and_codes_of_vexas_manifestations.csv')

dzs_of_cases_and_controls_w_vexas_man %>%
group_by(vexas_man) %>%
summarise(
    case_w = sum(case_or_control == 'case'),
    control_w = sum(case_or_control == 'control')) %>%
mutate(
    case_wo = 74L - case_w,
    control_wo = 740L - control_w
) ->
vexas_scan

ms = list()
ors = list()
for (i in 1:nrow(vexas_scan)) {
    ms[[i]] = matrix(
        c(vexas_scan$control_wo[i], vexas_scan$control_w[i],
          vexas_scan$case_wo[i], vexas_scan$case_w[i]),
        nrow = 2, ncol = 2)
    ors[[i]] = fisher.test(ms[[i]])
}

vexas_scan$OR = unlist(map(ors, 'estimate'))
vexas_scan$or_low = unlist(map(map(ors, 'conf.int'), function(x) { x[1] } ))
vexas_scan$or_hi = unlist(map(map(ors, 'conf.int'), function(x) { x[2] } ))
vexas_scan$p_val = unlist(map(ors, 'p.value'))

vexas_scan %>%
arrange(p_val) %>%
head(n = 10)

vexas_scan %>%
# to remove very poorly estimated
filter(or_hi < 50) %>%
# filter(case_w >= 2, control_w >= 2) %>%
# filter(p_val < 0.01) %>%
# mutate(icd_name = factor(icd10_name, levels = icd10_name[order(or)])) %>%
ggplot(aes(y = vexas_man, x = OR, xmin = or_low, xmax = or_hi)) +
geom_point(size = 3) +
geom_errorbarh(height = 0.2) +
geom_text(aes(label = paste0(format(round(OR, 1), nsmall = 1),
                             ', [', 
                             round(or_low, 1),
                             ', ', 
                             round(or_hi, 1), 
                             '], n = ',
                             case_w + control_w)),
         x = 1.2, nudge_y = 0.25, hjust = 0) +
scale_x_continuous(name = NULL) +
scale_y_discrete(name = 'odds ratio') +
geom_vline(xintercept = 1) +
theme_bw()

# ggsave(filename = 'M41L_paper/vexas_man_scan_forest_plot.pdf', height = 6, width = 6)

meds = read_csv('M41L_paper/meds_filtered.csv')

meds %>%
filter(person_id %in% c(cases_M41L$person_id, controls_M41L$person_id)) %>%
mutate(case_or_control = ifelse(person_id %in% cases_M41L$person_id, 'case', 'control')) %>%
glimpse() ->
meds_of_cases_and_controls

meds_of_cases_and_controls  %>%
mutate(
    med_type = case_when(
        str_detect(standard_concept_name, 
                   regex('predni|dexamethasone', ignore_case = TRUE)) ~ 'systemic steroid',
        str_detect(standard_concept_name, 
                   regex('cortisone|fluocinolone|budesonide|desonide|fluocinonide|triamcinolone|hydrocortisone|clobetasol|betamethasone', ignore_case = TRUE)) ~ 'topical steroid',
        str_detect(standard_concept_name, 
                   regex('beclomethasone|mometasone', ignore_case = TRUE)) ~ 'inhaled steroid',
        str_detect(standard_concept_name, 
                   regex('tacro|cyclosporine', ignore_case = TRUE)) ~ 'CNI',
        str_detect(standard_concept_name, 
                   regex('methotrexate|leflunomide|azathioprine|mycophen', ignore_case = TRUE)) ~ 'cDMARD',
        str_detect(standard_concept_name, 
                   regex('infliximab|adalimumab|golimumab|etanercept|tocilizumab|guselkumab|ustekinumab|abatacept|apremilast', ignore_case = TRUE)) ~ 'bDMARD',
        str_detect(standard_concept_name, 
                   regex('tofacitinib', ignore_case = TRUE)) ~ 'tsDMARD',
        TRUE ~ 'nomatch'
    )
) %>%
distinct(person_id, med_type, .keep_all = TRUE) %>%
filter(med_type != 'nomatch') %>%
glimpse() ->
meds_of_cases_and_controls_w_med_type

meds_of_cases_and_controls_w_med_type %>%
group_by(med_type) %>%
summarise(
    commonest_med = names(sort(table(standard_concept_name), decreasing = TRUE))[1],
    case_w = sum(case_or_control == 'case'),
    control_w = sum(case_or_control == 'control')) %>%
mutate(
    case_wo = 74L - case_w,
    control_wo = 740L - control_w
) ->
case_control_df_meds

ms = list()
ors = list()
for (i in 1:nrow(case_control_df_meds)) {
    ms[[i]] = matrix(
        c(case_control_df_meds$control_wo[i], case_control_df_meds$control_w[i],
          case_control_df_meds$case_wo[i], case_control_df_meds$case_w[i]),
        nrow = 2, ncol = 2)
    ors[[i]] = fisher.test(ms[[i]])
}

case_control_df_meds$or = unlist(map(ors, 'estimate'))
case_control_df_meds$or_low = unlist(map(map(ors, 'conf.int'), function(x) { x[1] } ))
case_control_df_meds$or_hi = unlist(map(map(ors, 'conf.int'), function(x) { x[2] } ))
case_control_df_meds$p_val = unlist(map(ors, 'p.value'))

case_control_df_meds %>%
arrange(p_val)

# case_control_df_meds %>%
# mutate(or_hi = ifelse(or_hi > 10, 10, or_hi)) %>%
# mutate(med_type = factor(
#     med_type,
#     levels = med_type[order(or)],
#     labels = c('inhaled steroid', 'bDMARD', 'csDMARD', 'topical steroid', 'CNI', 'systemic steroid')
# )) ->
# # glimpse() ->
# toplot

toplot %>%
filter(or_hi < 10) %>%
ggplot(aes(y = med_type, x = or, xmin = or_low, xmax = or_hi)) +
geom_point(size = 3) +
geom_errorbarh(height = 0.2) +
geom_text(aes(label = paste0(format(round(or, 1), nsmall = 1),
                             ', [', 
                             round(or_low, 1),
                             ', ', 
                             round(or_hi, 1), 
                             '], n = ',
                             case_w + control_w)),
         x = 1.2, nudge_y = 0.25, hjust = 0) +
scale_x_continuous(name = NULL) +
scale_y_discrete(name = 'odds ratio') +
geom_vline(xintercept = 1) +
theme_bw()

ggsave(filename = 'M41L_paper/vexas_meds_forest_plot.pdf', height = 3.5, width = 6)

meds_of_cases_and_controls_w_med_type %>%
filter(person_id %in% pid_w_highest_vaf)

meds_of_cases_and_controls_w_med_type %>%
filter(person_id %in% pid_w_chip)

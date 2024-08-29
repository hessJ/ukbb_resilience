require(plyr)
require(data.table)
require(dplyr)
 
# load UK BB compiled data
sub = readRDS("~/Google Drive/mac_storage/Manuscripts/Resilience/scz/scz-bull_special-issue/draft/manuscript_docs/circulate/edits/final_to_submit/revision/ukbb_compiled_data_qc_filtered-w-comorbidity-score.Rdata")
 

sub$dx = sub$ICD10_F20.Schizophrenia
sub_id = sub[!is.na(sub$dx), ]

fit = glm(dx ~ risk*resilience + comorbidity_score + GenomeWidePC_1 + GenomeWidePC_2 + 
            GenomeWidePC_3 + GenomeWidePC_4 + GenomeWidePC_5 +
            f.21022.0.0_Age.at.recruitment +  f.31.0.0_Sex, data = sub_id, family=binomial())

stats = data.frame(summary(fit)$coefficients)
stats = data.frame(pred = rownames(stats), stats, ntotal = nrow(sub_id), ncase = sum(sub_id$dx == 1))
stats$or = exp(stats$Estimate)
stats$lwr = exp(stats$Estimate - 1.96 * stats$Std..Error)
stats$upr = exp(stats$Estimate + 1.96 * stats$Std..Error)
stats
 
icd10_stats_f20 = stats

icd10_stats_f20 = data.frame(response = "ICD10_F20.Schizophrenia",
                                 pred = icd10_stats_f20$pred,
                                 estimate = icd10_stats_f20$Estimate,
                                 se = icd10_stats_f20$Std..Error,
                                 pval = icd10_stats_f20$Pr...z..,
                                 n_sample =icd10_stats_f20$ntotal,
                                n_case = icd10_stats_f20$ncase)
 
# repeat lifetime ICD code analysis using Elixhauser Comorbidity Index as additional covariates
# ids = ids[!grepl("F20", ids)]
# stats_list = list();
# plot_save = list()
# for(x in 1:length(ids)){
#   cat("\rVariable...", x, "of",length(ids))
#   # // non schizophrenia cases
#   dx = sub[,colnames(sub) %in% ids[[x]]]
#   sub$dx = dx
#   sub_id = sub[!is.na(sub$dx), ]
#   sub_id = sub_id[sub_id$ICD10_F20.Schizophrenia == 0, ]
#   
#   fit = glm(dx ~ risk*resilience + GenomeWidePC_1 + GenomeWidePC_2 + 
#               GenomeWidePC_3 + GenomeWidePC_4 + GenomeWidePC_5 +
#               f.21022.0.0_Age.at.recruitment +  f.31.0.0_Sex, data = sub_id, family=binomial())
#   
#   stats = data.frame(summary(fit)$coefficients)
#   stats = data.frame(pred = rownames(stats), stats, ntotal = nrow(sub_id), ncase = sum(dx == 1))
#   stats$response = ids[[x]]
#   stats_list[[x]] = stats
#   plot_save[[x]] = interactions::interact_plot(fit, pred='risk', modx = 'resilience')$data
# }
# res = ldply(stats_list)
# res = ddply(res, .(pred), transform ,fdr = p.adjust(Pr...z.., 'fdr'))
# all_res_narrow = res
# res = res[grepl("risk:resilience", res$pred), ]
# 
# sig = res[res$fdr < .05, ]
 
# -- include all subjects, incuding those with scz diagnosis for icd narrow analysis
# repeat lifetime ICD code analysis using Elixhauser Comorbidity Index as additional covariates
counts = colSums(sub[,grepl("ICD10_F|lupus|crohn|diabetes|asthma|Fibrosis.and.cirrhosis.of.liver", colnames(sub))])
ids = names(counts[counts >= 500])

ids = ids[!grepl("F20", ids)]
stats_list = list();
plot_save = list()
for(x in 1:length(ids)){
  cat("\rVariable...", x, "of",length(ids))
  # // non schizophrenia cases
  dx = sub[,colnames(sub) %in% ids[[x]]]
  sub$dx = dx
  sub_id = sub[!is.na(sub$dx), ]
  # sub_id = sub_id[sub_id$ICD10_F20.Schizophrenia == 0, ]
  
  fit = glm(dx ~ risk*resilience +  comorbidity_score + ICD10_F20.Schizophrenia + GenomeWidePC_1 + GenomeWidePC_2 + 
              GenomeWidePC_3 + GenomeWidePC_4 + GenomeWidePC_5 +
              f.21022.0.0_Age.at.recruitment +  f.31.0.0_Sex, data = sub_id, family=binomial())
  
  stats = data.frame(summary(fit)$coefficients)
  stats = data.frame(pred = rownames(stats), stats, ntotal = nrow(sub_id), ncase = sum(dx == 1))
  stats$response = ids[[x]]
  stats_list[[x]] = stats
  plot_save[[x]] = interactions::interact_plot(fit, pred='risk', modx = 'resilience')$data
}
res = ldply(stats_list)
res = ddply(res, .(pred), transform ,fdr = p.adjust(Pr...z.., 'fdr'))
ddply(res, .(pred), summarize, sum(fdr < .05))

all_res_narrow = res
 
all_res_narrow = data.frame(   response = all_res_narrow$response, 
                                 pred = all_res_narrow$pred,
                                 estimate = all_res_narrow$Estimate,
                                 se = all_res_narrow$Std..Error,
                                 pval = all_res_narrow$Pr...z..,
                                 n_sample = all_res_narrow$ntotal,
                                 n_case = all_res_narrow$ncase)
 


# == lifetime diagnoses in icd-10 chapters
counts = colnames(sub)[grepl("ICD10", colnames(sub))]
icd_no = gsub("ICD10_", "", unique(unlist(lapply(strsplit(counts, "[.]"), function(x) x[[1]]))))
cat_grab = unique(gsub("\\d+", "", icd_no)) # remove numbers from character string

# plot_int_combined = plot_int_combined[plot_int_combined$ICD10_F20.Schizophrenia == 0, ]

coef_list = list();
plot_list = list();
coef_main_list = list()
for(x in 1:length(cat_grab)){
  cat("\rAnalyzing category:", x)
  id_grab = paste0(paste0(paste0(cat_grab[[x]],0, 0:9), collapse="|"), paste0(paste0(cat_grab[[x]], 10:93), collapse="|"), collapse="|")
  disorders_icd = sub[,grepl(id_grab, colnames(sub))]
  disorders_icd = rowSums(disorders_icd)
  disorders_icd[disorders_icd > 1]  = 1
  table(disorders_icd)
  
  fit = glm(disorders_icd ~ resilience*risk +  comorbidity_score + ICD10_F20.Schizophrenia + GenomeWidePC_1 + GenomeWidePC_2 + GenomeWidePC_3 + GenomeWidePC_4 + GenomeWidePC_5 +
              f.31.0.0_Sex + f.21022.0.0_Age.at.recruitment, data = sub, family=binomial())
  
  coef = summary(fit)$coefficients
  coef = data.frame(pred = rownames(coef), coef)
  coef$n_case = sum(disorders_icd == 1)
  coef$n_control = sum(disorders_icd == 0)
  coef$n_total = length(disorders_icd)
  coef_list[[x]] = coef
  
  out = interactions::interact_plot(fit, pred = 'risk', modx = 'resilience')
  out = data.frame(out$data)
  plot_list[[x]] = out
  
  fit = glm(disorders_icd ~ resilience + risk + comorbidity_score + ICD10_F20.Schizophrenia + GenomeWidePC_1 + GenomeWidePC_2 + GenomeWidePC_3 + GenomeWidePC_4 + GenomeWidePC_5 +
              f.31.0.0_Sex + f.21022.0.0_Age.at.recruitment, data = sub, family=binomial())
  
  stats = summary(fit)$coefficients
  stats = data.frame(pred = rownames(stats), stats, n_case = sum(disorders_icd == 1), n_control = length(disorders_icd), n_total = nrow(sub))
  stats$name = cat_grab[[x]]
  coef_main_list[[x]] = stats
  
  
}

names(coef_list) = cat_grab
res = ldply(coef_list, .id='response')
res = ddply(res, .(pred), transform, fdr = p.adjust(Pr...z.., 'fdr'))
ddply(res, .(pred), summarize, sum(fdr < .05))
all_icd_broad = res

chapters = data.frame(response = cat_grab,
                      chapter = c("Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism",
                                  "Endocrine, nutritional and metabolic diseases",
                                  "Mental and behavioural disorders",
                                  "Diseases of the nervous system",
                                  "Diseases of the circulatory system",
                                  "Diseases of the respiratory system",
                                  "Diseases of the digestive system",
                                  "Diseases of the skin and subcutaneous tissue",
                                  "Diseases of the musculoskeletal system and connective tissue",
                                  "Diseases of the genitourinary system",
                                  "Symptoms, signs and abnormal clinical and laboratory findings, not elsewhere classified",
                                  "Injury, poisoning and certain other consequences of external causes (S)",
                                  "Injury, poisoning and certain other consequences of external causes (T)",
                                  "External causes of morbidity and mortality (V)",
                                  "External causes of morbidity and mortality (W)",
                                  "External causes of morbidity and mortality (X)",
                                  "External causes of morbidity and mortality (Y)"))

 
all_icd_broad = merge(chapters, all_icd_broad, by='response')
all_icd_broad$chapter = gsub("[,]", "/", all_icd_broad$chapter)
  

all_icd_broad = data.frame(   response = all_icd_broad$chapter, 
                              pred = all_icd_broad$pred,
                              estimate = all_icd_broad$Estimate,
                              se = all_icd_broad$Std..Error,
                              pval = all_icd_broad$Pr...z..,
                              n_sample = all_icd_broad$n_total,
                              n_case = all_icd_broad$n_case)

## education
grab = sub[!is.na(sub$isced_coding),]

fit2 = lm(RNOmni::RankNorm(isced_coding) ~ ICD10_F20.Schizophrenia + f.21022.0.0_Age.at.recruitment + f.31.0.0_Sex +  risk*resilience +  GenomeWidePC_1 +
            GenomeWidePC_2 + GenomeWidePC_3 + GenomeWidePC_4 + GenomeWidePC_5, data = grab)
coef_educ = data.frame(summary(fit2)$coefficients)
coef_educ$n_sample = nrow(grab)
coef_educ
 
## fluid intelligence
grab = sub[!is.na(sub$f.20016.0.0_Fluid.intelligence.score), ]

fit = lm(scale(f.20016.0.0_Fluid.intelligence.score) ~ ICD10_F20.Schizophrenia + f.21022.0.0_Age.at.recruitment + f.31.0.0_Sex +  risk*resilience +  GenomeWidePC_1 +
           GenomeWidePC_2 + GenomeWidePC_3 + GenomeWidePC_4 + GenomeWidePC_5, data = grab)
coef_fluid = data.frame(summary(fit)$coefficients)
coef_fluid$n_sample = nrow(grab)
coef_fluid
 

coef_fluid = data.frame(pred = rownames(coef_fluid), coef_fluid)
coef_educ = data.frame(pred = rownames(coef_educ), coef_educ)
all = ldply(list(intel = coef_fluid, educ = coef_educ))
all = ddply(all, .(pred), transform, fdr = p.adjust(Pr...t.., 'fdr'))
 
cognitive_stats = data.frame(   response = all$.id, 
                              pred = all$pred,
                              estimate = all$Estimate,
                              se = all$Std..Error,
                              pval = all$Pr...t..,
                              n_sample = all$n_sample)

all_stats_in_script = ldply(list(scz = icd10_stats_f20,
           icd_narrow = all_res_narrow,
           icd_broad = all_icd_broad,
           cognitive = cognitive_stats), .id='category')
all_stats_in_script = all_stats_in_script[grepl("risk|resilience", all_stats_in_script$pred), ]
all_stats_in_script = ddply(all_stats_in_script, .(pred, category), transform, fdr = p.adjust(pval, 'fdr'))
ddply(all_stats_in_script, .(category, pred), summarize, sum(fdr < .05))

fwrite(all_stats_in_script, file="~/Google Drive/mac_storage/Manuscripts/Resilience/scz/scz-bull_special-issue/draft/manuscript_docs/circulate/edits/final_to_submit/revision/stats_compile_ukbb-1_resilience_0.2.txt", quote=F,row.names=F, sep="\t")



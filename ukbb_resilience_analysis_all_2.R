
# mental health items
require(data.table);
require(plyr);
require(dplyr)

# load ukbb compiled data 
ukbb = readRDS("~/Google Drive/mac_storage/Manuscripts/Resilience/scz/scz-bull_special-issue/draft/manuscript_docs/circulate/edits/final_to_submit/revision/ukbb_compiled_data_qc_filtered-w-comorbidity-score.Rdata")
ukbb = data.frame(ukbb)

showcase = data.frame(fread("/Users/jonathanhess/Documents/ukbiobank/phenotypes/Data_Dictionary_Showcase.tsv"))
showcase = showcase[grepl("mental", ignore.case = T, showcase$Path), ]
showcase = showcase[grepl("Psychosocial factors", showcase$Path),]
showcase = showcase[showcase$ValueType %in% "Categorical single",]
showcase = showcase[grepl("ACE", showcase$Notes), ]
showcase = showcase[showcase$Participants >= 400e3, ]
range(showcase$Participants)

showcase$iter = paste0("f.",showcase$FieldID,".0.0")
showcase = showcase[showcase$iter %in% colnames(ukbb), ]
iter = showcase$iter


codings = data.frame(fread("~/Documents/ukbiobank/phenotypes/Codings.tsv", fill = T))
codings = codings[codings$Coding %in% showcase$Coding, ]


plot_list = list();
coef_list = list();
main_coef_list = list();
for(x in 1:length(iter)){
  cat("\rProgress:", x)
  grab = ukbb[,colnames(ukbb) %in% iter[[x]]]
  ukbb$outcome = grab
  grab = ukbb[ukbb$outcome >= 0, ]
  grab = grab[!is.na(grab$outcome), ]
  levels = length(unique(grab$outcome))
  if(levels > 2) next
  
  
  fit  = glm(outcome ~ resilience*risk + comorbidity_score + f.21022.0.0_Age.at.recruitment + f.31.0.0_Sex +  ICD10_F20.Schizophrenia + GenomeWidePC_1 + GenomeWidePC_2 + GenomeWidePC_3 + GenomeWidePC_4 + GenomeWidePC_5, data = grab, family=binomial())
  coef = summary(fit)$coefficients
  coef = data.frame(pred = rownames(coef), coef)
  coef$n_sample = nrow(grab)
  coef$outcome = iter[[x]]
  coef$n_case = sum(grab$outcome)
  coef_list[[x]] = coef

  plot_list[[x]] = data.frame(interactions::interact_plot(fit, pred='risk', modx = 'resilience', data = grab)$data)
  
  rm(grab)
  

  
}

names(plot_list) = iter
int_plot_mental_health = ldply(plot_list, .id='response')

showcase = data.frame(fread("/Users/jonathanhess/Documents/ukbiobank/phenotypes/Data_Dictionary_Showcase.tsv"))
showcase = showcase[grepl("mental", ignore.case = T, showcase$Path), ]
showcase = showcase[grepl("Psychosocial factors", showcase$Path),]
showcase = showcase[showcase$ValueType %in% "Categorical single",]
showcase = showcase[grepl("ACE", showcase$Notes), ]
showcase$FieldID = paste0("f.",showcase$FieldID,".0.0")
showcase = showcase[showcase$FieldID %in% int_plot_mental_health$response, ]
showcase = data.frame(response = showcase$FieldID,
                      id = showcase$Field)

int_plot_mental_health = merge(showcase, int_plot_mental_health, by='response')
int_plot_mental_health = int_plot_mental_health[,!colnames(int_plot_mental_health) %in% "response"]
colnames(int_plot_mental_health)[colnames(int_plot_mental_health) %in% "id"] = "response"


res = ldply(coef_list)
res = ddply(res, .(pred), transform, fdr = p.adjust(Pr...z.., 'fdr'))
ddply(res, .(pred), summarize, sum(fdr < .05))

range(res$n_sample) # number of participants that completed the item in the questionarie

colnames(res)[colnames(res) %in% "outcome"] = "id"

showcase = data.frame(fread("/Users/jonathanhess/Documents/ukbiobank/phenotypes/Data_Dictionary_Showcase.tsv"))
showcase = showcase[grepl("mental", ignore.case = T, showcase$Path), ]
showcase = showcase[grepl("Psychosocial factors", showcase$Path),]
showcase = showcase[showcase$ValueType %in% "Categorical single",]
showcase = showcase[grepl("ACE", showcase$Notes), ]
label = data.frame(id = paste0("f.",showcase$FieldID,".0.0"), Field = showcase$Field)
res = merge(label, res, by='id')


mental_health_stats = res

mental_health_stats = data.frame(fieldid = mental_health_stats$id, 
                                 response = mental_health_stats$Field, 
                                 pred = mental_health_stats$pred,
                                 estimate = mental_health_stats$Estimate,
                                 se = mental_health_stats$Std..Error,
                                 pval = mental_health_stats$Pr...z..,
                                 n_sample = mental_health_stats$n_sample,
                                 n_case = mental_health_stats$n_case)




# self-harm or contemplated 
sub = ukbb[ukbb$ever_contemplated_self_harm >= 0, ]
sub$outcome = sub$ever_contemplated_self_harm
sub$outcome[sub$outcome > 1] = 1
sub = sub[!is.na(sub$outcome), ]
fit = glm(outcome ~ ICD10_F20.Schizophrenia + comorbidity_score + f.21022.0.0_Age.at.recruitment + f.31.0.0_Sex +  risk*resilience +  GenomeWidePC_1 +
            GenomeWidePC_2 + GenomeWidePC_3 + GenomeWidePC_4 + GenomeWidePC_5, data = sub, family=binomial())
coef = summary(fit)$coefficients
coef_ideation = data.frame(pred = rownames(coef), coef)
coef_ideation$or = exp(coef_ideation$Estimate)
coef_ideation$lwr = round(exp(coef_ideation$Estimate - 1.96 * coef_ideation$Std..Error),3)
coef_ideation$upr = round(exp(coef_ideation$Estimate + 1.96 * coef_ideation$Std..Error),3)
coef_ideation$n_sample = nrow(sub)
coef_ideation$n_case = sum(sub$outcome)

int_plot_ideation = data.frame(interactions::interact_plot(fit, pred='risk', modx='resilience', data = sub)$data)
int_plot_ideation$response = "ever_contemplated_self_harm"

sub = ukbb[ukbb$ever_self_harmed >= 0, ]
sub$outcome = sub$ever_self_harmed
sub$outcome[sub$outcome > 1] = 1
sub = sub[!is.na(sub$outcome), ]
fit = glm(outcome ~ ICD10_F20.Schizophrenia + comorbidity_score + f.21022.0.0_Age.at.recruitment + f.31.0.0_Sex +  risk*resilience +  GenomeWidePC_1 +
            GenomeWidePC_2 + GenomeWidePC_3 + GenomeWidePC_4 + GenomeWidePC_5, data = sub, family=binomial())
coef = summary(fit)$coefficients
coef_harm = data.frame(pred = rownames(coef), coef)
coef_harm$or = exp(coef_harm$Estimate)
coef_harm$lwr = round(exp(coef_harm$Estimate - 1.96 * coef_harm$Std..Error),3)
coef_harm$upr = round(exp(coef_harm$Estimate + 1.96 * coef_harm$Std..Error),3)
coef_harm$n_sample = nrow(sub)
coef_harm$n_case = sum(sub$outcome)

int_plot_harm = data.frame(interactions::interact_plot(fit, pred="risk", modx="resilience", data = sub)$data)
int_plot_harm$response = "ever_engaged_in_self_harm"

res = ldply(list(`contemplated` = coef_ideation, `selfharm` = coef_harm), .id='response')
res = ddply(res, .(pred), transform, fdr = p.adjust(Pr...z.., 'fdr'))

harm_res = res

harm_res = data.frame(response = harm_res$response,
                     pred = harm_res$pred,
                     estimate = harm_res$Estimate,
                     se = harm_res$Std..Error,
                     pval = harm_res$Pr...z..,
                     n_sample = harm_res$n_sample)


# is self-reported ethnicity associated with ICD-10 code for schizophrenia (F20)?
# fit = glm(ICD10_F20.Schizophrenia ~  factor(f.21000.0.0) +  f.21022.0.0_Age.at.recruitment +  f.31.0.0_Sex, data = sub, family=binomial())
# car::Anova(fit, type='III')


# == run analysis of smri
ukb_smri = data.frame(fread("~/Documents/ukbiobank/ukb_smri_data.txt", nrows=1))
iters = colnames(ukb_smri)[-1]
iters = iters[grepl("\\.2.0", iters)]
iters = iters[iters %in% colnames(ukbb)]

grab = ukbb[,colnames(ukbb) %in% c(iters, "ICD10_F20.Schizophrenia", "GenomeWidePC_1", "GenomeWidePC_2","GenomeWidePC_3", "GenomeWidePC_4","GenomeWidePC_5","resilience","risk","f.25000.2.0",
                                   "f.25756.2.0", "f.25757.2.0", "f.25758.2.0", "f.25759.2.0", "f.21002.2.0", "f.21001.2.0", "f.12144.2.0", "f.31.0.0_Sex",
                                   "f.21003.2.0_Age.when.attended.assessment.centre")]
grab = grab[!rowSums(is.na(grab))/ncol(grab) > 0.5, ]
grab = grab[,!colSums(is.na(grab)/nrow(grab) > 0.5)]

grab$f.31.0.0_Sex = factor(grab$f.31.0.0_Sex)

coef_list = list();
main_coef_list = list();
plot_list = list();
for(x in 1:length(iters)){
  cat("\rProgress:",x)
  response = grab[,colnames(grab) %in% iters[[x]]]
  grab$outcome = scale(response)
  
  fit = lm(outcome ~ resilience * risk  + GenomeWidePC_1 + GenomeWidePC_2 + GenomeWidePC_3 + GenomeWidePC_4 + GenomeWidePC_5 + (f.31.0.0_Sex) + f.21003.2.0_Age.when.attended.assessment.centre + 
             f.25000.2.0 + f.25756.2.0 + f.25757.2.0 + f.25758.2.0 + f.25759.2.0 + f.21002.2.0 + f.21001.2.0 + f.12144.2.0, data = grab)
  
  coef = summary(fit)$coefficients
  coef = data.frame(pred = rownames(coef), coef)
  coef$n_total = nrow(grab)
  coef$response = iters[[x]]
  coef_list[[x]] = coef
  
  plot_list[[x]] = data.frame(interactions::interact_plot(fit, pred="risk", modx="resilience", data = grab)$data)
  plot_list[[x]]$response = iters[[x]]
  
  # main effects w/o interaction
  fit = lm(outcome ~ resilience + risk   + GenomeWidePC_1 + GenomeWidePC_2 + GenomeWidePC_3 + GenomeWidePC_4 + GenomeWidePC_5 +  (f.31.0.0_Sex) + f.21003.2.0_Age.when.attended.assessment.centre + 
             f.25000.2.0 + f.25756.2.0 + f.25757.2.0 + f.25758.2.0 + f.25759.2.0 + f.21002.2.0 + f.21001.2.0 + f.12144.2.0, data = grab)
  
  coef = summary(fit)$coefficients
  coef = data.frame(pred = rownames(coef), coef)
  coef$n_total = nrow(grab)
  coef$response = iters[[x]]
  main_coef_list[[x]] = coef
}

# summary of models w/ interaction of risk*resilience
res = ldply(coef_list)
res = ddply(res, .(pred), transform, fdr = p.adjust(Pr...t.., 'fdr'))
ddply(res, .(pred), summarize, sum(fdr < .05))
ddply(res, .(pred), summarize, sum(Pr...t.. < .05))


int_plot_smri = ldply(plot_list)


# label parcellations
showcase = data.frame(fread("/Users/jonathanhess/Documents/ukbiobank/phenotypes/Data_Dictionary_Showcase.tsv"))
showcase$FieldID = paste0("f.",showcase$FieldID, ".2.0")
showcase = showcase[showcase$FieldID %in% int_plot_smri$response, ]
label = data.frame(response = showcase$FieldID, name = showcase$Field)
int_plot_smri = merge(label, int_plot_smri ,by='response')
int_plot_smri = int_plot_smri[,!colnames(int_plot_smri) %in% "response"]
colnames(int_plot_smri)[colnames(int_plot_smri) %in% 'name'] = "response"

# label parcellations
showcase = data.frame(fread("/Users/jonathanhess/Documents/ukbiobank/phenotypes/Data_Dictionary_Showcase.tsv"))
showcase$FieldID = paste0("f.",showcase$FieldID, ".2.0")
showcase = showcase[showcase$FieldID %in% res$response, ]
label = data.frame(response = showcase$FieldID, name = showcase$Field)


res = merge(label, res, by='response')

smri_res = res

smri_res  = data.frame(fieldid = smri_res$response, 
                      response = smri_res$name, 
                      pred = smri_res$pred,
                      estimate = smri_res$Estimate,
                      se = smri_res$Std..Error,
                      pval = smri_res$Pr...t..,
                      n_sample = smri_res$n_total)


# combine all the statistics into one dataframe
harm_res$category[grepl("contemplated", harm_res$response)] = "ever_contemplated_self_harm"
harm_res$category[grepl("selfharm", harm_res$response)]  ="ever_engaged_self_harm"

smri_res$category = 'smri'
mental_health_stats$category = 'mental_health'

all_res = ldply(list( smri_res, harm_res, mental_health_stats))
all_res = all_res[grepl("risk|resilience", all_res$pred), ]
all_res$pred = gsub("risk:resilience","resilience:risk", all_res$pred)
table(all_res$pred)

all_res = ddply(all_res, .(pred), transform, fdr = p.adjust(pval, 'fdr'))

sum = ddply(all_res, .(pred, response), summarize, sum(fdr < .05))
sum[sum$`sum(fdr < 0.05)` == 1, ]

# export statistics
saveRDS(all_res, file="~/Google Drive/mac_storage/Manuscripts/Resilience/scz/scz-bull_special-issue/draft/manuscript_docs/circulate/edits/final_to_submit/revision/ukbb_analysis_stats_combined_part-2.Rdata")

# combine all plots into one data frame
all_plot = ldply(list(mental_health = int_plot_mental_health,
                      ever_contemplated_self_harm = int_plot_ideation,
                      ever_engaged_self_harm = int_plot_harm,
                      smri = int_plot_smri), .id='category')

saveRDS(all_plot, file="~/Google Drive/mac_storage/Manuscripts/Resilience/scz/scz-bull_special-issue/draft/manuscript_docs/circulate/edits/final_to_submit/revision/plot_intteractions_ukbb_analysis_combine_part-2.Rdata")


 

# mental health items
require(data.table);
require(ggpubr)
require(plyr);
require(dplyr)

# load ukbb compiled data 
ukbb = readRDS("~/Google Drive/mac_storage/Manuscripts/Resilience/scz/scz-bull_special-issue/draft/manuscript_docs/circulate/edits/final_to_submit/revision/ukbb_compiled_data_qc_filtered.Rdata")
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


# == run analysis of smri
ukb_smri = data.frame(fread("~/Documents/ukbiobank/ukb_smri_data.txt", nrows=1))
iters = colnames(ukb_smri)[-1]
iters = iters[grepl("\\.2.0", iters)]
iters = iters[iters %in% colnames(ukbb)]

grab = ukbb[,colnames(ukbb) %in% c(iters, "IID", "ICD10_F20.Schizophrenia", "GenomeWidePC_1", "GenomeWidePC_2","GenomeWidePC_3", "GenomeWidePC_4","GenomeWidePC_5","resilience","risk","f.25000.2.0",
                                   "f.25756.2.0", "f.25757.2.0", "f.25758.2.0", "f.25759.2.0", "f.21002.2.0", "f.21001.2.0", "f.12144.2.0", "f.31.0.0_Sex",
                                   "f.21003.2.0_Age.when.attended.assessment.centre")]
grab = grab[!rowSums(is.na(grab))/ncol(grab) > 0.5, ]
grab = grab[,!colSums(is.na(grab)/nrow(grab) > 0.5)]



# label parcellations
showcase = data.frame(fread("/Users/jonathanhess/Documents/ukbiobank/phenotypes/Data_Dictionary_Showcase.tsv"))
showcase$FieldID = paste0("f.",showcase$FieldID, ".2.0")
showcase = showcase[showcase$FieldID %in% colnames(grab), ]
label = data.frame(response = showcase$FieldID, name = showcase$Field)
# label = label[grepl("rostralmiddlefrontal|parsoperc", label$name),]
# label = label[grepl("thickness|Volume of", ignore.case = T, label$name), ]
label = label[label$name %in% c("Volume of rostralmiddlefrontal (right hemisphere)",
                        "Volume of rostralmiddlefrontal (left hemisphere)",
                        "Mean thickness of rostralmiddlefrontal (right hemisphere)",
                        "Mean thickness of rostralmiddlefrontal (left hemisphere)",
                        "Mean thickness of parsopercularis (left hemisphere)"),]

# enhanced well-being mediated in part by brain changes associated with resilience?
m1 = merge(grab, ukbb[,colnames(ukbb) %in% c("IID", iter)], by='IID')
label = label[label$response %in% colnames(m1), ]
 

# code "well-being" as absence of any of mental health items 
showcase = data.frame(fread("/Users/jonathanhess/Documents/ukbiobank/phenotypes/Data_Dictionary_Showcase.tsv"))
showcase = showcase[grepl("mental", ignore.case = T, showcase$Path), ]
showcase = showcase[grepl("Psychosocial factors", showcase$Path),]
showcase = showcase[showcase$ValueType %in% "Categorical single",]
showcase = showcase[grepl("ACE", showcase$Notes), ]
showcase = showcase[showcase$Participants >= 400e3, ]
grab = m1[,colnames(m1) %in% paste0("f.",showcase$FieldID,".0.0")]
rownames(grab) = m1$IID
grab[grab < 0] = NA
# grab[grab > 0] = -9
# grab[grab == 0] = 1
# grab[grab < 0] = 0
wellbeing_score = rowSums(grab , na.rm=T)

wellbeing_score = max(wellbeing_score) - (wellbeing_score) + min(wellbeing_score)
hist(wellbeing_score)

m1$wellbeing_score = wellbeing_score

hist(m1$wellbeing_score)

m1$wellbeing_score = RNOmni::RankNorm(m1$wellbeing_score)

stats_list = list();
for(x in 1:nrow(label)){
      # scale brain imaging measure
      item = m1[,label$response[[x]]]
      m1$smri = scale(item)
      
      frm = formula(paste("wellbeing_score ~ smri*resilience + f.21003.2.0_Age.when.attended.assessment.centre +  f.31.0.0_Sex + GenomeWidePC_1 +
                  GenomeWidePC_2 + GenomeWidePC_3 + GenomeWidePC_4 + GenomeWidePC_5 + f.25000.2.0 +
                                        f.25756.2.0 + f.25757.2.0 + f.25758.2.0 + f.25759.2.0 + f.21002.2.0 + f.21001.2.0 + f.12144.2.0"))
      
      fit = lm(frm, data = m1)
      stats = summary(fit)$coefficients
      stats = data.frame(pred = rownames(stats), stats)
      stats$response = "wellbeing"
      stats$model = x
      stats$smri = label$response[[x]]
      stats$n_total = nrow(m1)
      stats_list[[x]] = stats[grepl("smri|resilience|risk", stats$pred), ]
}
res = ldply(stats_list)
res = ddply(res, .(pred), transform, fdr = p.adjust(Pr...t.., 'fdr'))
sig = res[res$fdr < .05,]
sig

 

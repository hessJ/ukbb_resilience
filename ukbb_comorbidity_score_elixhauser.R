# compute co-morbidity score for ukbb participants
require(plyr)
require(data.table)
require(dplyr)

# load ukbb compiled data 
ukbb = readRDS("~/Google Drive/mac_storage/Manuscripts/Resilience/scz/scz-bull_special-issue/draft/manuscript_docs/circulate/edits/final_to_submit/revision/ukbb_compiled_data_qc_filtered.Rdata")
ukbb = data.frame(ukbb)

# grab columns with ICD-10 coded diagnoses 
ids = colnames(ukbb)[grepl("ICD", colnames(ukbb))]
grab = ukbb[,ids]
grab = data.frame(id = ukbb$IID, grab)
colnames(grab) = gsub("ICD10_", "", colnames(grab))
colnames(grab) = unlist(lapply(strsplit(colnames(grab), "[.]"), function(x) x[[1]]))

# munge data to compute comorbidity scores
require('comorbidity')
grab$id = factor(grab$id  )
melted = reshape2::melt(grab, .id='id')
melted = melted[melted$value == 1, ]
melted = melted[,colnames(melted) %in% c("id","variable")]

elixhauser <- comorbidity(x = melted, id = "id", code = "variable", map = "elixhauser_icd10_quan", assign0 = FALSE)
elixhauser$psycho = 0 # remove schizophrenia from score

coscores = score(elixhauser, assign0 = FALSE, weights = 'swiss')
coscores_unweighted = score(elixhauser, assign0 = FALSE, weights = NULL)

coscores_df = data.frame(IID = elixhauser$id, comorbidity_score = coscores, comorbidity_score_unweighted = coscores_unweighted)

missing_id = ukbb$IID[!ukbb$IID %in% coscores_df$IID]
missing_df = data.frame(IID = missing_id, comorbidity_score = 0)
all_scores = ldply(list(coscores_df, missing_df))

ukbb = merge(ukbb, all_scores, by='IID')

saveRDS(ukbb, file="~/Google Drive/mac_storage/Manuscripts/Resilience/scz/scz-bull_special-issue/draft/manuscript_docs/circulate/edits/final_to_submit/revision/ukbb_compiled_data_qc_filtered-w-comorbidity-score.Rdata")


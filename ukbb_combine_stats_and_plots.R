
# compile stats into a single table for supplement
l1 = readRDS("~/Google Drive/mac_storage/Manuscripts/Resilience/scz/scz-bull_special-issue/draft/manuscript_docs/circulate/edits/final_to_submit/revision/ukbb_analysis_stats_combined_part-1.Rdata")
l2 = readRDS("~/Google Drive/mac_storage/Manuscripts/Resilience/scz/scz-bull_special-issue/draft/manuscript_docs/circulate/edits/final_to_submit/revision/ukbb_analysis_stats_combined_part-2.Rdata")

all = ldply(list(l1, l2))
all$pred = gsub("risk:resilience","resilience:risk",all$pred)

# apply FDR correction to each X variable, per category
all = ddply(all, .(pred, category), transform, fdr = p.adjust(pval,'fdr'))
ddply(all, .(pred), summarize, sum(fdr < .05))

c1 = all[all$pred %in% "resilience:risk",]
c1 = c1[order(c1$category),]
fwrite(c1, file="~/Desktop/stats.txt", quote=F, row.names=F, sep="\t")

c1_sig = c1[c1$fdr < .05, ]
table(c1_sig$response)
table(c1_sig$category)

length(unique(c1_sig$response))

# load plot data -- plot only associatons that are significant
p1 = readRDS("~/Google Drive/mac_storage/Manuscripts/Resilience/scz/scz-bull_special-issue/draft/manuscript_docs/circulate/edits/final_to_submit/revision/plot_intteractions_ukbb_analysis_combine_part-1.Rdata")
p2 = readRDS("~/Google Drive/mac_storage/Manuscripts/Resilience/scz/scz-bull_special-issue/draft/manuscript_docs/circulate/edits/final_to_submit/revision/plot_intteractions_ukbb_analysis_combine_part-2.Rdata")
all_plots = ldply(list(p1,p2))
all_plots$response = ifelse(is.na(all_plots$response) == T, all_plots$id, all_plots$response)
unique(all_plots$response)
all_plots$response[all_plots$response %in% "icd10_scz"]  = "ICD10_F20.Schizophrenia"
all_plots$response[all_plots$response %in% "ever_contemplated_self_harm"] ="contemplated"
all_plots$response[all_plots$response %in% "ever_engaged_in_self_harm"] ="selfharm"
all_plots$category = as.character(all_plots$category)
all_plots$response = ifelse(all_plots$category %in% "icd10_broad", gsub("[,]", "/", all_plots$response), all_plots$response)


all_plots$response = gsub("fluid_intelligence",'intel', all_plots$response)
all_plots$response = gsub("education",'educ', all_plots$response)

sig_plots = all_plots[all_plots$response %in% c1_sig$response, ]
length(table(sig_plots$response))

# missing = c1_sig[!c1_sig$response %in% sig_plots$response, ]
# missing


require(ggplot2)
sig_plots$response = gsub("ICD10_F20.Schizophrenia","F20 Schizophrenia", sig_plots$response)

png("~/Desktop/scz_plot.png",res=300,units="in",height=2.5,width=3.5)
ggplot(sig_plots[sig_plots$response %in% "F20 Schizophrenia",], aes(x = risk, col = modx_group, lty = modx_group, y = outcome)) +
  geom_line() +
  scale_color_manual("Resilience", values=c("dodgerblue3","firebrick3","forestgreen")) +
  scale_linetype_manual(NULL, values=c(1,2,3)) +
  guides(linetype="none") +
  facet_wrap(~response,scale='free_y') +
  theme_bw() +
  xlab("Schizophrenia PRS") +
  ylab("Probability of diagnosis") +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=11, color='black'))
dev.off()


sig_plots$label = gsub("ICD10_", "", sig_plots$response,30)
sig_plots$label = gsub("[.]", " ", sig_plots$label)

png("~/Desktop/icd_narrow.png",res=300,units="in",height=4,width=8)
ggplot(sig_plots[sig_plots$category %in% "icd10_narrow",], aes(x = risk, col = modx_group, lty = modx_group, y = outcome)) +
  geom_line() +
  scale_color_manual("Resilience", values=c("dodgerblue3","firebrick3","forestgreen")) +
  scale_linetype_manual(NULL, values=c(1,2,3)) +
  facet_wrap(~stringr::str_wrap(label,30),scale='free_y') +
  theme_bw() +
  guides(linetype="none") +
  xlab("Schizophrenia PRS") +
  ylab("Probability of diagnosis") +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=11, color='black'))
dev.off()

png("~/Desktop/icd_broad.png",res=300,units="in",height=4,width=11)
ggplot(sig_plots[sig_plots$category %in% "icd10_broad",], aes(x = risk, col = modx_group, lty = modx_group, y = outcome)) +
  geom_line() +
  scale_color_manual("Resilience", values=c("dodgerblue3","firebrick3","forestgreen")) +
  scale_linetype_manual(NULL, values=c(1,2,3)) +
  facet_wrap(~stringr::str_wrap(response,25),scale='free_y', ncol = 5) +
  theme_bw() +
  guides(linetype="none") +
  xlab("Schizophrenia PRS") +
  ylab("Probability of diagnosis") +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=11, color='black'))
dev.off()


sig_plots$response = gsub("educ", "Educational attainment", sig_plots$response)
sig_plots$response = gsub("intel", "Fluid intelligence", sig_plots$response)

png("~/Desktop/cognitive.png",res=300,units="in",height=2.5,width=6.25)
ggplot(sig_plots[sig_plots$category %in% c("education","intelligence"),], aes(x = risk, col = modx_group, lty = modx_group, y = outcome)) +
  geom_line() +
  scale_color_manual("Resilience", values=c("dodgerblue3","firebrick3","forestgreen")) +
  scale_linetype_manual(NULL, values=c(1,2,3)) +
  facet_wrap(~response,scale='free_y') +
  theme_bw() +
  guides(linetype="none") +
  xlab("Schizophrenia PRS") +
  ylab("Fitted values") +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=11, color='black'))
dev.off()


png("~/Desktop/mental_health.png",res=300,units="in",height=5.25,width=7.5)
ggplot(sig_plots[sig_plots$category %in% "mental_health",], aes(x = risk, col = modx_group, lty = modx_group, y = outcome)) +
  geom_line() +
  scale_color_manual("Resilience", values=c("dodgerblue3","firebrick3","forestgreen")) +
  scale_linetype_manual(NULL, values=c(1,2,3)) +
  facet_wrap(~stringr::str_wrap(response,25),scale='free_y') +
  theme_bw() +
  guides(linetype="none") +
  xlab("Schizophrenia PRS") +
  ylab("Fitted values") +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=11, color='black'))
dev.off()


sig_plots$response = gsub("contemplated", "Ever contemplated self-harm", sig_plots$response)
sig_plots$response = gsub("selfharm", "Ever self-harmed", sig_plots$response)

png("~/Desktop/self_harm.png",res=300,units="in",height=4,width=3.5)
ggplot(sig_plots[grepl("harm", sig_plots$category),], aes(x = risk, col = modx_group, lty = modx_group, y = outcome)) +
  geom_line() +
  scale_color_manual("Resilience", values=c("dodgerblue3","firebrick3","forestgreen")) +
  scale_linetype_manual(NULL, values=c(1,2,3)) +
  facet_wrap(~stringr::str_wrap(response,30),scale='free_y', ncol=1) +
  theme_bw() +
  guides(linetype="none") +
  xlab("Schizophrenia PRS") +
  ylab("Probability of self-harm measure") +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=11, color='black'))
dev.off()


l1 = readRDS("~/Google Drive/mac_storage/Manuscripts/Resilience/scz/scz-bull_special-issue/draft/manuscript_docs/circulate/edits/final_to_submit/revision/ukbb_analysis_stats_combined_part-1.Rdata")
l2 = readRDS("~/Google Drive/mac_storage/Manuscripts/Resilience/scz/scz-bull_special-issue/draft/manuscript_docs/circulate/edits/final_to_submit/revision/ukbb_analysis_stats_combined_part-2.Rdata")

all = ldply(list(l1, l2))
all$pred = gsub("risk:resilience","resilience:risk",all$pred)

all = ddply(all, .(pred), transform, fdr = p.adjust(pval,'fdr'))
ddply(all, .(pred), summarize, sum(fdr < .05))

c1 = all[all$pred %in% "resilience",]
c1 = c1[order(c1$category),]

c1_sig = c1[c1$fdr < .05 & c1$category %in% "smri",]
c1_sig

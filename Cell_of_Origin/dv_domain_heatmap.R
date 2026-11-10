## =============================================================================
## dv_domain_heatmap.R
##
## Heat map of the probability with which malignant cells of the pons-like (PL) and
## thalamus-like (TL) tumours map onto the 13 dorsoventral (DV) progenitor domains of the
## developing human spinal cord (reference: GEO GSE171890, Rayon et al., 2021).
##
## Run from the Cell_of_Origin/ folder:   Rscript dv_domain_heatmap.R
##
## Input   data/dv_domain_cellstate_counts.tsv   written by dv_domain_label_transfer.R
##         (sample, group, predicted_domain, mp_state, n)
## Output  dv_domain_mapping.pdf
##
## Probability: within a group (PL or TL) and a malignant state, the cells are counted per
## domain and divided by the number of cells of that state; the six states are then averaged,
## so every state weighs the same whatever its size and each group's probabilities sum to 1.
## =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

f_counts <- 'data/dv_domain_cellstate_counts.tsv'
f_pdf    <- 'dv_domain_mapping.pdf'

## ---- constants --------------------------------------------------------------
dv_levels    <- c('RP','dp1','dp2','dp3','dp4','dp5','dp6','p0','p1','p2','pMN','p3','FP')  # dorsal -> ventral
tumor_states <- c('AC-like','OPC-like','MES-like','OC-like','S','G2M')
group_labels <- c(PL = 'Pons-like', TL = 'Thalamus-like')
dv_colors    <- setNames(c("#8C1A11","#C6402D","#DE6E51","#ED9364","#ED9364","#F3BE8C","#FCF1DC",
                           "#F3C7C2","#EDA3B5","#E671A0","#CC4494","#9F207B","#6F1473"), dv_levels)

# mp_state is written with dots (AC.like), the state names above with hyphens (AC-like)
canon <- function(x) tolower(gsub('[-_. ]', '', as.character(x)))

## ---- summary: counts -> probability of each domain per group -----------------
counts <- read.delim(f_counts, stringsAsFactors = FALSE)
resdf <- data.frame(Var1 = counts$predicted_domain, Var2 = counts$mp_state, Freq = counts$n,
                    tumorGroup = counts$group, stringsAsFactors = FALSE)

stopifnot(all(canon(resdf$Var2) %in% canon(tumor_states)))     # every state in the table is a malignant state
bad <- setdiff(unique(resdf$Var1), dv_levels)
if (length(bad)) stop('predicted domains outside dv_levels: ', paste(bad, collapse = ', '))
resdf$Var1 <- factor(resdf$Var1, levels = dv_levels)

message('malignant cells per group: ',
        paste(names(tapply(resdf$Freq, resdf$tumorGroup, sum)), tapply(resdf$Freq, resdf$tumorGroup, sum), collapse = ', '))

tmp <- resdf %>%
  group_by(tumorGroup, Var2) %>%
  filter(sum(Freq) > 0) %>%                       # a state absent from a group is dropped (avoids 0/0)
  mutate(Perc = Freq / sum(Freq)) %>%             # within (group, state): every state contributes 1
  group_by(tumorGroup) %>%
  mutate(n_state = n_distinct(Var2)) %>%          # denominator = number of states present in the group
  group_by(tumorGroup, Var1) %>%
  summarise(Probability = sum(Perc) / first(n_state), .groups = 'drop') %>%
  tidyr::complete(tumorGroup, Var1, fill = list(Probability = 0))   # domains without cells -> 0

# dorsal at the top, ventral at the bottom
tmp$Var1 <- factor(as.character(tmp$Var1), levels = rev(dv_levels))

# sanity check: the probabilities of each group sum to 1
chk <- tmp %>% group_by(tumorGroup) %>% summarise(check_sum = sum(Probability), n_domain = n())
print(chk)
stopifnot(all(abs(chk$check_sum - 1) < 1e-9), all(chk$n_domain == length(dv_levels)))

## ---- plot -------------------------------------------------------------------
p_dv <- tmp %>%
  ggplot(aes(y = Var1, fill = Probability, x = tumorGroup)) +
  geom_tile(color = 'white', linewidth = 0.5) +
  scale_fill_viridis_c() +
  theme_classic() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  # labels by name, independent of the order of PL / TL on the axis
  scale_x_discrete(labels = function(x) unname(ifelse(x %in% names(group_labels), group_labels[x], x))) +
  # label colours by the actual order of the domains
  theme(axis.text.y = element_text(color = dv_colors[levels(tmp$Var1)], face = 'bold')) +
  labs(x = '', fill = 'Probability', y = '')

ggsave(f_pdf, p_dv, width = 3.35, height = 4)
message('wrote ', f_pdf)

print(as.data.frame(tmp), digits = 4)

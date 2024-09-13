
library(pheatmap)
library(viridis)
library(reshape) 
library(ggpubr)

library(ggplot2)
library(viridis)
library(plyr)
library( gtools) 
library(dplyr)
out_path = r"(D:\OneDrive - HKUST Connect\Work\SingleCell\Multiplex-IHC\spatial_analysis\0623\csv)"
setwd(out_path)

meta = read.csv(r"(D:\OneDrive - HKUST Connect\Work\SingleCell\Multiplex-IHC\spatial_analysis\objects\meta.csv)")
name_lst = meta$ID
rownames(meta) = meta$ID
st_color_list = c('GPM'='#b62526','MTC'='#e97550','PPR'='#fb9e4f','B'='#439030','T_NK'='#2f609f',
                  'Myeloid'='#57b894','Neutrophil'='#97d68f','Olig'='#ffffa4',
                  'Vascular'='#b0aad1','Neural_cell'='#00a6ff','Tumor_unknown'='#fb9a99','Other'='#b9b9b9')

celltype_order =  c('GPM','MTC','PPR','Tumor_unknown','Neural_cell','T_NK','Myeloid','Neutrophil','Vascular','B','Other')
all_df = read.csv('all_df.csv')
all_df_filtered = read.csv('all_df.csv')
all_df_bak = all_df
all_df_filtered_bak = all_df_filtered
table(all_df$kmeans9)
table(all_df$celltype_new)
table(all_df_mapped$celltype_agg)
stat_df = data.frame(table(all_df_mapped$celltype_agg, all_df_mapped$group))
colnames(stat_df) = c('Celltype','Group','Freq')

stat_df$Celltype = factor(stat_df$Celltype, levels=celltype_order)
ggplot(stat_df, aes(fill=Celltype, y=Freq, x=Group)) + 
  geom_bar(position="fill", stat="identity") + theme_classic() + scale_fill_manual(values=st_color_list) 
ggsave('celltype_group_stat_0623.pdf', width=5, height=5, dpi=300)

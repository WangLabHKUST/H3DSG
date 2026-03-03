library(tidyverse)
library(viridis)     # 更好看的颜色（可选）
# library(MetBrewer) # 如果想要更学术的配色也可以用

# -----------------------
# 读入数据（以 Group1 为例）
# -----------------------
base_path = 'all_samples_full_27_5/'
df1 = read_csv(file.path(base_path, "grid_data_entropy_tumor_vs_non_malig_entropy_Group1.csv"))
df2 = read_csv(file.path(base_path, "grid_data_entropy_tumor_vs_non_malig_entropy_Group2.csv"))
density_lst = c(df1$density, df2$density)
# max_density = max(df1$density, df2$density)
# max_density = quantile(density_lst, 0.9999)
max_density = 20
head(sort(density_lst, decreasing = TRUE))
# group = 'Group1'
for(cutoff in c('auto','max')){
for(group in c('Group1', 'Group2')){
contour_num_mapping = c('Group1'=5, 'Group2'=10)
contour_num = contour_num_mapping[group]
df <- read_csv(file.path(base_path, sprintf("grid_data_entropy_tumor_vs_non_malig_entropy_%s.csv", group)))
head(df)
df = df[df$non_malig_entropy < 1.2,]
# 确认列名（根据你导出的字段调整）
# 假设列名为：entropy_tumor, non_malig_entropy, density

# 如果你的变量名不同，请改这里
var_x <- "entropy_tumor"
var_y <- "non_malig_entropy"
var_z <- "density"

# -----------------------
# 画图 - 两种主要风格，任选其一
# -----------------------

# 风格1：接近你原图的 pcolormesh + contour（推荐）
if(cutoff == 'max'){
  df[df$density > max_density, 'density'] = max_density
}
p <- ggplot(df, aes(x = .data[[var_x]], y = .data[[var_y]], z = .data[[var_z]])) +
  
  # 热力图主体（类似 pcolormesh）
  geom_raster(aes(fill = .data[[var_z]]), interpolate = TRUE) +   # interpolate ≈ gouraud
  
  # 等高线（contour）
  geom_contour(aes(z = .data[[var_z]]), 
               color = "black", 
               linewidth = 0.05, 
               bins = 15
              #  breaks = seq(0, max_density, length.out = 15)
               ) +          # 你原图第二个图用了20级，第一个用了10级
  
  
  # 坐标轴范围（根据你的数据调整）
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1.2), expand = FALSE) +
  
  # 主题调整
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right",
    legend.key.height = unit(1.2, "cm")
  ) +
  
  labs(
    title = "Group1",
    x = "entropy_tumor",
    y = "non_malig_entropy"
  )
  ### 
  this_max = max(df$density)
  if(cutoff == 'auto'){
    p = p + scale_fill_distiller(
      palette   = "YlOrRd",       # 核心：就是这里写 "RdYlBu"
      direction = 1,              # 1 = 低值蓝 → 高值红（默认）
      name      = "Density",       # 图例标题
      breaks    = seq(0, this_max, by=round(this_max/5))
    )
  }
  if(cutoff == 'max'){
    p = p + scale_fill_distiller(
      palette   = "YlOrRd",       # 核心：就是这里写 "RdYlBu"
      direction = 1,              # 1 = 低值蓝 → 高值红（默认）
                                  # -1 = 反转 → 低值红 → 高值蓝
      limits    = c(0, max_density),        # 匹配你原来的 vmin/vmax
      # breaks    = 0:max_density,
      name      = "Density"       # 图例标题
    )
  }


# 保存（可选）
ggsave(file.path(base_path, sprintf("density_r.%s.%s.png", group, cutoff)), p, width = 6, height = 5.2, dpi = 300)
}
}
library(ggplot2)


#######################
# set working directory
setwd("/home/tcigene/jade/ncounter/data/1216/")

folder = "/home/tcigene/jade/ncounter/result/plot/2746_51/"
plant_name = "余甘子化合物"
team_name = c("控制組", "UVA","余甘子化合物-6 6h", "余甘子化合物-6 24h")
team = c("mock_nor", "UVA1", "X2751.6h1", "X2751.24h1")
color = c("cornflowerblue", "#148f77", "#f5b7b1", "#e74c3c")
team_num = length(team)

############################
# get the ncounter table csv
ibd = read.csv("20191209_余甘子化合物-1光老化基因檢測.csv")
ibd = cbind(ibd[, 1], ibd[, 8:dim(ibd)[2]])
ibd = ibd[3:173,]           
colnames(ibd)[1] = "gene"

sample_name = colnames(ibd)[seq(2, length(ibd), 3)]
##############
# Select data
# 1. Select 2 of 3 replicate
# 2. Assign 0 if 2 of the replicate lower than 30

# Select 2 of 3 replicate: mock 

mock = matrix(NA, nrow = dim(ibd)[1], ncol = 2)

for (row in 1:dim(ibd)[1]){
  
  a = ibd[row, 2:4]
  
  r1 = a[1]
  r2 = a[2]
  r3 = a[3]
  
  d12 = abs(r1 - r2)
  d23 = abs(r2 - r3)
  d13 = abs(r1 - r3)
  
  # Select the 2 close value from 3 replicates
  
  if (min(d12, d23, d13) == d12){
    mock[row, 1] = as.numeric(r1)
    mock[row, 2] = as.numeric(r2)
    
  } else if (min(d12, d23, d13) == d23){
 
    mock[row, 1] = as.numeric(r2)
    mock[row, 2] = as.numeric(r3)
    
  } else if (min(d12, d23, d13) == d13){
    
    mock[row, 1] = as.numeric(r1)
    mock[row, 2] = as.numeric(r3)
    
  }
  
}


# Select 2 of 3 replicate: other samples

data = matrix(NA, nrow = dim(ibd)[1], ncol = (length(sample_name) - 1) * 2)

for (col in seq(5, dim(ibd)[2], 3)){
  
  for (row in 1:dim(ibd)[1]){
    
    data_col = ((col-2)/3) * 2
      
    a = ibd[row, col:(col + 2)]
    
    r1 = a[1]
    r2 = a[2]
    r3 = a[3]
    
    d12 = abs(r1 - r2)
    d23 = abs(r2 - r3)
    d13 = abs(r1 - r3)
    
    # Select the 2 close value from 3 replicates
    
    if (sum((r1 < 30), (r2 < 30), (r3 < 30)) >= 2){

      data[row, (data_col - 1)] = as.numeric(0)
      data[row, data_col] = as.numeric(0)
      
    } else if (min(d12, d23, d13) == d12){
      
      data[row, (data_col - 1)] = as.numeric(r1)
      data[row, data_col] = as.numeric(r2)
      
    } else if (min(d12, d23, d13) == d23){
      
      data[row, (data_col - 1)] = as.numeric(r2)
      data[row, data_col] = as.numeric(r3)
      
    } else if (min(d12, d23, d13) == d13){
      
      data[row, (data_col - 1)] = as.numeric(r1)
      data[row, data_col] = as.numeric(r3)
      
    }
    
  }
  
}


#################
# Normalized data

# calculate the mean of mock

mean_mock = matrix(NA, nrow  = 171, ncol = 1)

colnames(mean_mock) = "mean_mock"
rownames(mean_mock) = ibd$gene

for (row in 1:171){
  
  mean_mock[row ,1] = mean(mock[row,])
  
}

# every data / mean_mock

nor_data = matrix(NA, nrow = 171, ncol = dim(data)[2])

rownames(nor_data) = ibd$gene

for (col in 1: dim(data)[2]){
  
  for (row in 1:171){
    
    nor_data[row, col] = data[row,col] / mean_mock[row, 1]
    
  }
}

######################
# Mean of each samples

mean = matrix(NA, nrow = 171, ncol = (dim(nor_data)[2] / 2))

colnames(mean) = sample_name[2:length(sample_name)]
rownames(mean) = ibd$gene

for (col in 1:dim(mean)[2]){
  
  for (row in 1:171){
    
    mean[row, col] = mean(nor_data[row, ((col*2 - 1):(col*2))])
      
  }
}

mock_nor = rep(1, 181)

mean = cbind(mock_nor, mean)

######################
# Std of each samples

std = matrix(NA, nrow = 171, ncol = (dim(nor_data)[2] / 2))

colnames(std) = sample_name[2:length(sample_name)]
rownames(std) = ibd$gene

for (col in 1:dim(std)[2]){
  
  for (row in 1:171){
    
    std[row, col] = sd(nor_data[row, ((col*2 - 1):(col*2))])
    
  }
}

#########################
# p-value of each samples 

# delete the UNA columns
p = matrix(NA, nrow = 171, ncol = (dim(nor_data)[2] / 2) - 1)

rownames(p) = ibd$gene

# delete UVA colnames
colnames(p) = sample_name[3:length(sample_name)]
  
for (col in 1:dim(p)[2]){
  
  for (row in 1:171){
    
    # compare with UVA team
    p[row, col] = t.test(c(nor_data[row,1], nor_data[row,2] + 0.00001), c(nor_data[row, (col*2 + 1)], nor_data[row, (col*2 + 2)] + 0.00001), alternative = "two.sided", paired = F)$p.value
    
  }
}


#######################
# Sort by genetic group

# 抗氧化
g1 = c('SOD1','SOD2','GPX1','CAT')
g1_sign = c('up', 'up', 'up', 'up')
g1_cat = rep("Anti-oxidation", 4)

# 抗老
g2_1 = c('CCT2','CCT5','CCT6A','CCT7','CCT8')
g2_2 = c('Pink1', 'Atg1','Atg8','SIRT1','FOXO')
g2_3 = c('PARP1','PARP2','NADSYN','MRPS5','Ubl-5','SOD3')
g2 = c('CCT2','CCT5','CCT6A','CCT7','CCT8','Pink1', 'Atg1','Atg8','SIRT1','FOXO','PARP1','PARP2','NADSYN','MRPS5','Ubl-5','SOD3') #,'Parkin'
g2_sign = c('up', 'up', 'up', 'up', 'up', 'up', 'up', 'up', 'up', 'up', 'down', 'down', 'up', 'up', 'up', 'up')
g2_cat = rep("Anti-aging", 16)

# DNA修復
g3_1 = c('UNG','OGG1','MPG','APEX1','ERCC1','ERCC6')
g3_2 = c('XPA','XRCC1','XRCC5','MSH2','MLH1','MSH6')
g3 = c('UNG','OGG1','MPG','APEX1','ERCC1','ERCC6','XPA','XRCC1','XRCC5','MSH2','MLH1','MSH6')
g3_sign = c('up', 'up', 'up', 'up', 'up', 'up', 'up', 'up', 'up', 'up', 'up', 'up')
g3_cat = rep("DNA repair", 12)

# 免疫
g4 = c('IL-1B','IL-8','IL-6','IL-10','IL-18','TNF-a')
g4_sign = c('up', 'up', 'up', 'up', 'up', 'up')
g4_cat = rep("Immunity", 6)

# 美白-抗黑色素生成
g5 = c('TYR','TYRP1','MC1R','MITF')
g5_sign = c('down', 'down', 'down', 'down')
g5_cat = rep("Whitening", 4)

# 膠原蛋白 彈力蛋白 玻尿酸 合成 降解
g6_1 = c('COL1A1','COL1A2','COL4A1','COL4A4','COL4A5')
g6_2 = c('MMP1','MMP9','MMP2','TIMP1')
g6_3 = c('ELN','FBN1','LOX','HAS2','HAS3')   
g6 = c('COL1A1','COL1A2','COL4A1','COL4A4','COL4A5','MMP1','MMP9','MMP2','TIMP1','ELN','FBN1','LOX','HAS2','HAS3')
g6_sign = c('up', 'up', 'up', 'up', 'up', 'down', 'down', 'down', 'up', 'up', 'up', 'up', 'up', 'up')
g6_cat = rep("Collagen", 14)

# 抗發炎
g7_1 = c('IL-1B','IL-8','IL-6','IL-10','IL-18')
g7_2 = c('TNF-a','IL-16','IL23','IL12A')
g7_3 = c('IFNG','TGFB','IL3','IL4')
g7 = c('IL-1B','IL-8','IL-6','IL-10','IL-18','TNF-a','IL-16','IL23','IL12A','IFNG','TGFB','IL3','IL4')
g7_sign = c('down', 'down', 'down', 'up', 'down', 'down', 'down', 'down', 'down', 'down', 'up', 'down', 'up')
g7_cat = rep("Anti-inflammatory", 13)

# 心血管保健
g8_1 = c('PTGIS','NOS3','EDN1','PLAT','PROC','VWF')
g8_2 = c('F3','SERPINE1','PDGFC','FGF2','IGF2BP3','IGF1R')
g8_3 = c('IL-8','IL-6','ICAM1','VCAM1','CASP8')
g8 = c('PTGIS','NOS3','EDN1','PLAT','PROC','VWF','F3','SERPINE1','PDGFC','FGF2','IGF2BP3','IGF1R','IL-8','IL-6','ICAM1','VCAM1','CASP8')
g8_sign = c('up', 'up', 'down', 'up', 'up', 'down', 'down', 'down', 'down', 'down', 'down', 'down', 'down', 'down', 'down', 'down', 'down')
g8_cat = rep("Cardiovascular care", 17)

# 皮膚角質保濕
g9_1 = c('Tgm1','Krt1','Keratin 10','Keratin 14','AQP3')
g9_2 = c('FLG-F','SMPD1','GBA','HAS2','HAS3')
g9 = c('Tgm1','Krt1','Keratin 10','Keratin 14','AQP3','FLG-F','SMPD1','GBA','HAS2','HAS3')
g9_sign = c('up', 'up', 'up', 'up', 'up', 'up', 'up', 'up', 'up', 'up')
g9_cat = rep("Skin moisturizing", 10)

# 脂肪肝
g10 = c('SREBP-1c (SREBF1)','PPAR-g','PPAR-a','SCD1 (SCD)','ACC (ACACA)')
g10_sign = c('down', 'down', 'up', 'down', 'down')
g10_cat = rep("Fatty liver", 5)

# 提升HDL
g11 = c('CETP','SCARB1','apoA-I (APA1)','LDLR','ABCA1')
g11_sign = c('up', 'up', 'up', 'up', 'up')
g11_cat = rep("Increase HDL", 5)

# 健髮
g12_1 = c('SRD5A1','SRD5A2','AR','KROX20')
g12_2 = c('SCF','VEGFA','IGF1','TGFB','BDNF')
g12 = c('SRD5A1','SRD5A2','AR','KROX20','SCF','VEGFA','IGF1','TGFB','BDNF')
g12_sign = c('down', 'down', 'down', 'up', 'up', 'up', 'up', 'down', 'down')
g12_cat = rep("Hair", 9)

# 端粒酶活性
g13 = c('TERT','TERC','RTEL1')
g13_sign = c('up', 'up', 'up')
g13_cat = rep("Telomerase activity", 3)

# 呼吸道過敏
g14_1 = c('CD40','ERBB2','LIF','MALT1','NCK1','PAF1')
g14_2 = c('DYNLL2','GRK5','PSMD4','RDH10','RELB','SCARF1')
g14_3 = c('TNFSF14','ABR','IL13','IL4R','IL5RA','RELA')
g14 = c('CD40','ERBB2','LIF','MALT1','NCK1','PAF1','DYNLL2','GRK5','PSMD4','RDH10','RELB','SCARF1','TNFSF14','ABR','IL13','IL4R','IL5RA','RELA')
g14_sign = c('down', 'down', 'down', 'down', 'down', 'down', 'down', 'down', 'down', 'down', 'down', 'down', 'down', 'down', 'down', 'down', 'down', 'down')
g14_cat = rep("Allergy", 18)

# combine g_1 to g_14  
g = c(g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11, g12, g13, g14)
g_sign = c(g1_sign, g2_sign, g3_sign, g4_sign, g5_sign, g6_sign, g7_sign, g8_sign, g9_sign, g10_sign, g11_sign, g12_sign, g13_sign, g14_sign)
g_cat = c(g1_cat, g2_cat, g3_cat, g4_cat, g5_cat, g6_cat, g7_cat, g8_cat, g9_cat, g10_cat, g11_cat, g12_cat, g13_cat, g14_cat)

g_result = cbind(g_cat, g, g_sign)




index = matrix(NA, nrow = length(g), ncol = 1)

for (i in 1:length(g)){
  
  index[i] = which(rownames(p) == g[i])
  
}

mean = mean[index,]
std = std[index, ]
p = p[index, ]


p_star = matrix(NA, nrow = dim(p)[1], ncol = dim(p)[2])

rownames(p_star) = rownames(p)
colnames(p_star) = colnames(p)

for (col in 1:dim(p)[2]){
  for(row in 1:dim(p)[1]){
    
    if (mean[row, (col + 1)] == 0){
      p_star[row, col] = ""
    } else if (p[row, col] < 0.001){
      p_star[row, col] = "***"
    } else if(p[row, col] < 0.01){
      p_star[row, col] = "**"
    } else if(p[row, col] < 0.05){
      p_star[row, col] = "*"
    } else{
      p_star[row, col] = ""
    }
    
  }
}

result = cbind(g_result, mean,  rep("", length(g)), std, rep("", length(g)), p)

path = paste(folder, plant_name, "_result.csv", sep = "")
write.csv(result, path)

# combine the category, gene, mean, std, p-value, p_significant
mean_team = mean[,team]
std_team = std[,team[2:length(team)]]
# delete the UNA columns
p_team = p[,team[3:length(team)]]
p_star_team = p_star[,team[3:length(team)]]

result = cbind(g_result, mean_team,  rep("", length(g)), std_team, rep("", length(g)), p_team, rep("", length(g)), p_star_team)

path = paste(folder, plant_name, "_result.csv", sep = "")
write.csv(result, path)

# select 1) mean without 0 ; 2) at least one with up or down regulate(also with significant) 

# delete mean with 0
a1 = which(mean_team[,2] == 0)
a2 = which(mean_team[,3] == 0)
a3 = which(mean_team[,4] == 0)

u1 = union(a1, a2)
u3 = union(u1, a3)

r1 = result[-u3, ]

mean_delete0 = mean_team[-u3,]
std_delete0 = std_team[-u3,]
p_delete0 = p_team[-u3,]
p_star_delete0 = p_star_team[-u3,]

# select at least one with up or down regulate(also with significant)
# 1. get mean col1 and p col1  intersect from col 1, 2, 3, 4;  get the intersection row number

select = matrix(NA, nrow = dim(r1)[1], ncol = 1)

for (row in 1:dim(r1)[1]){
  
  if (r1[row,3] == "up"){
    
    if(length(intersect(which(p_delete0[row,] <= 0.05), which(mean_delete0[row,2:dim(mean_delete0)[2]] > 1))) > 0){
      
      select[row, 1] = row
      
    }
    
  } else if(r1[row,3] == "down"){
    
    if(length(intersect(which(p_delete0[row,] <= 0.05), which(mean_delete0[row,2:dim(mean_delete0)[2]] < 1))) > 0){
      
      select[row, 1] = row
      
    }
    
  }
  
}

select = na.omit(select)
result_sign = r1[select, ]

mean_sign = mean_delete0[select,]
std_sign = std_delete0[select,]
p_sign = p_delete0[select,]
p_star_sign = p_star_delete0[select,]


path = paste(folder, plant_name, "_result_significant.csv", sep = "")
write.csv(result_sign, path)




###############################################
###############################################
# group bar plot for SELECTED gene and category

############
# G1: 抗氧化

# set variables (1) gene_group (2) team_number (3) team_name
group_name = "抗氧化"
#group = g1
group = result_sign[which(result_sign[,1] == "Anti-oxidation"), 2]
gene_num = length(group)


# get team mean
value = mean_sign[group,]

# check if only one gene
if (gene_num == 1){
  
  # one gene
  avg = matrix(NA, nrow = 1*dim(mean_sign)[2], ncol = 1)
  
} else {
  
  # > 1 genes
  avg = matrix(NA, nrow = dim(value)[1] * dim(value)[2], ncol = 1)
}



# add the value the avg
if (gene_num == 1){
  
  # one gene
  for(col in 1:dim(mean_sign)[2]){
    avg[col, 1] = value[col]
    
  }
  
} else{
  
  # > 1 genes
  for (row in 1:dim(value)[1]){
    for(col in 1:dim(value)[2]){
      avg[i, 1] = value[row, col]
      i = i + 1
    }
  }
}


# get team std
value_std = std_sign[group,]

if (gene_num == 1){
  
  # 1 gene
  mock_std = rep(0, 1)
  
} else{
  
  # > 1 genes
  mock_std = rep(0, rep(dim(value_std)[1]))
}



if (gene_num == 1){
  
  # one gene
  value_std = c(mock_std, value_std)
} else{
  
  # > 1 gene
  value_std = cbind(mock_std, value_std)
}

value_std = as.data.frame(value_std)

sd = matrix(NA, nrow = dim(value_std)[1] * dim(value_std)[2], ncol = 1)

i = 1

for (row in 1:dim(value_std)[1]){
  for(col in 1:dim(value_std)[2]){
    sd[i, 1] = value_std[row, col]
    i = i + 1
  }
}


# get team p-value
value_p = p_star_sign[group,]

if (gene_num == 1){
  
  # 1 gene
  mock_p = rep("", 1)
} else {
  
  # > 1 genes
  mock_p = rep("", rep(dim(value_p)[1]))
}

# UVA p_value
if (gene_num == 1){
  
  # 1 gene
  uva_p = rep("", 1)
} else {
  
  # > 1 genes
  uva_p = rep("", rep(dim(value_p)[1]))
}


# combine the mock_p, uva_p, value_p
if (gene_num == 1){
  
  # 1 gene
  value_p = c(mock_p, uva_p,value_p)
  
} else {
  
  # > 1 genes
  value_p = cbind(mock_p, uva_p,value_p)
  
}




value_p = as.data.frame(value_p)

p_team = matrix(NA, nrow = dim(value_p)[1] * dim(value_p)[2], ncol = 1)

i = 1

for (row in 1:dim(value_p)[1]){
  for(col in 1:dim(value_p)[2]){
    p_team[i, 1] = as.character(value_p[row, col])
    i = i + 1
  }
}

# get the team name
team = rep(team_name, gene_num)

# get the gene name
gene = matrix(NA, nrow = gene_num*team_num, ncol = 1)

for (i in 1:length(group)){
  
  gene[(i*team_num-(team_num-1)):(i*team_num),1] = rep(group[i], team_num)
  
}


# add the gene, team, avg for barplot
data = data.frame(gene, team, avg, sd, p_team)

data = within(data, team <- factor(team, levels = team_name))


path = paste(folder, plant_name, "-", group_name, ".png", sep = "")


ggplot(data, aes(fill = team, y = avg, x = gene)) + 
  geom_bar(position = "dodge", stat = "identity") +
  scale_fill_manual(values = color) +
  geom_errorbar(aes(ymin = avg, ymax = avg + sd), width = .2, position = position_dodge(.9)) +
  theme(text = element_text(size = 20), plot.title = element_text(hjust = 0.5)) +
  xlab("Gene") + 
  ylab("Relative expression ratio") + 
  ggtitle(paste(plant_name, "-", group_name, spe = "")) +
  scale_y_continuous(breaks = seq(0, max(data$avg + data$sd + 0.2), 0.2)) +
  geom_text(aes(y = (avg + sd), label = p_team), position = position_dodge(width = 0.9), vjust = -0.1, size = 5)

ggsave(path, width = 7, height = 7, dpi = 300)





##################
# G3: DNA修復 g3

# set variables (1) gene_group (2) team_number (3) team_name
group_name = "DNA修復"
group = result_sign[which(result_sign[,1] == "DNA repair"), 2]
gene_num = length(group)


# get team mean
value = mean_sign[group,]

# check if only one gene
if (gene_num == 1){
  
  # one gene
  avg = matrix(NA, nrow = 1*dim(mean_sign)[2], ncol = 1)
  
} else {
  
  # > 1 genes
  avg = matrix(NA, nrow = dim(value)[1] * dim(value)[2], ncol = 1)
}



# add the value the avg
if (gene_num == 1){
  
  # one gene
  for(col in 1:dim(mean_sign)[2]){
    avg[col, 1] = value[col]
    
  }
  
} else{
  
  # > 1 genes
  for (row in 1:dim(value)[1]){
    for(col in 1:dim(value)[2]){
      avg[i, 1] = value[row, col]
      i = i + 1
    }
  }
}


# get team std
value_std = std_sign[group,]

if (gene_num == 1){
  
  # 1 gene
  mock_std = rep(0, 1)
  
} else{
  
  # > 1 genes
  mock_std = rep(0, rep(dim(value_std)[1]))
}



if (gene_num == 1){
  
  # one gene
  value_std = c(mock_std, value_std)
} else{
  
  # > 1 gene
  value_std = cbind(mock_std, value_std)
}

value_std = as.data.frame(value_std)

sd = matrix(NA, nrow = dim(value_std)[1] * dim(value_std)[2], ncol = 1)

i = 1

for (row in 1:dim(value_std)[1]){
  for(col in 1:dim(value_std)[2]){
    sd[i, 1] = value_std[row, col]
    i = i + 1
  }
}


# get team p-value
value_p = p_star_sign[group,]

if (gene_num == 1){
  
  # 1 gene
  mock_p = rep("", 1)
} else {
  
  # > 1 genes
  mock_p = rep("", rep(dim(value_p)[1]))
}


# UVA p_value
if (gene_num == 1){
  
  # 1 gene
  uva_p = rep("", 1)
} else {
  
  # > 1 genes
  uva_p = rep("", rep(dim(value_p)[1]))
}


# combine the mock_p, uva_p, value_p
if (gene_num == 1){
  
  # 1 gene
  value_p = c(mock_p, uva_p,value_p)
  
} else {
  
  # > 1 genes
  value_p = cbind(mock_p, uva_p,value_p)
  
}


value_p = as.data.frame(value_p)

p_team = matrix(NA, nrow = dim(value_p)[1] * dim(value_p)[2], ncol = 1)

i = 1

for (row in 1:dim(value_p)[1]){
  for(col in 1:dim(value_p)[2]){
    p_team[i, 1] = as.character(value_p[row, col])
    i = i + 1
  }
}

# get the team name
team = rep(team_name, gene_num)

# get the gene name
gene = matrix(NA, nrow = gene_num*team_num, ncol = 1)

for (i in 1:length(group)){
  
  gene[(i*team_num-(team_num-1)):(i*team_num),1] = rep(group[i], team_num)
  
}


# add the gene, team, avg for barplot
data = data.frame(gene, team, avg, sd, p_team)

data = within(data, team <- factor(team, levels = team_name))


path = paste(folder, plant_name, "-", group_name, ".png", sep = "")


ggplot(data, aes(fill = team, y = avg, x = gene)) + 
  geom_bar(position = "dodge", stat = "identity") +
  scale_fill_manual(values = color) +
  geom_errorbar(aes(ymin = avg, ymax = avg + sd), width = .2, position = position_dodge(.9)) +
  theme(text = element_text(size = 20), plot.title = element_text(hjust = 0.5)) +
  xlab("Gene") + 
  ylab("Relative expression ratio") + 
  ggtitle(paste(plant_name, "-", group_name, spe = "")) +
  scale_y_continuous(breaks = seq(0, max(data$avg + data$sd + 0.2), 0.2)) +
  geom_text(aes(y = (avg + sd), label = p_team), position = position_dodge(width = 0.9), vjust = -0.1, size = 5)

ggsave(path, width = 11.8, height = 7, dpi = 300)

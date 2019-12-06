library(ggplot2)


#######################
# set working directory
setwd("/home/tcigene/jade/ncounter/data/1205/")

############################
# get the ncounter table csv
ibd = read.csv("20191204_鳳梨釋迦萃取抗老化基因檢測V1.csv")
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

p = matrix(NA, nrow = 171, ncol = (dim(nor_data)[2] / 2))

rownames(p) = ibd$gene
colnames(p) = sample_name[2:length(sample_name)]
  
for (col in 1:dim(p)[2]){
  
  for (row in 1:171){
    
    p[row, col] = t.test(c(mock[row,1], mock[row,2] + 0.00001), c(nor_data[row, (col*2 - 1)], nor_data[row, (col*2)] + 0.00001), alternative = "two.sided", paired = F)$p.value
    
  }
}


#######################
# Sort by genetic group

# 抗氧化
g1 = c('SOD1','SOD2','GPX1','CAT')

# 抗老
g2_1 = c('CCT2','CCT5','CCT6A','CCT7','CCT8')
g2_2 = c('Pink1', 'Atg1','Atg8','SIRT1','FOXO')
g2_3 = c('PARP1','PARP2','NADSYN','MRPS5','Ubl-5','SOD3')
g2 = c('CCT2','CCT5','CCT6A','CCT7','CCT8','Pink1', 'Atg1','Atg8','SIRT1','FOXO','PARP1','PARP2','NADSYN','MRPS5','Ubl-5','SOD3') #,'Parkin'



# DNA修復
g3 = c('UNG','OGG1','MPG','APEX1','ERCC1','ERCC6','XPA','XRCC1','XRCC5','MSH2','MLH1','MSH6')

# 免疫
g4 = c('IL-1B','IL-8','IL-6','IL-10','IL-18','TNF-a')

# 美白-抗黑色素生成
g5 = c('TYR','TYRP1','MC1R','MITF')

# 膠原蛋白 彈力蛋白 玻尿酸 合成/降解
g6 = c('COL1A1','COL1A2','COL4A1','COL4A4','COL4A5','MMP1','MMP9','MMP2','TIMP1','ELN','FBN1','LOX','HAS2','HAS3')

# 抗發炎
g7 = c('IL-1B','IL-8','IL-6','IL-10','IL-18','TNF-a','IL-16','IL23','IL12A','IFNG','TGFB','IL3','IL4')

# 心血管保健
g8 = c('PTGIS','NOS3','EDN1','PLAT','PROC','VWF','F3','SERPINE1','PDGFC','FGF2','IGF2BP3','IGF1R','IL-8','IL-6','ICAM1','VCAM1','CASP8')

# 皮膚角質保濕
g9 = c('Tgm1','Krt1','Keratin 10','Keratin 14','AQP3','FLG-F','SMPD1','GBA','HAS2','HAS3')

# 脂肪肝
g10 = c('SREBP-1c (SREBF1)','PPAR-g','PPAR-a','SCD1 (SCD)','ACC (ACACA)')

# 提升HDL
g11 = c('CETP','SCARB1','apoA-I (APA1)','LDLR','ABCA1')

# 健髮
g12 = c('SRD5A1','SRD5A2','AR','KROX20','SCF','VEGFA','IGF1','TGFB','BDNF')

# 端粒酶活性
g13 = c('TERT','TERC','RTEL1')

# 呼吸道過敏
g14 = c('CD40','ERBB2','LIF','MALT1','NCK1','PAF1','DYNLL2','GRK5','PSMD4','RDH10','RELB','SCARF1','TNFSF14','ABR','IL13','IL4R','IL5RA','RELA')


g = c(g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11, g12, g13, g14)




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


  ################
  ################
  # group bar plot

    # Selected gene(kick out not significant gene )
    g2_1 = c('CCT2','CCT5','CCT6A','CCT7','CCT8')
    g2_2 = c('Pink1', 'Atg1','Atg8','SIRT1')
    g2_3 = c('PARP1','PARP2','NADSYN','MRPS5','Ubl-5') 
  
   
    # concentration 1
    folder = "/home/tcigene/jade/ncounter/result/plot/2579/"
    plant_name = "鳳梨釋迦萃取"
    team_name = c("Mock", "0.25 mg/ml 24h", "0.25 mg/ml 48h","0.5 mg/ml 24h", "0.5 mg/ml 48h")
    team = c("mock_nor", "X2579L.24h1", "X2579L.48h1", "X2579H.24h1", "X2579H.48h1")
    color = c("#ABB2B9", "#F24D98", "#59D044", "#F3A002", "#F2F44D")
    team_num = length(team)
    
    # concentration 2
    #team_name2 = c("Mock", "17% 蕎麥種皮萃取液 0.5% 24h", "17% 蕎麥種皮萃取液 0.5% 48h", "微脂粒 17% 蕎麥種皮 0.5% 24h", "微脂粒 17% 蕎麥種皮 0.5% 48h")
    #team2 = c("mock_nor", "X2603H.24h1", "X2603H.48h1", "X2604H.24h1", "X2604H.48h1")
    #color2 = c("#78909C", "#F8BBD0", "#D32F2F", "#90CAF9", "#1565C0")
    #team_num2 = length(team2)
    
    
    
    
    # sample ex:2649
    sample = mean[,team]
    sample_std = std[,team[2:length(team)]]
    sample_p = p_star[,team[2:length(team)]]
    
    
    
    
    ###############
    # G2: 抗老 g2_1
    
    # set variables (1) gene_group (2) team_number (3) team_name
    group_name = "抗老1"
    group = g2_1
    gene_num = length(group)
    
    
    # get team mean
    value = sample[group,]
    
    avg = matrix(NA, nrow = dim(value)[1] * dim(value)[2], ncol = 1)
    
    i = 1
    
    for (row in 1:dim(value)[1]){
      for(col in 1:dim(value)[2]){
        avg[i, 1] = value[row, col]
        i = i + 1
      }
    }
    
    # get team std
    value_std = sample_std[group,]
    mock_std = rep(0, rep(dim(value_std)[1]))
    
    value_std = cbind(mock_std, value_std)
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
    value_p = sample_p[group,]
    mock_p = rep("", rep(dim(value_p)[1]))
    
    value_p = cbind(mock_p, value_p)
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
    
    color = c("#ABB2B9", "#2e9fff", "#ffcd00", "#ae63e4", "#ff3c41")
    
    ggplot(data, aes(fill = team, y = avg, x = gene)) + 
      geom_bar(position = "dodge", stat = "identity") +
      scale_fill_manual(values = color) +
      geom_errorbar(aes(ymin = avg, ymax = avg + sd), width = .2, position = position_dodge(.9)) +
      theme(text = element_text(size = 20), plot.title = element_text(hjust = 0.5)) +
      xlab("Gene") + 
      ylab("Relative expression ratio") + 
      ggtitle(paste(plant_name, "-", group_name, spe = "")) +
      geom_text(aes(y = (avg + sd), label = p_team), position = position_dodge(width = 0.9), vjust = -0.1, size = 5)
    
    
    ggsave(path)
    
    
    ###########
    # 抗老 g2_2
    
    # set variables (1) gene_group (2) team_number (3) team_name
    group_name = "抗老2"
    group = g2_2
    gene_num = length(group)
    
    
    # get team mean
    value = sample[group,]
    
    avg = matrix(NA, nrow = dim(value)[1] * dim(value)[2], ncol = 1)
    
    i = 1
    
    for (row in 1:dim(value)[1]){
      for(col in 1:dim(value)[2]){
        avg[i, 1] = value[row, col]
        i = i + 1
      }
    }
    
    # get team std
    value_std = sample_std[group,]
    mock_std = rep(0, rep(dim(value_std)[1]))
    
    value_std = cbind(mock_std, value_std)
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
    value_p = sample_p[group,]
    mock_p = rep("", rep(dim(value_p)[1]))
    
    value_p = cbind(mock_p, value_p)
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
    
    color = c("#ABB2B9", "#2e9fff", "#ffcd00", "#ae63e4", "#ff3c41")
    
    ggplot(data, aes(fill = team, y = avg, x = gene)) + 
      geom_bar(position = "dodge", stat = "identity") +
      scale_fill_manual(values = color) +
      geom_errorbar(aes(ymin = avg, ymax = avg + sd), width = .2, position = position_dodge(.9)) +
      theme(text = element_text(size = 20), plot.title = element_text(hjust = 0.5)) +
      xlab("Gene") + 
      ylab("Relative expression ratio") + 
      ggtitle(paste(plant_name, "-", group_name, spe = "")) +
      geom_text(aes(y = (avg + sd), label = p_team), position = position_dodge(width = 0.9), vjust = -0.1, size = 5)
    
    
    ggsave(path)
    
    
    ###########
    # 抗老 g2_3
    
    # set variables (1) gene_group (2) team_number (3) team_name
    group_name = "抗老3"
    group = g2_3
    gene_num = length(group)
    
    
    # get team mean
    value = sample[group,]
    
    avg = matrix(NA, nrow = dim(value)[1] * dim(value)[2], ncol = 1)
    
    i = 1
    
    for (row in 1:dim(value)[1]){
      for(col in 1:dim(value)[2]){
        avg[i, 1] = value[row, col]
        i = i + 1
      }
    }
    
    # get team std
    value_std = sample_std[group,]
    mock_std = rep(0, rep(dim(value_std)[1]))
    
    value_std = cbind(mock_std, value_std)
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
    value_p = sample_p[group,]
    mock_p = rep("", rep(dim(value_p)[1]))
    
    value_p = cbind(mock_p, value_p)
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
    
    color = c("#ABB2B9", "#2e9fff", "#ffcd00", "#ae63e4", "#ff3c41")
    
    ggplot(data, aes(fill = team, y = avg, x = gene)) + 
      geom_bar(position = "dodge", stat = "identity") +
      scale_fill_manual(values = color) +
      geom_errorbar(aes(ymin = avg, ymax = avg + sd), width = .2, position = position_dodge(.9)) +
      theme(text = element_text(size = 20), plot.title = element_text(hjust = 0.5)) +
      xlab("Gene") + 
      ylab("Relative expression ratio") + 
      ggtitle(paste(plant_name, "-", group_name, spe = "")) +
      geom_text(aes(y = (avg + sd), label = p_team), position = position_dodge(width = 0.9), vjust = -0.1, size = 5)
    
    
    ggsave(path)
    
    
library("readxl")
library("openxlsx")
library("hash")
library("ggplot2")
library("stringr")
library(UpSetR)
library("dplyr")
library(plotly)
# packageVersion('plotly')
# library(htmlwidgets)
# source("http://peterhaschke.com/Code/multiplot.R")

bulkInfo <- function(dt_df, lipid_cl){
  dt_df$Bulk_c <- sub(":.*\\)", '', sub("[a-zA-Z]+[0-9]*[a-zA-Z]*\\([a-zA-Z\\-]*", '', dt_df$Bulk)) 
  dt_df$Bulk_l <- sub("\\-*[0-9]+:[0-9]+\\+*O*\\)", '', sub("[a-zA-Z]+[0-9]*[a-zA-Z]*\\(", '', dt_df$Bulk))
  dt_df$Bulk_l <- sub("[0-9]+\\:[0-9]\\)", 'A', dt_df$Bulk_l)
  #dt_df$Bulk_l <- sub("[a-zA-Z]+\\(", '', dt_df$Bulk)
  dt_df$Bulk_db <- sub("\\)", '', sub(".*:", '', dt_df$Bulk))
  
  # This section needs to change in different situations because sometimes we may have different name
  if (lipid_cl == 'HexCer' | lipid_cl == 'Cer'){
    pa_s <- "[a-zA-Z]+[0-9]*[a-zA-Z]*\\("
    pa_s2 <- "[a-zA-Z]+[0-9]*[a-zA-Z]*\\([a-zA-Z]*[0-9]+\\:[0-9]+_"
    dt_df$fa1 <- dt_df$Discrete
    dt_df$fa1 <- sub("_.*\\)", '', sub(pa_s, '', dt_df$fa1))
    dt_df$fa2 <- dt_df$Discrete
    dt_df$fa2 <- sub("\\+*O*\\)", '', sub(pa_s2, '', dt_df$fa2))
    dt_df$fa2_db <- sub("[0-9]+\\:", '', dt_df$fa2)
    dt_df$fa2_c <- sub("\\:[0-9]+", '', dt_df$fa2)
    
    dt_df$fa1_db <- sub("[a-z]*[0-9]+\\:", '', dt_df$fa1)
    #dt_df$fa1_c <- sub("[A-Z]*\\-*", "", sub("\\:[0-9]+", '', dt_df$fa1))
    dt_df$fa1_c <- sub('[a-z]*', '', sub("\\:[0-9]+", '', dt_df$fa1))
  }else if (lipid_cl == 'All'){
    dt_df$Class <- sub("\\(.+\\)", '', dt_df$Bulk)
  } else if (lipid_cl == 'TG'){
    pa_s <- paste(lipid_cl,  "\\(", sep = '')
    pa_s2 <- paste(lipid_cl,  ".*\\([0-9][0-9]*\\:[0-9][0-9]*_", sep = '')
    pa_s3 <- paste(lipid_cl,  ".*_", sep = '')
    dt_df$fa1 <- dt_df$Discrete
    dt_df$fa1 <- sub("_.*\\)", '', sub(pa_s, '', dt_df$fa1))
    dt_df$fa2 <- dt_df$Discrete
    dt_df$fa2 <- sub("_[0-9][0-9]*\\:[0-9][0-9]*\\)", '', sub(pa_s2, '', dt_df$fa2))
    dt_df$fa2_db <- sub("[0-9]+\\:", '', dt_df$fa2)
    dt_df$fa2_c <- sub("\\:[0-9]+", '', dt_df$fa2)
    
    dt_df$fa3 <- dt_df$Discrete
    dt_df$fa3 <- sub("\\)", '', sub(pa_s3, '', dt_df$fa3))
    dt_df$fa3_db <- sub("[0-9]+\\:", '', dt_df$fa3)
    dt_df$fa3_c <- sub("\\:[0-9]+", '', dt_df$fa3)
    
    
    dt_df$fa1_db <- sub("[A-Z]*\\-*[0-9]+\\:", '', dt_df$fa1)
    dt_df$fa1_c <- sub("[A-Z]*\\-*", "", sub("\\:[0-9]+", '', dt_df$fa1))
    dt_df$fa1 <- sub("[A-Z]*\\-*", "", dt_df$fa1)
    #dt_df$fa1_c <- sub("\\:[0-9]+", '', dt_df$fa1)
  }else if(lipid_cl %in% c('LPC', 'LPE')){
    pa_s <- paste(lipid_cl,  "\\(", sep = '')
    pa_s2 <- paste(lipid_cl,  ".*_", sep = '')
    dt_df$fa1 <- dt_df$Discrete
    dt_df$fa1 <- sub("\\)", '', sub(pa_s, '', dt_df$fa1))
    
    dt_df$fa1_db <- sub("[A-Z]*\\-*[0-9]+\\:", '', dt_df$fa1)
    dt_df$fa1_c <- sub("[A-Z]*\\-*", "", sub("\\:[0-9]+\\)", '', dt_df$fa1))
    dt_df$fa1 <- sub("[A-Z]*\\-*", "", dt_df$fa1)
    #dt_df$fa1_c <- sub("\\:[0-9]+", '', dt_df$fa1)
  }else{
    pa_s <- paste(lipid_cl,  "\\(", sep = '')
    pa_s2 <- paste(lipid_cl,  ".*_", sep = '')
    dt_df$fa1 <- dt_df$Discrete
    dt_df$fa1 <- sub("_.*\\)", '', sub(pa_s, '', dt_df$fa1))
    dt_df$fa2 <- dt_df$Discrete
    dt_df$fa2 <- sub("\\)", '', sub(pa_s2, '', dt_df$fa2))
    dt_df$fa2_db <- sub("[0-9]+\\:", '', dt_df$fa2)
    dt_df$fa2_c <- sub("\\:[0-9]+", '', dt_df$fa2)
    
    dt_df$fa1_db <- sub("[A-Z]*\\-*[0-9]+\\:", '', dt_df$fa1)
    dt_df$fa1_c <- sub("[A-Z]*\\-*", "", sub("\\:[0-9]+", '', dt_df$fa1))
    dt_df$fa1 <- sub("[A-Z]*\\-*", "", dt_df$fa1)
    #dt_df$fa1_c <- sub("\\:[0-9]+", '', dt_df$fa1)
  }
  
  return(dt_df)
  
}

calc_elem_n_f <- function(mass, ion, charge){
  element_df <- hash()
  element_df[['h']] <- 1.0078250321
  element_df[['c']] <- 12.0
  element_df[['n']] <- 14.0030740052
  element_df[['o']] <- 15.9949146221
  element_df[['p']] <- 30.97376151
  if (ion == '[M+H]+'){
    if (charge == "[M+HCOO]-"){
      ion_mass <- mass + element_df[['c']] + 2*element_df[['o']]
    }
  }else if (ion == '[M-H]-'){
    if (charge == '[M+H]+'){
      ion_mass <- mass + 2 * element_df[['h']]
    }
  }
  return(ion_mass)
}

col_palet <- function(i){
  library("scales")
  l_col <- scales::hue_pal()(i)
  return(l_col)
}

kmd_calculation <- function(d_df, elem){
  element_df <- hash()
  element_df[['h']] <- c(1.0078250321, 1)
  element_df[['c']] <- c(12.0, 12)
  element_df[['n']] <- c(14.0030740052, 14)
  element_df[['o']] <- c(15.9949146221, 16)
  element_df[['p']] <- c(30.97376151, 31)
  d_df$C <-  ''
  d_df$C <- str_match(d_df$Formula_Neutral, "C(\\d+)")[,2]
  d_df$C[is.na(d_df$C)] <- 0
  d_df$H <- ''
  d_df$H <- str_match(d_df$Formula_Neutral, "H(\\d+)")[,2]
  d_df$H[is.na(d_df$H)] <- 0
  d_df$O <- ''
  d_df$O <- str_match(d_df$Formula_Neutral, "O(\\d+)")[,2]
  d_df$O[is.na(d_df$O)] <- 0
  d_df$P <- ''
  d_df$P <- str_match(d_df$Formula_Neutral, "P(\\d+)")[,2]
  d_df$P[is.na(d_df$P)] <- 0
  d_df$N <- str_match(d_df$Formula_Neutral, "N(\\d+)")[,2]
  d_df$N[is.na(d_df$N)] <- 0
  #if (str_detect(d_df$Formula_Neutral, 'C')){
  #  d_df$C <- str_match(d_df$Formula_Neutral, "C(\\d+)")[,2]
  #}
  d_df$theo <- (as.integer(d_df$C)*element_df[['c']][1])+(as.integer(d_df$H)*element_df[['h']][1]) +
    (as.integer(d_df$O) * element_df[['o']][1]) + (as.integer(d_df$P) * element_df[['p']][1]) +
    (as.integer(d_df$N) * element_df[['n']][1])
  d_df$nominal <- (as.integer(d_df$C)*element_df[['c']][2])+(as.integer(d_df$H)*element_df[['h']][2]) +
    (as.integer(d_df$O) * element_df[['o']][2]) + (as.integer(d_df$P) * element_df[['p']][2]) +
    (as.integer(d_df$N) * element_df[['n']][2])
  km <- paste("km(", elem, ')', sep = '')
  d_df$V1 <- d_df$theo*(element_df[[elem]][2]/element_df[[elem]][1])
  kmd <- paste("kmd(", elem, ')', sep = '')
  d_df$V2 <- d_df$nominal - d_df$V1
  colnames(d_df)[which(names(d_df) == "V1")] <- km
  colnames(d_df)[which(names(d_df) == "V2")] <- kmd
  return(d_df)
}

####################################### LPC ##############################
#
# lpc boxplot ang ggplot info


cols<- c("0" = "#F8766D", "1" = "#E38900", "2" = "darkgoldenrod1", "3" = "#99A800", 
         "4" = "aquamarine3", "5" = "#00BC56", "6" = "dodgerblue3", "7" = "#00BFC4", 
         "8" = "darkslateblue", "9" ="#06A4FF", "10"=  "#A58AFF", "11"="hotpink3",
         "12"="#FB61D7", "13"="#FF66A8")

wb <- loadWorkbook("./AdiposeAtlas_lipidIdentification_05062020.xlsx")
data_info <- readWorkbook(wb, sheet=1)
lpc_df <- data_info[data_info$Class == "LPC", ]
lpc_df <- bulkInfo(lpc_df, 'LPC')
lpc_df <- kmd_calculation(lpc_df, 'h')
lpc_df <- kmd_calculation(lpc_df, 'o')

ggplot(data = lpc_df, mapping = aes(x =`RT_P` , y = `Exact_mass`, color = Bulk_c))  + geom_point(aes(shape=Bulk_l), size = 4) +
  scale_shape_manual(values = c(19, 17, 15))+
  labs(title = "LPC RT information", x = "RT", y = "Exact mass")  +
  theme_bw()

ggplot(data = lpc_df, mapping = aes(x = `RT_P`, y = `Exact_mass`, color = Bulk_db))  + geom_point(aes(shape=Bulk_l), size = 4) +
  scale_shape_manual(values = c(19, 17, 15))+
  labs(title = "LPC RT information", x = "RT", y = "Exact mass")  +
  theme_bw()

ggplot(data = lpc_df[lpc_df$Bulk_l == "", ], mapping = aes(x = `RT_P`, y = `kmd(h)`, color = Bulk_db))  + geom_point(aes(shape=Bulk_l), size = 4) +
  scale_shape_manual(values = c(19, 17, 15))+
  labs(title = "LPC RT information", x = "RT", y = "Exact mass")  +
  theme_bw()

a <- range(lpc_df[lpc_df$Bulk_l =='', 'kmd(h)'])
b <- range(lpc_df[lpc_df$Bulk_l =='', 'Bulk_c'])
scale_factor <- diff(as.numeric(a))/diff(as.numeric(b))
trans <- ~((. -as.numeric(a[1]))/ scale_factor) + as.numeric(b[1])
lpc_df[lpc_df$Bulk_l == '', 'Bulk_l'] <- 'A'

ggplot(data = lpc_df, mapping = aes(x = RT_P, y = `kmd(h)`, color = Bulk_db))  + geom_point(aes(shape=Bulk_l), size = 4) +
  scale_shape_manual(values = c(19,17, 15))+
  scale_color_manual(values = cols, aesthetics = c("colour", "fill"),  breaks= c(0, 1,2,3,4,5,6,7,8,9,10,11,12,13)) +
  scale_y_continuous(sec.axis = sec_axis(trans=trans, name="# number of carbons", breaks = c(14, 15, 16, 17, 18, 19, 20, 22, 24))) +
  labs(title = "LPC RT information", x = "Retention Time", y = "KMD(H)") +
  theme_bw()


####################################### LPE ##############################
#
# lpe ggplot info
#
lpe_df <- data_info[data_info$Class == "LPE", ]
lpe_df <- bulkInfo(lpe_df, 'LPE')
lpe_df <- kmd_calculation(lpe_df, 'h')
lpe_df <- kmd_calculation(lpe_df, 'o')



a <- range(lpe_df[lpe_df$Bulk_l =='', 'kmd(h)'])
b <- range(lpe_df[lpe_df$Bulk_l =='', 'Bulk_c'])
scale_factor <- diff(as.numeric(a))/diff(as.numeric(b))
trans <- ~((. -as.numeric(a[1]))/ scale_factor) + as.numeric(b[1])
lpe_df[lpe_df$Bulk_l == '','Bulk_l'] <- 'A'

ggplot(data = lpe_df, mapping = aes(x = RT_P, y = `kmd(h)`, color = Bulk_db))  + geom_point(aes(shape=Bulk_l), size = 4) +
  scale_shape_manual(values = c(19,17, 15))+
  scale_color_manual(values = cols, aesthetics = c("colour", "fill"),  breaks= c(0, 1,2,3,4,5,6,7,8,9,10,11,12,13)) +
  scale_y_continuous(sec.axis = sec_axis(trans=trans, name="# number of carbons", breaks = c(16, 17, 18, 20, 22))) +
  labs(title = "LPE RT information", x = "Retention Time", y = "KMD(H)") +
  theme_bw()

################################################################ PC ################################################################
#
#
#


pc_df <- data_info[data_info['Class'] == 'PC', ]
pc_df <- bulkInfo(pc_df, 'PC')
pc_df <- fa_all_info(pc_df)

pc_df <- kmd_calculation(pc_df, 'h')
pc_df <- kmd_calculation(pc_df, 'o')


ggplot(data = pc_df, mapping = aes(x = RT_P, y = `Exact_mass`, color = Bulk_db))  + geom_point(aes(shape= Bulk_l), size =4) +
  facet_wrap(vars(Fa_a))+
  labs(title = "PC information", x = "Retention Time", y = "Exact mass") +
  theme_bw()

ggplot(data = pc_df, mapping = aes(x = RT_P, y = `Exact_mass`, color = Bulk_c))  + geom_point(aes(shape= Bulk_l), size =4) +
  labs(title = "PC information", x = "Retention Time", y = "Exact mass") +
  theme_bw()

a <- range(pc_df[pc_df$Bulk_l =='', 'kmd(h)'])
b <- range(pc_df[pc_df$Bulk_l =='', 'Bulk_c'])
scale_factor <- diff(as.numeric(a))/diff(as.numeric(b))
trans <- ~((. -as.numeric(a[1]))/ scale_factor) + as.numeric(b[1])

ggplot(data = pc_df, mapping = aes(x = RT_P, y = `kmd(h)`, color = Bulk_db, shape= Bulk_l))  + geom_point(size = 4) +
#ggplot(data = pc_df[pc_df$Bulk_l == "", ], mapping = aes(x = RT_P, y = `kmd(h)`, color = Bulk_db))  + geom_point(size = 4) +
  scale_shape_manual(values = c(19,17, 15))+
  scale_color_manual(values = cols, aesthetics = c("colour", "fill"),  breaks= c(0, 1,2,3,4,5,6,7,8,9,10,11,12,13)) +
  scale_y_continuous(sec.axis = sec_axis(trans=trans, name="# number of carbons", breaks = c(28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 40, 44))) +
  labs(title = "PC RT information", x = "Retention Time", y = "KMD(H)") +
  theme_bw()

a <- range(pc_df[pc_df$Bulk_l =='P', 'kmd(h)'])
b <- range(pc_df[pc_df$Bulk_l =='P', 'Bulk_c'])
scale_factor <- diff(as.numeric(a))/diff(as.numeric(b))
trans <- ~((. -as.numeric(a[1]))/ scale_factor) + as.numeric(b[1])


ggplot(data = pc_df[pc_df$Bulk_l %in% c('P', 'O'), ], mapping = aes(x = RT_P, y = `kmd(h)`, color = Bulk_db))  + geom_point(aes(shape=Bulk_l), size = 4) +
  scale_shape_manual(values = c(17, 15))+
  scale_color_manual(values = cols, aesthetics = c("colour", "fill"),  breaks= c(0, 1,2,3,4,5,6,7,8,9,10,11,12,13)) +
  scale_y_continuous(sec.axis = sec_axis(trans=trans, name="# number of carbons", breaks = c( 30,32, 34, 35, 36, 37, 38, 44))) +
  labs(title = "PC RT information", x = "Retention Time", y = "KMD(H)") +
  theme_bw()


ggplot(data = pc_df, mapping = aes(x = RT_P, y = `kmd(o)`, color = Bulk_c))  + geom_point(aes(shape=Bulk_l), size = 4) +
  scale_shape_manual(values = c(19,17, 15))+
  labs(title = "PC RT information", x = "Retention Time", y = "KMD(O)") +
  theme_bw()



############################################################ PE ######################################


pe_df <- data_info[data_info['Class'] == 'PE', ]
pe_df <- bulkInfo(pe_df, 'PE')
pe_df <- fa_all_info(pe_df)

pe_df <- kmd_calculation(pe_df, 'h')
pe_df <- kmd_calculation(pe_df, 'o')

ggplot(data = pe_df, mapping = aes(x = RT_P, y = `Exact_mass`, color = Bulk_c))  + geom_point(aes(shape= Bulk_l), size =4) +
  facet_wrap(vars(Fa_a))+
  labs(title = "PE information", x = "Retention Time", y = "Exact mass") +
  theme_bw()

a <- range(pe_df[pe_df$Bulk_l =='', 'kmd(h)'])
b <- range(pe_df[pe_df$Bulk_l =='', 'Bulk_c'])
scale_factor <- diff(as.numeric(a))/diff(as.numeric(b))
trans <- ~((. -as.numeric(a[1]))/ scale_factor) + as.numeric(b[1])


#ggplot(data = pe_df[pe_df$Bulk_l == "", ], mapping = aes(x = RT_P, y = `kmd(h)`, color = Bulk_db))  + geom_point(aes(shape=Bulk_l), size = 4) +
ggplot(data = pe_df[pe_df$Bulk_l == "", ], mapping = aes(x = RT_P, y = `kmd(h)`, color = Bulk_db))  + geom_point(size = 4) +
  scale_shape_manual(values = c(19,17, 15))+
  scale_color_manual(values = cols, aesthetics = c("colour", "fill"),  breaks= c(0, 1,2,3,4,5,6,7,8,9,10,11,12,13)) +
  scale_y_continuous(sec.axis = sec_axis(trans=trans, name="# number of carbons", breaks = c(32, 34, 35, 36, 37, 38, 40, 42))) +
  labs(title = "PE RT information", x = "Retention Time", y = "KMD(H)") +
  theme_bw()



a <- range(pe_df[pe_df$Bulk_l =='P', 'kmd(h)'])
b <- range(pe_df[pe_df$Bulk_l =='P', 'Bulk_c'])
scale_factor <- diff(as.numeric(a))/diff(as.numeric(b))
trans <- ~((. -as.numeric(a[1]))/ scale_factor) + as.numeric(b[1])

ggplot(data = pe_df[pe_df$Bulk_l %in% c('O','P'), ], mapping = aes(x = RT_P, y = `kmd(h)`, color = Bulk_db, shape=Bulk_l))  + geom_point(size = 4) +
  scale_shape_manual(values = c(17, 15))+
  scale_color_manual(values = cols, aesthetics = c("colour", "fill"),  breaks= c(0, 1,2,3,4,5,6,7,8,9,10,11,12,13)) +
  scale_y_continuous(sec.axis = sec_axis(trans=trans, name="# number of carbons", breaks = c(32, 33, 34, 35, 36, 37, 38,39, 40, 42))) +
  labs(title = "PE RT information", x = "Retention Time", y = "KMD(H)") +
  theme_bw()

################################################# PG ########################################################
#
#
#

pg_df <- data_info[data_info['Class'] == 'PG', ]
pg_df <- bulkInfo(pg_df, 'PG')
pg_df <- fa_all_info(pg_df)

pg_df <- kmd_calculation(pg_df, 'h')
pg_df <- kmd_calculation(pg_df, 'o')

ggplot(data = pg_df, mapping = aes(x = RT_P, y = `Exact_mass`, color = Bulk_db))  + geom_point(aes(shape= Bulk_l), size =4) +
  labs(title = "PG information", x = "Retention Time", y = "Exact mass") +
  theme_bw()

a <- range(pg_df[pg_df$Bulk_l =='', 'kmd(h)'])
b <- range(pg_df[pg_df$Bulk_l =='', 'Bulk_c'])
scale_factor <- diff(as.numeric(a))/diff(as.numeric(b))
trans <- ~((. -as.numeric(a[1]))/ scale_factor) + as.numeric(b[1])


ggplot(data = pg_df, mapping = aes(x = RT_P, y = `kmd(h)`, color = Bulk_db))  + geom_point( size = 4) +
  scale_shape_manual(values = c(19,17, 15))+
  scale_color_manual(values = cols, aesthetics = c("colour", "fill"),  breaks= c(0, 1,2,3,4,5,6,7,8,9,10,11,12,13)) +
  scale_y_continuous(sec.axis = sec_axis(trans=trans, name="# number of carbons", breaks = c(34, 36, 38))) +
  labs(title = "PG RT information", x = "Retention Time", y = "KMD(H)") +
  theme_bw()







################################################# PI ########################################################
#
#
#

pi_df <- data_info[data_info$Class == 'PI',]
pi_df <- bulkInfo(pi_df, 'PI')
pi_df <- fa_all_info(pi_df)
pi_df <- kmd_calculation(pi_df, 'h')

ggplot(data = pi_df, mapping = aes(x = RT_P, y = `Exact_mass`, color = Bulk_c))  + geom_point(size = 4) +
  labs(title = "PI RT information", x = "Retention Time", y = "Exact mass") +
  theme_bw()

ggplot(data = pi_df, mapping = aes(x = RT_P, y = `Exact_mass`, color = Bulk_db))  + geom_point(size = 4) +
  facet_wrap(vars(Fa_a))+
  labs(title = "PI RT information", x = "Retention Time", y = "Exact mass") +
  theme_bw()

a <- range(pi_df[pi_df$Bulk_l =='', 'kmd(h)'])
b <- range(pi_df[pi_df$Bulk_l =='', 'Bulk_c'])
scale_factor <- diff(as.numeric(a))/diff(as.numeric(b))
trans <- ~((. -as.numeric(a[1]))/ scale_factor) + as.numeric(b[1])


ggplot(data = pi_df, mapping = aes(x = RT_P, y = `kmd(h)`, color = Bulk_db))  + geom_point( size = 4) +
  scale_y_continuous(sec.axis = sec_axis(trans=trans, name="# number of carbons", breaks = c(32, 34, 36, 37,  38, 40))) +
  scale_color_manual(values = cols, aesthetics = c("colour", "fill"),  breaks= c(0, 1,2,3,4,5,6,7,8,9,10,11,12,13)) +
  labs(title = "PI RT information", x = "Retention Time", y = "KMD(H)") +
  theme_bw()


################################################################# PS ########################################
#
#
#


ps_df <- data_info[data_info$Class == 'PS',]
ps_df <- bulkInfo(ps_df, 'PS')
ps_df <- kmd_calculation(ps_df, 'h')


ggplot(data = ps_df, mapping = aes(x = RT_P, y = Exact_mass, color = Bulk_c))  + geom_point(size = 4) +
  labs(title = "PS RT information", x = "Retention time", y = "Exact mass") +
  theme_bw()

a <- range(ps_df[ps_df$Bulk_l =='', 'kmd(h)'])
b <- range(ps_df[ps_df$Bulk_l =='', 'Bulk_c'])
scale_factor <- diff(as.numeric(a))/diff(as.numeric(b))
trans <- ~((. -as.numeric(a[1]))/ scale_factor) + as.numeric(b[1])


ggplot(data = ps_df, mapping = aes(x = RT_P, y = `kmd(h)`, color = Bulk_db))  + geom_point( size = 4) +
  scale_y_continuous(sec.axis = sec_axis(trans=trans, name="# number of carbons", breaks = c(34, 36, 38, 40))) +
  scale_color_manual(values = cols, aesthetics = c("colour", "fill"),  breaks= c(0, 1,2,3,4,5,6,7,8,9,10,11,12,13)) +
  labs(title = "PS RT information", x = "Retention Time", y = "KMD(H)") +
  theme_bw()


########################################## DG ####################################################


dg_df <- data_info[data_info$Class == 'DG',]
dg_df <- bulkInfo(dg_df, 'DG')
dg_df <- fa_all_info(dg_df)
dg_df <- kmd_calculation(dg_df, 'h')


ggplot(data =dg_df, mapping = aes(x = RT_P, y = `Exact_mass`, color = Bulk_c))  + geom_point(size =4) +
  labs(title = "DG RT Information", x = "Retention time", y = "Exact mass") +
  theme_bw()


ggplot(data =dg_df[dg_df$Fa_a == '18:1',], mapping = aes(x = RT_P, y = `Exact_mass`, color = Bulk_db))  + geom_point(size =4) +
  labs(title = "DG RT Information", x = "Retention time", y = "Exact mass") +
  theme_bw()
a <- range(dg_df[dg_df$Bulk_db ==2 & dg_df$Fa_a == '18:1', 'Exact_mass'])
b <- range(dg_df[dg_df$Bulk_db == 2 & dg_df$Fa_a == '18:1', 'Bulk_c'])
scale_factor <- diff(as.numeric(a))/diff(as.numeric(b))
trans <- ~((. -as.numeric(a[1]))/ scale_factor) + as.numeric(b[1])


ggplot(data =dg_df[dg_df$Fa_a == '18:1',], mapping = aes(x = RT_P, y = `Exact_mass`, color = Bulk_db))  + geom_point( size = 4) +
  scale_y_continuous(sec.axis = sec_axis(trans=trans, name="# number of carbons", breaks = c(26, 28, 30, 31, 32, 33, 34, 35, 36, 37,  38, 40))) +
  scale_color_manual(values = cols, aesthetics = c("colour", "fill"),  breaks= c(0, 1,2,3,4,5,6,7,8,9,10,11,12,13)) +
  labs(title = "DG RT information", x = "Retention Time", y = "Exact mass") +
  theme_bw()

ggplot(data = dg_df, mapping = aes(x = RT_P, y = Exact_mass, color = Bulk_db))  + geom_point(size =4) +
  facet_wrap(vars(Fa_a))+
  labs(title = "DG RT Information", x = "Retention time", y = "Exact_mass") +
  theme_bw()

a <- range(dg_df[dg_df$Bulk_l =='', 'kmd(h)'])
b <- range(dg_df[dg_df$Bulk_l =='', 'Bulk_c'])
scale_factor <- diff(as.numeric(a))/diff(as.numeric(b))
trans <- ~((. -as.numeric(a[1]))/ scale_factor) + as.numeric(b[1])


ggplot(data = dg_df, mapping = aes(x = RT_P, y = `kmd(h)`, color = Bulk_db))  + geom_point( size = 4) +
  scale_y_continuous(sec.axis = sec_axis(trans=trans, name="# number of carbons", breaks = c(26, 28, 30, 31, 32, 33, 34, 35, 36, 37,  38, 40))) +
  scale_color_manual(values = cols, aesthetics = c("colour", "fill"),  breaks= c(0, 1,2,3,4,5,6,7,8,9,10,11,12,13)) +
  labs(title = "DG RT information", x = "Retention Time", y = "KMD(H)") +
  theme_bw()


########################################## TG ####################################################


tg_df2 <- data_info[data_info$Class == 'TG', ]
tg_df2 <- bulkInfo(tg_df2, 'TG')
tg_df2 <- kmd_calculation(tg_df2, 'h')
tg_df2$fa1_fa2 <- paste(tg_df2$fa1, tg_df2$fa2, sep = ' ')
tg_df2$fa1_fa3 <- paste(tg_df2$fa1, tg_df2$fa3, sep = ' ')
tg_df2$fa2_fa3 <- paste(tg_df2$fa2, tg_df2$fa3, sep = ' ')


ggplot(data =tg_df2, mapping = aes(x = RT_P, y = `kmd(h)`, color = Bulk_c))  + geom_point(size =4) +
  labs(title = "TG RT Information RPC18", x = "Retention time", y = "KMD(H)") +
  theme_bw()
ggplot(data =tg_df2, mapping = aes(x = RT_IT, y = `Exact_mass`, color = Bulk_c))  + geom_point(size =4) +
  labs(title = "TG RT Information", x = "Retention time", y = "Exact mass") +
  theme_bw()

ggplot(data =tg_df2, mapping = aes(x = RT_U, y = `kmd(h)`, color = Bulk_c))  + geom_point(size =4) +
  labs(title = "TG RT Information RPC30", x = "Retention time", y = "KMD(H)") +
  theme_bw()

ggplot(data =tg_df2, mapping = aes(x = RT_A, y = `Exact_mass`, color = Bulk_c))  + geom_point(size =4) +
  labs(title = "TG RT Information", x = "Retention time", y = "Exact mass") +
  theme_bw()

ggplot(data =tg_df2, mapping = aes(x = RT_P, y = `Exact_mass`, color = Bulk_db))  + geom_point(size =4) +
  labs(title = "TG RT Information", x = "Retention time", y = "Exact mass") +
  theme_bw()
ggplot(data =tg_df2, mapping = aes(x = RT_IT, y = `Exact_mass`, color = Bulk_db))  + geom_point(size =4) +
  labs(title = "TG RT_IT Information", x = "Retention time", y = "Exact mass") +
  theme_bw()

ggplot(data =tg_df2, mapping = aes(x = RT_U, y = `Exact_mass`, color = Bulk_db))  + geom_point(size =4) +
  labs(title = "TG RT Information", x = "Retention time", y = "Exact mass") +
  theme_bw()

ggplot(data =tg_df2, mapping = aes(x = RT_A, y = `Exact_mass`, color = Bulk_db))  + geom_point(size =4) +
  labs(title = "TG RT_A Information", x = "Retention time", y = "Exact mass") +
  theme_bw()

tg_df2 <- kmd_calculation(tg_df2, 'h')

a <- range(tg_df2[!is.na(tg_df2$RT_P) ,'kmd(h)'])
b <- range(tg_df2[!is.na(tg_df2$RT_P), 'Bulk_c'])
scale_factor <- diff(as.numeric(a))/diff(as.numeric(b))
trans <- ~((. -as.numeric(a[1]))/ scale_factor) + as.numeric(b[1])


ggplot(data = tg_df2[!is.na(tg_df2$RT_P),], mapping = aes(x = RT_P, y = `kmd(h)`, color = Bulk_db))  + geom_point( size = 4) +
  scale_y_continuous(sec.axis = sec_axis(trans=trans, name="# number of carbons", breaks = c(24, 26, 28, 30, 34, 36,  38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 60, 62))) +
  scale_color_manual(values = cols, aesthetics = c("colour", "fill"),  breaks= c(0, 1,2,3,4,5,6,7,8,9,10,11,12,13)) +
  labs(title = "TG RT_P information", x = "Retention Time", y = "KMD(H)") +
  theme_bw()

a <- range(tg_df2[!is.na(tg_df2$RT_U) ,'kmd(h)'])
b <- range(tg_df2[!is.na(tg_df2$RT_U), 'Bulk_c'])
scale_factor <- diff(as.numeric(a))/diff(as.numeric(b))
trans <- ~((. -as.numeric(a[1]))/ scale_factor) + as.numeric(b[1])


ggplot(data = tg_df2[!is.na(tg_df2$RT_U),], mapping = aes(x = RT_U, y = `kmd(h)`, color = Bulk_db))  + geom_point( size = 4) +
  scale_y_continuous(sec.axis = sec_axis(trans=trans, name="# number of carbons", breaks = c(24, 26, 28, 34, 36,  38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 62))) +
  scale_color_manual(values = cols, aesthetics = c("colour", "fill"),  breaks= c(0, 1,2,3,4,5,6,7,8,9,10,11,12,13)) +
  labs(title = "TG RT_U information", x = "Retention Time", y = "KMD(H)") +
  theme_bw()

a <- range(tg_df2[!is.na(tg_df2$RT_A) ,'kmd(h)'])
b <- range(tg_df2[!is.na(tg_df2$RT_A), 'Bulk_c'])
scale_factor <- diff(as.numeric(a))/diff(as.numeric(b))
trans <- ~((. -as.numeric(a[1]))/ scale_factor) + as.numeric(b[1])


ggplot(data = tg_df2[!is.na(tg_df2$RT_A),], mapping = aes(x = RT_A, y = `kmd(h)`, color = Bulk_db))  + geom_point( size = 4) +
  scale_y_continuous(sec.axis = sec_axis(trans=trans, name="# number of carbons", breaks = c(24, 26,32, 34, 36,  38, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61,62))) +
  scale_color_manual(values = cols, aesthetics = c("colour", "fill"),  breaks= c(0, 1,2,3,4,5,6,7,8,9,10,11,12,13)) +
  labs(title = "TG RT_A information", x = "Retention Time", y = "KMD(H)") +
  theme_bw()

a <- range(tg_df2[!is.na(tg_df2$RT_IT) ,'kmd(h)'])
b <- range(tg_df2[!is.na(tg_df2$RT_IT), 'Bulk_c'])
scale_factor <- diff(as.numeric(a))/diff(as.numeric(b))
trans <- ~((. -as.numeric(a[1]))/ scale_factor) + as.numeric(b[1])


ggplot(data = tg_df2[!is.na(tg_df2$RT_IT),], mapping = aes(x = RT_IT, y = `kmd(h)`, color = Bulk_db))  + geom_point( size = 4) +
  scale_y_continuous(sec.axis = sec_axis(trans=trans, name="# number of carbons", breaks = c(24, 26, 28, 30, 36,  38, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 60, 61,62))) +
  scale_color_manual(values = cols, aesthetics = c("colour", "fill"),  breaks= c(0, 1,2,3,4,5,6,7,8,9,10,11,12,13)) +
  labs(title = "TG RT_IT information", x = "Retention Time", y = "KMD(H)") +
  theme_bw()


a <- range(tg_df2[!is.na(tg_df2$RT_A) ,'kmd(h)'])
b <- range(tg_df2[!is.na(tg_df2$RT_A), 'Bulk_c'])
scale_factor <- diff(as.numeric(a))/diff(as.numeric(b))
trans <- ~((. -as.numeric(a[1]))/ scale_factor) + as.numeric(b[1])


ggplot(data = tg_df2[!is.na(tg_df2$RT_A),], mapping = aes(x = RT_A, y = `kmd(h)`, color = Bulk_db))  + geom_point( size = 4) +
  scale_y_continuous(sec.axis = sec_axis(trans=trans, name="# number of carbons", breaks = c(24, 26, 28, 30, 34, 36,  38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 60, 62))) +
  scale_color_manual(values = cols, aesthetics = c("colour", "fill"),  breaks= c(0, 1,2,3,4,5,6,7,8,9,10,11,12,13)) +
  labs(title = "TG RT_P information", x = "Retention Time", y = "KMD(H)") +
  theme_bw()

f_i = "17:0"
f_i2 = '17:1 19:2'
fa18_df2 <- tg_df2[tg_df2$fa1 == f_i | tg_df2$fa2 == f_i | tg_df2$fa3 == f_i, ]
fa18_df2 <- tg_df2[tg_df2$fa1_fa2 == f_i2 | tg_df2$fa1_fa3 == f_i2 | tg_df2$fa2_fa3 == f_i2, ]
ggplot(data = fa18_df2, mapping = aes(x = RT_IT, y = `kmd(h)`, color = Bulk_db))  + geom_point( size = 4) +
  scale_color_manual(values = cols, aesthetics = c("colour", "fill"),  breaks= c(0, 1,2,3,4,5,6,7,8,9,10,11,12,13)) +
  labs(title = "TG RT_IT information", x = "Retention Time", y = "KMD(H)") +
  theme_bw()


# upset graph


tg_df2$RT_P[tg_df2$RT_P > 0] <- 1
tg_df2$RT_P[is.na(tg_df2$RT_P)] <- 0
tg_df2$RT_U[tg_df2$RT_U > 0] <- 1
tg_df2$RT_U[is.na(tg_df2$RT_U)] <- 0
tg_df2$RT_A[tg_df2$RT_A > 0] <- 1
tg_df2$RT_A[is.na(tg_df2$RT_A)] <- 0
tg_df2$RT_IT[tg_df2$RT_IT > 0] <- 1
tg_df2$RT_IT[is.na(tg_df2$RT_IT)] <- 0
tg_df <- tg_df2[c("Discrete", "RT_P", "RT_IT", "RT_U", "RT_A")]

upset(tg_df2, nsets = 4, number.angles = 30, point.size = 3.5, line.size = 2, 
      mainbar.y.label = "Genre Intersections", sets.x.label = "Discrete Lipids", 
      text.scale = c(2, 2, 1, 1, 2, 0.75))

upset(tg_df, nsets = 4, number.angles = 30, point.size = 3.5, line.size = 2, 
      mainbar.y.label = "Genre Intersections", sets.x.label = "Discrete Lipids")


tg_df <- tg_df2[c("Bulk", "RT_P", "RT_IT", "RT_U", "RT_A")]
tg_df <- tg_df[!duplicated(tg_df$Bulk), ]
upset(tg_df, nsets = 4, number.angles = 30, point.size = 3.5, line.size = 2, 
      mainbar.y.label = "Genre Intersections", sets.x.label = "Bulk Lipids")


################################################### All lipids Plot

ggplot(data = data_info, mapping = aes(x = RT_P, y = `Exact_mass`, color = Class))  + geom_point(size =4) +
  labs(title = "All Lipids RT Information", x = "Retention time", y = "Exact mass") +
  theme_bw()


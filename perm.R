library(vegan)
library(tidyverse)
library(reshape2)
library(scales)
otu <-
  read.table(
    file = "./metaphlan.tsv",
    sep = "\t",
    header = T,
    row.names = 1,
    quote = ""
  ) 
group <-
  read.table(
    file = "./group_info.tsv",
    sep = "\t",
    header = T,
    # row.names = 1,
    quote = ""
  ) %>% filter(.,!grepl("control", Group)) #remove control
#import clinical (without sex and age)
clinic_original <-
  readxl::read_xlsx("./clincal.xlsx",
                    sheet = 1)
#import sex and age
clinic_sexage <-
  readxl::read_xlsx("./clincal.xlsx",
                    sheet = 2)  %>% filter(!grepl("C.*", m0)) %>% pivot_longer(
                      cols = 1:3,
                      names_to = "Group",
                      values_to = "samples",
                      values_drop_na = T
                    )
path <-
  read.table(
    file = "./path.csv",
    sep = ",",
    header = T,
    row.names = 1,
    quote = ""
  ) %>% t()
ko <-
  read.table(
    file = "./ko.anno.tsv",
    sep = ",",
    header = T,
    row.names = 1,
    quote = ""
  ) %>% t()
#to a list
clinic <- list(
  TG = clinic_original[, 1:2] %>% fill(2),
  Chol = clinic_original[, c(1, 3)] %>% fill(2),
  HDL = clinic_original[, c(1, 4)] %>% fill(2),
  LDL = clinic_original[, c(1, 5)] %>% fill(2),
  ALB = clinic_original[, c(1, 6)] %>% fill(2),
  Cr = clinic_original[, c(1, 7)] %>% fill(2),
  CTnT = clinic_original[, c(1, 8)] %>% fill(2),
  BNP = clinic_original[, c(1, 9)] %>% fill(2),
  HbA1 = clinic_original[, c(1, 11)] %>% fill(2, .direction = "downup"),
  sex = clinic_sexage[, c(4, 1)] %>% mutate(gender = ifelse(sex == "male", 1, 2)) %>% select(c(-2)),
  age = clinic_sexage[, c(4, 2)]
)
#select species from otu
tax_s <-
  slice(otu, str_which(rownames(otu), ".*s__.*")) %>% select(., contains(group$Sample)) %>% t()
#select genus from otu
tax_g <-
  slice(otu, str_which(rownames(otu), ".*g__((?!s__).)*$")) %>% select(., contains(group$Sample)) %>% t()
#prepare a matrix 
adonis_result <-
  matrix(nrow = 11,
         ncol = 8,
         dimnames = list(
           names(clinic),
           c(
             "tax_s_R2",
             "tax_s_p",
             "tax_g_R2",
             "tax_g_p",
             "path_R2",
             "path_p",
             "ko_R2",
             "ko_p"
           )
         ))

{
  adonis_result[1, 1:2] <-
    adonis2(tax_s ~ TG,
            clinic[["TG"]],
            permutations = 999,
            distance = 'bray')[1, c(3, 5)] %>% as.matrix()
  adonis_result[2, 1:2] <-
    adonis2(tax_s ~ Chol,
            clinic[["Chol"]],
            permutations = 999,
            distance = 'bray')[1, c(3, 5)] %>% as.matrix()
  adonis_result[3, 1:2] <-
    adonis2(tax_s ~ HDL,
            clinic[["HDL"]],
            permutations = 999,
            distance = 'bray')[1, c(3, 5)] %>% as.matrix()
  adonis_result[4, 1:2] <-
    adonis2(tax_s ~ LDL,
            clinic[["LDL"]],
            permutations = 999,
            distance = 'bray')[1, c(3, 5)] %>% as.matrix()
  adonis_result[5, 1:2] <-
    adonis2(tax_s ~ ALB,
            clinic[["ALB"]],
            permutations = 999,
            distance = 'bray')[1, c(3, 5)] %>% as.matrix()
  adonis_result[6, 1:2] <-
    adonis2(tax_s ~ Cr,
            clinic[["Cr"]],
            permutations = 999,
            distance = 'bray')[1, c(3, 5)] %>% as.matrix()
  adonis_result[7, 1:2] <-
    adonis2(tax_s ~ cTnT,
            clinic[["CTnT"]],
            permutations = 999,
            distance = 'bray')[1, c(3, 5)] %>% as.matrix()
  adonis_result[8, 1:2] <-
    adonis2(tax_s ~ BNP,
            clinic[["BNP"]],
            permutations = 999,
            distance = 'bray')[1, c(3, 5)] %>% as.matrix()
  
  adonis_result[9, 1:2] <-
    adonis2(tax_s ~ HbA1c,
            clinic[["HbA1"]],
            permutations = 999,
            distance = 'bray')[1, c(3, 5)] %>% as.matrix()
  adonis_result[10, 1:2] <-
    adonis2(tax_s ~ gender,
            clinic[["sex"]][clinic[["sex"]]$samples %in% rownames(tax_s),],
            permutations = 999,
            distance = 'bray')[1, c(3, 5)] %>% as.matrix()
  adonis_result[11, 1:2] <-
    adonis2(tax_s ~ age,
            clinic[["age"]][clinic[["sex"]]$samples %in% rownames(tax_s),],
            permutations = 999,
            distance = 'bray')[1, c(3, 5)] %>% as.matrix()
  }
{
  adonis_result[1, 3:4] <-
    adonis2(tax_g ~ TG,
            clinic[["TG"]],
            permutations = 999,
            distance = 'bray')[1, c(3, 5)] %>% as.matrix()
  adonis_result[2, 3:4] <-
    adonis2(tax_g ~ Chol,
            clinic[["Chol"]],
            permutations = 999,
            distance = 'bray')[1, c(3, 5)] %>% as.matrix()
  adonis_result[3, 3:4] <-
    adonis2(tax_g ~ HDL,
            clinic[["HDL"]],
            permutations = 999,
            distance = 'bray')[1, c(3, 5)] %>% as.matrix()
  adonis_result[4, 3:4] <-
    adonis2(tax_g ~ LDL,
            clinic[["LDL"]],
            permutations = 999,
            distance = 'bray')[1, c(3, 5)] %>% as.matrix()
  adonis_result[5, 3:4] <-
    adonis2(tax_g ~ ALB,
            clinic[["ALB"]],
            permutations = 999,
            distance = 'bray')[1, c(3, 5)] %>% as.matrix()
  adonis_result[6, 3:4] <-
    adonis2(tax_g ~ Cr,
            clinic[["Cr"]],
            permutations = 999,
            distance = 'bray')[1, c(3, 5)] %>% as.matrix()
  adonis_result[7, 3:4] <-
    adonis2(tax_g ~ cTnT,
            clinic[["CTnT"]],
            permutations = 999,
            distance = 'bray')[1, c(3, 5)] %>% as.matrix()
  adonis_result[8, 3:4] <-
    adonis2(tax_g ~ BNP,
            clinic[["BNP"]],
            permutations = 999,
            distance = 'bray')[1, c(3, 5)] %>% as.matrix()
  
  adonis_result[9, 3:4] <-
    adonis2(tax_g ~ HbA1c,
            clinic[["HbA1"]],
            permutations = 999,
            distance = 'bray')[1, c(3, 5)] %>% as.matrix()
  adonis_result[10, 3:4] <-
    adonis2(tax_g ~ gender,
            clinic[["sex"]][clinic[["sex"]]$samples %in% rownames(tax_g),],
            permutations = 999,
            distance = 'bray')[1, c(3, 5)] %>% as.matrix()
  adonis_result[11, 3:4] <-
    adonis2(tax_g ~ age,
            clinic[["age"]][clinic[["sex"]]$samples %in% rownames(tax_g),],
            permutations = 999,
            distance = 'bray')[1, c(3, 5)] %>% as.matrix()
}
{
  adonis_result[1, 5:6] <-
    adonis2(path[rownames(path) %in% clinic[["TG"]]$sample,] ~ TG,
            clinic[["TG"]][clinic[["TG"]]$sample %in% rownames(path),],
            permutations = 999,
            distance = 'bray')[1, c(3, 5)] %>% as.matrix()
  adonis_result[2, 5:6] <-
    adonis2(path[rownames(path) %in% clinic[["Chol"]]$sample,] ~ Chol,
            clinic[["Chol"]][clinic[["Chol"]]$sample %in% rownames(path),],
            permutations = 999,
            distance = 'bray')[1, c(3, 5)] %>% as.matrix()
  adonis_result[3, 5:6] <-
    adonis2(path[rownames(path) %in% clinic[["HDL"]]$sample,] ~ HDL,
            clinic[["HDL"]][clinic[["HDL"]]$sample %in% rownames(path),],
            permutations = 999,
            distance = 'bray')[1, c(3, 5)] %>% as.matrix()
  adonis_result[4, 5:6] <-
    adonis2(path[rownames(path) %in% clinic[["LDL"]]$sample,] ~ LDL,
            clinic[["LDL"]][clinic[["LDL"]]$sample %in% rownames(path),],
            permutations = 999,
            distance = 'bray')[1, c(3, 5)] %>% as.matrix()
  adonis_result[5, 5:6] <-
    adonis2(path[rownames(path) %in% clinic[["ALB"]]$sample,] ~ ALB,
            clinic[["ALB"]][clinic[["ALB"]]$sample %in% rownames(path),],
            permutations = 999,
            distance = 'bray')[1, c(3, 5)] %>% as.matrix()
  adonis_result[6, 5:6] <-
    adonis2(path[rownames(path) %in% clinic[["Cr"]]$sample,] ~ Cr,
            clinic[["Cr"]][clinic[["Cr"]]$sample %in% rownames(path),],
            permutations = 999,
            distance = 'bray')[1, c(3, 5)] %>% as.matrix()
  adonis_result[7, 5:6] <-
    adonis2(path[rownames(path) %in% clinic[["CTnT"]]$sample,] ~ cTnT,
            clinic[["CTnT"]][clinic[["CTnT"]]$sample %in% rownames(path),],
            permutations = 999,
            distance = 'bray')[1, c(3, 5)] %>% as.matrix()
  adonis_result[8, 5:6] <-
    adonis2(path[rownames(path) %in% clinic[["BNP"]]$sample,] ~ BNP,
            clinic[["BNP"]][clinic[["BNP"]]$sample %in% rownames(path),],
            permutations = 999,
            distance = 'bray')[1, c(3, 5)] %>% as.matrix()
  
  adonis_result[9, 5:6] <-
    adonis2(path[rownames(path) %in% clinic[["HbA1"]]$sample,] ~ HbA1c,
            clinic[["HbA1"]][clinic[["HbA1"]]$sample %in% rownames(path),],
            permutations = 999,
            distance = 'bray')[1, c(3, 5)] %>% as.matrix()
  adonis_result[10, 5:6] <-
    adonis2(path[rownames(path) %in% clinic[["sex"]]$samples,] ~ gender,
            clinic[["sex"]][clinic[["sex"]]$samples %in% rownames(path),],
            permutations = 999,
            distance = 'bray')[1, c(3, 5)] %>% as.matrix()
  adonis_result[11, 5:6] <-
    adonis2(path[rownames(path) %in% clinic[["age"]]$samples,] ~ age,
            clinic[["age"]][clinic[["sex"]]$samples %in% rownames(path),],
            permutations = 999,
            distance = 'bray')[1, c(3, 5)] %>% as.matrix()
}
{
  adonis_result[1, 7:8] <-
    adonis2(ko[rownames(ko) %in% clinic[["TG"]]$sample,] ~ TG,
            clinic[["TG"]][clinic[["TG"]]$sample %in% rownames(ko),],
            permutations = 999,
            distance = 'bray')[1, c(3, 5)] %>% as.data()
  adonis_result[2, 7:8] <-
    adonis2(ko[rownames(ko) %in% clinic[["Chol"]]$sample,] ~ Chol,
            clinic[["Chol"]][clinic[["Chol"]]$sample %in% rownames(ko),],
            permutations = 999,
            distance = 'bray')[1, c(3, 5)] %>% as.matrix()
  adonis_result[3, 7:8] <-
    adonis2(ko[rownames(ko) %in% clinic[["HDL"]]$sample,] ~ HDL,
            clinic[["HDL"]][clinic[["HDL"]]$sample %in% rownames(ko),],
            permutations = 999,
            distance = 'bray')[1, c(3, 5)] %>% as.matrix()
  adonis_result[4, 7:8] <-
    adonis2(ko[rownames(ko) %in% clinic[["LDL"]]$sample,] ~ LDL,
            clinic[["LDL"]][clinic[["LDL"]]$sample %in% rownames(ko),],
            permutations = 999,
            distance = 'bray')[1, c(3, 5)] %>% as.matrix()
  adonis_result[5, 7:8] <-
    adonis2(ko[rownames(ko) %in% clinic[["ALB"]]$sample,] ~ ALB,
            clinic[["ALB"]][clinic[["ALB"]]$sample %in% rownames(ko),],
            permutations = 999,
            distance = 'bray')[1, c(3, 5)] %>% as.matrix()
  adonis_result[6, 7:8] <-
    adonis2(ko[rownames(ko) %in% clinic[["Cr"]]$sample,] ~ Cr,
            clinic[["Cr"]][clinic[["Cr"]]$sample %in% rownames(ko),],
            permutations = 999,
            distance = 'bray')[1, c(3, 5)] %>% as.matrix()
  adonis_result[7, 7:8] <-
    adonis2(ko[rownames(ko) %in% clinic[["CTnT"]]$sample,] ~ cTnT,
            clinic[["CTnT"]][clinic[["CTnT"]]$sample %in% rownames(ko),],
            permutations = 999,
            distance = 'bray')[1, c(3, 5)] %>% as.matrix()
  adonis_result[8, 7:8] <-
    adonis2(ko[rownames(ko) %in% clinic[["BNP"]]$sample,] ~ BNP,
            clinic[["BNP"]][clinic[["BNP"]]$sample %in% rownames(ko),],
            permutations = 999,
            distance = 'bray')[1, c(3, 5)] %>% as.matrix()
  
  adonis_result[9, 7:8] <-
    adonis2(ko[rownames(ko) %in% clinic[["HbA1"]]$sample,] ~ HbA1c,
            clinic[["HbA1"]][clinic[["HbA1"]]$sample %in% rownames(ko),],
            permutations = 999,
            distance = 'bray')[1, c(3, 5)] %>% as.matrix()
  adonis_result[10, 7:8] <-
    adonis2(ko[rownames(ko) %in% clinic[["sex"]]$samples,] ~ gender,
            clinic[["sex"]][clinic[["sex"]]$samples %in% rownames(ko),],
            permutations = 999,
            distance = 'bray')[1, c(3, 5)] %>% as.matrix()
  adonis_result[11, 7:8] <-
    adonis2(ko[rownames(ko) %in% clinic[["age"]]$samples,] ~ age,
            clinic[["age"]][clinic[["sex"]]$samples %in% rownames(ko),],
            permutations = 999,
            distance = 'bray')[1, c(3, 5)] %>% as.matrix()
}
#R2 value
adonis_R2 <- melt(adonis_result[, c(1, 3, 5)])
#create a fake R2 value matrix to adjust color
adonis_R2_f <- melt(adonis_result[, c(1, 3, 5)])
adonis_R2_f$value[adonis_R2$value > 0.04] <-
  adonis_R2$value[adonis_R2$value > 0.04] %>% rescale(to = c(0.04, 0.05)) #让颜色分布均匀
#p value
adonis_p <- melt(adonis_result[, c(2, 4, 6)])
adonis_p$value <-
  ifelse(adonis_p$value <= 0.05,
         ifelse(
           adonis_p$value <= 0.01,
           ifelse(adonis_p$value <= 0.001, "***", "**"),
           "*"
         ),
         NA)
#join R2 and p value tables
adonis_p$Var2 <- str_replace_all(adonis_p$Var2, "_p", "_R2")
adonis_R2 <-
  left_join(adonis_R2,
            adonis_p,
            by = c("Var1", "Var2"),
            suffix = c("_R", "_p"))
ggplot(adonis_R2_f, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  geom_text(
    aes(label = label_percent(accuracy = 0.1)(round(
      adonis_R2$value_R, 3
    ))),
    color =
      ifelse(adonis_R2$value_R > 0.025, "white", "black"),
    size = 4
  ) +
  geom_text(
    aes(label = adonis_R2$value_p),
    color =
      ifelse(adonis_R2$value_R > 0.025, "white", "black"),
    size = 4.5,
    vjust = -0.2
  ) +
  scale_fill_distiller(direction = 1) +
  scale_y_discrete(position = "right") +
  coord_fixed() + theme_minimal() + theme(legend.position = "none")
# adon <- function(x) {
#   temp <- matrix(ncol = 2,dimnames = list(NA,c("R2","p")))
#   # temp_name <- c(TG,Chol,HDL,LDL,ALB,Cr,CTnT,BNP,EF,HbA1c)
#   for (i in 1:length(clinic)) {
# print(i)
#     temp[i,] <-
#       adonis2(x ~ c(TG,Chol,HDL,LDL,ALB,Cr,CTnT,BNP,EF,HbA1c)[i],
#               clinic[[i]],
#               permutations = 999,
#               distance = 'bray')
#   }
#   temp
# }
# sapply(c(tax_s,tax_g), adon)

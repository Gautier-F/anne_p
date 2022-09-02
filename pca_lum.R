library(openxlsx)
library(FactoMineR)
library(factoextra)
library(missMDA)
library(readxl)
library(writexl)
library(httpgd)
hgd()



meta = readRDS("meta_data.RDS")
dat = read.xlsx("data/MCL1-THEREX-J2S63 vs NT-Données brutes normalisées-MAJ 25.06.2022.xlsx", sheet=2)

# l'échantillon NT_8A_2028 est retiré de l'étude car clivage spontané de PARP
dat = dat[-6,]

# char to num
for(i in 2:dim(dat)[2]){
  dat[,i]=as.numeric(dat[,i])
}
str(dat)
dim(dat)
tnt = c(rep("NT", 26), rep("S63", 27))
dat = cbind(tnt, dat)
dat[1:5,1:5]
colnames(dat)[1] = "N_tum"

# rename tumor's name
list_name = strsplit(dat[,1], "-", fixed = T)
name = c()
for (n in seq_along(list_name)){
    name[n] = as.numeric(list_name[[n]][3])
}
name
dat[,1] = name
str(dat[,1])
# changer HER2 en lum
meta$Classement[c(7,17)] = "Lum"
meta$Lum = meta$Classement == "Lum"


saveRDS(meta, "meta_data.rds")
N_tum_Lum = meta$N_tum[meta$Lum]
N_tum_Lum = na.omit(N_tum_Lum)
dat_lum = dat[dat[,2] %in% N_tum_Lum, ]
dim(dat_lum)
dat_lum[dat_lum[,1] =="S63",2] = paste("_", dat_lum[dat_lum[,1] =="S63",2], sep = "")
dat_lum[, 1:2]

# pca sur la totalité des lum

## retirer prot avec na >5

del_na = function(dat){

    l=length(colnames(dat))
    ind = c()
    nbNA= c()
    for ( i in 1:l){
      s = sum(is.na(dat[,i]))
      nbNA = append(nbNA, s)
      if (s >=6){
        ind = append(ind, i)
      }
    }
    print(length(ind)) # 3245 gènes ayant plus de 5 NAs
    dat = dat[,-ind]
    return(dat)
}
sum(is.na(dat_lum))
dat_lum = del_na(dat_lum)
sum(is.na(dat_lum))

# imputation
imp_dat_lum = imputePCA(dat_lum[, 4:dim(dat_lum)[2]], ncp=10)
sum(is.na(imp_dat_lum))

# pca
pca_lum= PCA(imp_dat_lum$completeObs, graph= FALSE)

rownames(pca_lum$ind$coord) = dat_lum[,2]

# Screeplot
scr_plt = fviz_screeplot(pca_lum, addlabels=TRUE) + labs(title="Screeplot")
scr_plt

pdf(file = "screeplot_Lum.pdf")
scr_plt
dev.off()

# données PCA
var = get_pca_var(pca_lum)

cos21 = sort(var$cos2[,1], decreasing = TRUE)
cosdim1 = as.data.frame(names(cos21))
cosdim1$cos2 = cos21


cos22 = sort(var$cos2[,2], decreasing = TRUE)
cosdim2 = as.data.frame(names(cos22))
cosdim2$cos2 = cos22

cos23 = sort(var$cos2[,3], decreasing = TRUE)
cosdim3 = as.data.frame(names(cos22))
cosdim3$cos2 = cos23

coord1 = sort(var$coord[,1], decreasing = TRUE)
coorddim1= as.data.frame(names(coord1))
coorddim1$coord=coord1

coord2 = sort(var$coord[,2], decreasing = TRUE)
coorddim2= as.data.frame(names(coord2))
coorddim2$coord=coord2

coord3 = sort(var$coord[,3], decreasing = TRUE)
coorddim3= as.data.frame(names(coord3))
coorddim3$coord=coord3

write_xlsx(list("Dim1_cos2"=cosdim1, 
                "Dim2_cos2"=cosdim2, "Dim3_cos2"=cosdim3,"Dim1_coord"=coorddim1,
                "Dim2_coord"=coorddim2, "Dim3_coord"=coorddim3), "variables_pca_Lum.xlsx")

hgd()
bi_plt_1_2 = fviz_pca_biplot(pca_lum, 
                          mean.point = FALSE, col.ind = dat_lum[,1], legend.title = "Tnt",
                          pointsize = dat_lum[,3],
                         select.var = list(name = c(names(cos21)[1:10]), 
                                           names(cos22)[1:10]), 
                         repel = TRUE) + labs(title = "Lum")
bi_plt_1_2

bi_plt_1_3 = fviz_pca_biplot(pca_lum, axes = c(1,3),  
                          mean.point = FALSE, col.ind = dat_lum[,1], legend.title = "Tnt",
                         select.var = list(name = c(names(cos21)[1:10]), 
                                           names(cos23)[1:10]), 
                         repel = TRUE) + labs(title = "Lum")

bi_plt_1_3

bi_plt_2_3 = fviz_pca_biplot(pca_lum, axes = c(2,3), 
                           
                        mean.point = FALSE, col.ind = dat_lum[,1], legend.title = "Tnt",
 
                        select.var = list(name = c(names(cos22)[1:10]), 
                                           names(cos23)[1:10]), 
                         repel = TRUE) + labs(title = "Lum")
bi_plt_2_3

pdf("biplots_lum.pdf")
bi_plt_1_2 
bi_plt_1_3
bi_plt_2_3
dev.off()



# #################### PAR CONDITION ###################################################

dat_lum_nt = dat_lum[dat_lum[,1] =="NT",]
dat_lum_nt[1:5,1:5]
# imputation
imp_dat_lum_nt = imputePCA(dat_lum_nt[, 4:dim(dat_lum_nt)[2]], ncp=10)
sum(is.na(imp_dat_lum))

# pca
pca_lum_nt= PCA(imp_dat_lum_nt$completeObs, graph= FALSE)

rownames(pca_lum_nt$ind$coord) = dat_lum_nt[,2]

# Screeplot
scr_plt = fviz_screeplot(pca_lum_nt, addlabels=TRUE) + labs(title="Screeplot")
scr_plt

pdf(file = "screeplot_Lum_nt.pdf")
scr_plt
dev.off()

# données PCA
var = get_pca_var(pca_lum_nt)

cos21 = sort(var$cos2[,1], decreasing = TRUE)
cosdim1 = as.data.frame(names(cos21))
cosdim1$cos2 = cos21


cos22 = sort(var$cos2[,2], decreasing = TRUE)
cosdim2 = as.data.frame(names(cos22))
cosdim2$cos2 = cos22

cos23 = sort(var$cos2[,3], decreasing = TRUE)
cosdim3 = as.data.frame(names(cos22))
cosdim3$cos2 = cos23

coord1 = sort(var$coord[,1], decreasing = TRUE)
coorddim1= as.data.frame(names(coord1))
coorddim1$coord=coord1

coord2 = sort(var$coord[,2], decreasing = TRUE)
coorddim2= as.data.frame(names(coord2))
coorddim2$coord=coord2

coord3 = sort(var$coord[,3], decreasing = TRUE)
coorddim3= as.data.frame(names(coord3))
coorddim3$coord=coord3

write_xlsx(list("Dim1_cos2"=cosdim1, 
                "Dim2_cos2"=cosdim2, "Dim3_cos2"=cosdim3,"Dim1_coord"=coorddim1,
                "Dim2_coord"=coorddim2, "Dim3_coord"=coorddim3), "variables_pca_Lum_nt.xlsx")

hgd()
bi_plt_1_2 = fviz_pca_biplot(pca_lum_nt, 
                          mean.point = FALSE, 
                          pointsize = dat_lum_nt[,3],
                         select.var = list(name = c(names(cos21)[1:10]), 
                                           names(cos22)[1:10]), 
                         repel = TRUE) + labs(title = "Lum_nt", subtitle = "pointsize: %parp clvd")
bi_plt_1_2

bi_plt_1_3 = fviz_pca_biplot(pca_lum_nt, axes = c(1,3),  
                          mean.point = FALSE, pointsize = dat_lum_nt[,3],
                         select.var = list(name = c(names(cos21)[1:10]), 
                                           names(cos23)[1:10]), 
                         repel = TRUE) + labs(title = "Lum_nt", subtitle = "pointsize: %parp clvd")

bi_plt_1_3

bi_plt_2_3 = fviz_pca_biplot(pca_lum_nt, axes = c(2,3), 
                           
                        mean.point = FALSE, pointsize = dat_lum_nt[,3],
 
                        select.var = list(name = c(names(cos22)[1:10]), 
                                           names(cos23)[1:10]), 
                         repel = TRUE) + labs(title = "Lum_nt", subtitle = "pointsize: %parp clvd")
bi_plt_2_3

pdf("biplots_lum_nt.pdf")
bi_plt_1_2 
bi_plt_1_3
bi_plt_2_3
dev.off()




############### S63 ####################################################################



hgd()

dat_lum_s63 = dat_lum[dat_lum[,1] =="S63",]
dat_lum_s63[1:5,1:5]
# imputation
imp_dat_lum_s63 = imputePCA(dat_lum_s63[, 4:dim(dat_lum_s63)[2]], ncp=10)
sum(is.na(imp_dat_lum))

# pca
pca_lum_s63= PCA(imp_dat_lum_s63$completeObs, graph= FALSE)

rownames(pca_lum_s63$ind$coord) = dat_lum_s63[,2]

# Screeplot
scr_plt = fviz_screeplot(pca_lum_s63, addlabels=TRUE) + labs(title="Screeplot")
scr_plt

pdf(file = "screeplot_Lum_s63.pdf")
scr_plt
dev.off()

# données PCA
var = get_pca_var(pca_lum_s63)

cos21 = sort(var$cos2[,1], decreasing = TRUE)
cosdim1 = as.data.frame(names(cos21))
cosdim1$cos2 = cos21


cos22 = sort(var$cos2[,2], decreasing = TRUE)
cosdim2 = as.data.frame(names(cos22))
cosdim2$cos2 = cos22

cos23 = sort(var$cos2[,3], decreasing = TRUE)
cosdim3 = as.data.frame(names(cos22))
cosdim3$cos2 = cos23

coord1 = sort(var$coord[,1], decreasing = TRUE)
coorddim1= as.data.frame(names(coord1))
coorddim1$coord=coord1

coord2 = sort(var$coord[,2], decreasing = TRUE)
coorddim2= as.data.frame(names(coord2))
coorddim2$coord=coord2

coord3 = sort(var$coord[,3], decreasing = TRUE)
coorddim3= as.data.frame(names(coord3))
coorddim3$coord=coord3

write_xlsx(list("Dim1_cos2"=cosdim1, 
                "Dim2_cos2"=cosdim2, "Dim3_cos2"=cosdim3,"Dim1_coord"=coorddim1,
                "Dim2_coord"=coorddim2, "Dim3_coord"=coorddim3), "variables_pca_Lum_s63.xlsx")

hgd()
bi_plt_1_2 = fviz_pca_biplot(pca_lum_s63, 
                          mean.point = FALSE, 
                          pointsize = dat_lum_s63[,3],
                         select.var = list(name = c(names(cos21)[1:10]), 
                                           names(cos22)[1:10]), 
                         repel = TRUE) + labs(title = "Lum_s63", subtitle = "pointsize: %parp clvd")
bi_plt_1_2

bi_plt_1_3 = fviz_pca_biplot(pca_lum_s63 , axes = c(1,3),  
                          mean.point = FALSE, pointsize = dat_lum_s63[,3],
                         select.var = list(name = c(names(cos21)[1:10]), 
                                           names(cos23)[1:10]), 
                         repel = TRUE) + labs(title = "Lum_s63", subtitle = "pointsize: %parp clvd")

bi_plt_1_3

bi_plt_2_3 = fviz_pca_biplot(pca_lum_s63, axes = c(2,3), 
                           
                        mean.point = FALSE, pointsize = dat_lum_s63[,3],
 
                        select.var = list(name = c(names(cos22)[1:10]), 
                                           names(cos23)[1:10]), 
                         repel = TRUE) + labs(title = "Lum_s63", subtitle = "pointsize: %parp clvd")
bi_plt_2_3

pdf("biplots_lum_s63.pdf")
bi_plt_1_2 
bi_plt_1_3
bi_plt_2_3
dev.off()


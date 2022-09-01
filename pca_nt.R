
library(openxlsx)
library(FactoMineR)
library(factoextra)
library(missMDA)
library(readxl)
library(writexl)
library(httpgd)
hgd()

dat = read.xlsx("data/MCL1-THEREX-J2S63 vs NT-Données brutes normalisées-MAJ 25.06.2022.xlsx", sheet=2)

# l'échantillon NT_8A_2028 est retiré de l'étude car clivage spontané de PARP
dat = dat[-6,]

# char to num
for(i in 2:dim(dat)[2]){
  dat[,i]=as.numeric(dat[,i])
}
str(dat)
dim(dat)

rownames(dat)[1] = "N_tum"

# rename tumor's name
list_name = strsplit(dat[,1], "-", fixed = T)
name = c()
for (n in seq_along(list_name)){
    name[n] = as.numeric(list_name[[n]][3])
}
name
dat[,1] = name

dat_nt = dat[1:26,]
dat_s63 = dat[27:53,]
dim(dat_s63)

dat_nt = dat_nt[order(dat_nt[,1]),]
rownames(dat_nt) = 1:26
dat_s63 = dat_s63[order(dat_s63[,1]),]
rownames(dat_s63) = 1:27

meta_nt = dat_nt[,c(1,2)]
meta_s63 = dat_s63[,c(1,2)]

meta = cbind(meta_nt, meta_s63[meta_s63$X1 %in% meta_nt$X1,])
meta =  meta[,-3]
colnames(meta) = c("N_tum", "sens_nt", "sens_s63")
meta$delta_sens = meta$sens_s63 - meta$sens_nt

saveRDS(meta, "meta_clvd_parp.rds")




# recherche des gènes qui ont plus de 5 NA
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

dat_nt_na = del_na(dat_nt)
dat_s63_na = del_na(dat_s63)






# imputation
est_ncp = estim_ncp(dat_s63_na, ncp.max = 5)
imp_dat_nt = imputePCA(dat_nt_na, ncp=10)
sum(is.na(imp_dat))
sum(is.na(dat_nt_na))
dim(imp_dat$completeObs)

saveRDS(imp_dat_nt, "imputed_data_NT.rds")

imp_dat_s63 = readRDS('imputed_data_S63.rds')
dim(imp_dat_s63$completeObs)

#PCA
# enlever la tumeur 2028 de la série S63 qui n'a pas de correspondance avec nt
pca_s63 = PCA(imp_dat_s63$completeObs[-8,], graph= FALSE)

pca_nt = PCA(imp_dat_nt$completeObs, graph= FALSE)

# Screeplot
scr_plt = fviz_screeplot(pca, addlabels=TRUE) + labs(title="Screeplot")
scr_plt

pdf(file = "screeplot_CS.pdf")
scr_plt
dev.off()

# données PCA
var = get_pca_var(pca_nt)

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
                "Dim2_coord"=coorddim2, "Dim3_coord"=coorddim3), "variables_pca_NT.xlsx")

# plot pca
# # dim 1 / 2
meta = readRDS("meta_clvd_parp.rds")
str(meta)
bi_plt_1_2 = fviz_pca_biplot(pca_nt, geom.ind = "point", pointsize = meta$delta_sens,
                         mean.point = FALSE, 
                         select.var = list(name = c(names(cos21)[1:10]), 
                                           names(cos22)[1:10]), 
                         repel = TRUE) + labs(title = "NT", subtitle = "pointsize: delta_parp_clvd")
bi_plt_1_2

bi_plt_1_3 = fviz_pca_biplot(pca_nt, axes = c(1,3), geom.ind = "point", pointsize =meta$delta_sens,
                         mean.point = FALSE, 
                         select.var = list(name = c(names(cos21)[1:10]), 
                                           names(cos23)[1:10]), 
                         repel = TRUE) + labs(title = "NT", subtitle = "pointsize: delta_parp_clvd")

bi_plt_1_3

bi_plt_2_3 = fviz_pca_biplot(pca_nt, axes = c(2,3), geom.ind = "point", pointsize = meta$delta_sens,
                         mean.point = FALSE, 
                         select.var = list(name = c(names(cos22)[1:10]), 
                                           names(cos23)[1:10]), 
                         repel = TRUE) + labs(title = "NT", subtitle = "pointsize: delta_parp_clvd")
bi_plt_2_3

pdf("biplots_NT.pdf")
bi_plt_1_2 
bi_plt_1_3
bi_plt_2_3
dev.off()


# indication du nom des tumeurs
org_name = rownames(imp_dat$completeObs)
org_name
list_name = strsplit(org_name, "-", fixed = T)
name = c()
for (l in seq_along(list_name)){
  name[l]=list_name[[l]][3]
}
name = trimws(name)
name =  paste(c(rep('', 26), rep('_', 27)), name, sep = "")
name
rownames(pca$ind$coord) = name
 p_1_2 = fviz_pca_ind(pca, pointsize = clvd_parp,
                         col.ind = tnt, pointshape = 19, legend.title = "Treatment",
                         mean.point = FALSE, 
                         
                         repel = TRUE) + labs(title = "pca CS")
p_2_3 = fviz_pca_ind(pca, axes = c(2,3),pointsize = clvd_parp,
                         col.ind = tnt, pointshape = 19, legend.title = "Treatment",
                         mean.point = FALSE, 
                         repel = TRUE) + labs(title = "pca CS")
                         
p_1_3 = fviz_pca_ind(pca, axes = c(1,3),pointsize = clvd_parp,
                         col.ind = tnt, pointshape = 19, legend.title = "Treatment",
                         mean.point = FALSE,  
                         repel = TRUE) + labs(title = "pca CS")
pdf("pca_ind_named.pdf")
p_1_2
p_1_3
p_2_3
dev.off()



dim(pca$ind)
p_1_2

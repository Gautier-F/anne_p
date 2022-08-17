library(openxlsx)
library(FactoMineR)
library(factoextra)
library(missMDA)
library(readxl)
library(writexl)
library(httpgd)
hgd()

dat = read.xlsx("data/MCL1-THEREX-J2S63 vs NT-Données brutes normalisées-MAJ 25.06.2022.xlsx", sheet=2)
dim(dat)
head(dat)
colnames(dat)[1:10]

rownames(dat) = dat[,1]
dat = dat[, -1]
rownames(dat)

tnt = c(rep("NT", 27), rep("S63", 26))
clvd_parp = dat[,1]

dat = dat[, -1]

# l'échantillon NT_8A_2028 est retiré de l'étude car clivage spontané de PARP
dat = dat[-6,]
rownames(dat)

# char to num
for(i in 1:dim(dat)[2]){
  dat[,i]=as.numeric(dat[,i])
}
str(dat)

#qq stats
s_na = sum(is.na(dat))
tot = dim(dat)[1]*dim(dat)[2]
pct_na = s_na/tot*100
pct_na # 21% de na

# recherche des gènes qui ont plus de 5 NA
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
length(ind) # 3245 gènes ayant plus de 5 NAs
dat = dat[,-ind]
dim(dat)

# imputation
est_ncp = estim_ncp(dat, ncp.max = 5)
imp_dat = imputePCA(dat, ncp=10)
sum(is.na(imp_dat))
sum(is.na(dat))
dim(imp_dat$completeObs)

#PCA
pca = PCA(imp_dat$completeObs, graph= FALSE)

tnt = as.factor(tnt)
levels(tnt)

# Screeplot
scr_plt = fviz_screeplot(pca, addlabels=TRUE) + labs(title="Screeplot")
scr_plt

pdf(file = "screeplot_CS.pdf")
scr_plt
dev.off()

# données PCA
var = get_pca_var(pca)

cos21 = sort(var$cos2[,1], decreasing = TRUE)
cosdim1 = as.data.frame(names(cos21))
cosdim1$cos2 = cos21


cos22 = sort(var$cos2[,2], decreasing = TRUE)
cosdim2 = as.data.frame(names(cos22))
cosdim2$cos2 = cos22

coord1 = sort(var$coord[,1], decreasing = TRUE)
coorddim1= as.data.frame(names(coord1))
coorddim1$coord=coord1

coord2 = sort(var$coord[,2], decreasing = TRUE)
coorddim2= as.data.frame(names(coord2))
coorddim2$coord=coord2

write_xlsx(list("Dim1_cos2"=cosdim1, 
                "Dim2_cos2"=cosdim2, "Dim1_coord"=coorddim1,
                "Dim2_coord"=coorddim2), "CS_variables_pca.xlsx")

# plot pca
# # dim 1 / 2
bi_plt_1_2 = fviz_pca_biplot(pca, geom.ind = "point", pointsize = 2,
                         col.ind = tnt, pointshape = 19, legend.title = "Treatment",
                         mean.point = FALSE, 
                         select.var = list(name = c(names(cos21)[1:10]), 
                                           names(cos22)[1:10]), 
                         repel = TRUE) + labs(title = "pca CS")
bi_plt_1_2

bi_plt_1_3 = fviz_pca_biplot(pca, axes = c(1,3), geom.ind = "point", pointsize = 2,
                         col.ind = tnt, pointshape = 19, legend.title = "Treatment",
                         mean.point = FALSE, 
                         select.var = list(name = c(names(cos21)[1:10]), 
                                           names(cos22)[1:10]), 
                         repel = TRUE) + labs(title = "pca CS")
bi_plt_1_3

bi_plt_2_3 = fviz_pca_biplot(pca, axes = c(2,3), geom.ind = "point", pointsize = 2,
                         col.ind = tnt, pointshape = 19, legend.title = "Treatment",
                         mean.point = FALSE, 
                         select.var = list(name = c(names(cos21)[1:10]), 
                                           names(cos22)[1:10]), 
                         repel = TRUE) + labs(title = "pca CS")
bi_plt_2_3

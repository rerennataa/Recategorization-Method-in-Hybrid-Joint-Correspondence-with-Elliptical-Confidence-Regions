#========= LIBRARY ==========
library(FactoMineR)
library(stats)
library(datasets)
library(magrittr)
library(dplyr)
library(tidyverse)
#library(writexl)
library(ca)
library("factoextra")
library(ggplot2)


#========== MULTIPLE CORRESPONDENCE ANALYSIS ============
##== Input Data ====
data <- read.csv(file.choose(),header=TRUE,sep=";",dec=",") 


##== Tabel Kontingensi ====
x1=data[,c(2,3)]
x2=data[,c(2,4)]
x3=data[,c(2,5)]
x4=data[,c(2,6)]
x5=data[,c(2,7)]
x6=data[,c(2,8)]
x7=data[,c(2,9)]
x8=data[,c(2,10)]


contingency1 = caconv(x1, from="rpm", to="freq")
cont1 = as.table(contingency1)
names(dimnames(cont1))<-c("Kecamatan","X1")

chisq.test(cont1)

contingency2 = caconv(x2, from="rpm", to="freq")
cont2 = as.table(contingency2)
names(dimnames(cont2))<-c("Kecamatan","X2")

chisq.test(cont2)

contingency3 = caconv(x3, from="rpm", to="freq")
cont3 = as.table(contingency3)
names(dimnames(cont3))<-c("Kecamatan","X3")

chisq.test(cont3)

contingency4 = caconv(x4, from="rpm", to="freq")
cont4 = as.table(contingency4)
names(dimnames(cont4))<-c("Kecamatan","X4")

chisq.test(cont4)

contingency5 = caconv(x5, from="rpm", to="freq")
cont5 = as.table(contingency5)
names(dimnames(cont5))<-c("Kecamatan","X5")

chisq.test(cont5)

contingency6 = caconv(x6, from="rpm", to="freq")
cont6 = as.table(contingency6)
names(dimnames(cont6))<-c("Kecamatan","X6")

chisq.test(cont6)

contingency7 = caconv(x7, from="rpm", to="freq")
cont7 = as.table(contingency7)
names(dimnames(cont7))<-c("Kecamatan","X7")

chisq.test(cont7)

contingency8 = caconv(x8, from="rpm", to="freq")
cont8 = as.table(contingency8)
names(dimnames(cont7))<-c("Kecamatan","X8")

chisq.test(cont8)

write.table(cont1, file="cont1.csv", 
            row.names=T,sep=";")

databaru = data.frame(c(data[2],data[4:9]))
JCA_real <- mjca(databaru, nd=2, lambda = "JCA")
summary(JCA_real) ###inersia total
JCA_real$colinertia

singular=JCA_real$sv
# Hitung total inersia
eigen = singular^2
total_inertia <- sum(eigen)
diagonal = diag(eigen)
# Hitung persentase inersia untuk setiap dimensi
inertia_percentage <- (eigen / total_inertia) * 100
cumulative_inertia_percentage <- cumsum(inertia_percentage)


eig.val <- get_eigenvalue(JCA_real)
eig.val$eigenvalue



###### CONFIDENCE REGIONS #####
regions <- function(coord, xcoord, ycoord, col = col){
  t <- seq(0, 2*pi, length = 1000)
  pcoord1 <- coord[1] + xcoord*cos(t)
  pcoord2 <- coord[2] + ycoord*sin(t)
  lines(pcoord1, pcoord2, col = col)
}

ca.regions.exe <- function (N, a1 = 1, a2 = 2, alpha = 0.1, cols = c(2, 4), M =
                              min(nrow(cont2), ncol(cont2)) - 1, region = 2, scaleplot = 1.2) {
  #############################################################
  # #
  # Defining features of the contingency table for CA #
  # #
  #############################################################
  I <- nrow(N) # Number of rows of table
  J <- ncol(N) # Number of columns of table
  Inames <- dimnames(N)[[1]] # Row category names
  Jnames <- dimnames(N)[[2]] # Column category names
  n <- sum(N) # Total number classified in the table
  p <- N *(1/n) # Matrix of joint relative proportions
  Imass <- as.matrix(apply(p, 1, sum))
  Jmass <- as.matrix(apply(p, 2, sum))
  ItJ <- Imass %*% t(Jmass)
  y <- p - ItJ
  dI <- diag(Imass[1:I])
  dJ <- diag(Jmass[1:J])
  Ih <- Imass^-0.5
  Jh <- Jmass^-0.5
  dIh <- diag(Ih[1:I])
  dJh <- diag(Jh[1:J])
  x <- dIh%*%y%*%dJh
  sva <- svd(x)
  a <- dIh%*%sva$u
  b <- dJh%*%sva$v
  dmu <- diag(sva$d) # Diagonal matrix of singular values
  f <- a %*% dmu # Row coordinates for Classical CA
  g <- b %*% dmu # Column coordinates for Classical CA
  dimnames(f)[[1]] <- Inames
  dimnames(g)[[1]] <- Jnames
  Principal.Inertia <- diag(t(f[, 1:min(I-1, J-1)])%*%dI%*%
                              f[, 1:min(I-1,J-1)])
  Total.Inertia <- sum(Principal.Inertia)
  Percentage.Inertia <- (Principal.Inertia/Total.Inertia) * 100
  Total.Perc.Inertia.M <- sum(Principal.Inertia[1:M])
  chisq.val <- qchisq(1-alpha, df = (I - 1) * (J - 1) )
  #############################################################
  # #
  # Construction of correspondence plot #
  # #
  #############################################################
  par(pty = "s")
  plot(0, 0, pch = " ", xlim = scaleplot * range(f[, 1:M], g[, 1:M]),
       ylim = scaleplot * range(f[, 1:M], g[, 1:M]),
       xlab = paste("Principal Axis ", a1, "(", round(Percentage.Inertia[a1],
                                                      digits = 2), "%)"), ylab = paste("Principal Axis ", a2, "(",
                                                                                       round(Percentage.Inertia[a2], digits = 2), "%)"))
  text(f[,1], f[,2], labels = Inames, adj = 0, col = cols[1])
  points(f[, a1], f[, a2], pch = "*", col = cols[1])
  text(g[,1], g[,2], labels = Jnames, adj = 1, col = cols[2])
  points(g[, a1], g[, a2], pch = "#", col = cols[2])
  abline(h = 0, v = 0)
  title(main = paste(100 * (1 - alpha), "% Confidence Regions"))
  #############################################################
  # #
  # Calculating the row and column radii length for a #
  # confidence circle #
  # #
  #############################################################
  radii <- sqrt(qchisq(1 - alpha, 2)/(n * Imass))
  radij <- sqrt(qchisq(1 - alpha, 2)/(n * Jmass))
  #############################################################
  # #
  # Calculating the semi-axis lengths for the confidence #
  # ellipses #
  # #
  #############################################################
  hlax1.row <- vector(mode = "numeric", length = I)
  hlax2.row <- vector(mode = "numeric", length = I)
  hlax1.col <- vector(mode = "numeric", length = J)
  hlax2.col <- vector(mode = "numeric", length = J)
  if (M > 2){
    # Semi-axis lengths for the row coordinates in an optimal plot
    for (i in 1:I){
      hlax1.row[i] <- dmu[1,1] * sqrt((chisq.val/(n*Total.Inertia))*
                                        (1/Imass[i] - sum(a[i, 3:M]^2)))
      hlax2.row[i] <- dmu[2,2] * sqrt((chisq.val/(n*Total.Inertia))*
                                        (1/Imass[i] - sum(a[i, 3:M]^2)))
    }
    # Semi-axis lengths for the column coordinates in an optimal plot
    for (j in 1:J){
      hlax1.col[j] <- dmu[1,1] * sqrt((chisq.val/(n * Total.Inertia))*
                                        (1/Jmass[j] - sum(b[j, 3:M]^2)))
      hlax2.col[j] <- dmu[2,2] * sqrt((chisq.val/(n * Total.Inertia))*
                                        (1/Jmass[j] - sum(b[j, 3:M]^2)))
    }
  } else {
    # Semi-axis lengths for the row coordinates in a two-dimensional plot
    for (i in 1:I){
      hlax1.row[i] <- dmu[1,1] * sqrt((chisq.val/(n * Total.Inertia))*
                                        (1/Imass[i]))
      hlax2.row[i] <- dmu[2,2] * sqrt((chisq.val/(n * Total.Inertia))*
                                        (1/Imass[i]))
    }
    # Semi-axis lengths for the column coordinates in a two-dimensional plot
    for (j in 1:J){
      hlax1.col[j] <- dmu[1,1] * sqrt((chisq.val/(n * Total.Inertia))*
                                        (1/Jmass[j]))
      hlax2.col[j] <- dmu[2,2] * sqrt((chisq.val/(n * Total.Inertia))*
                                        (1/Jmass[j]))
    }
  }
  #############################################################
  # #
  # Eccentricity #
  # #
  #############################################################
  eccentricity <- sqrt(1-(dmu[a2,a2]/dmu[a1,a1])^2)
  #############################################################
  # #
  # Approximate P-values #
  # #
  #############################################################
  pvalrow <- vector(mode = "numeric", length = I)
  pvalrowcircle <- vector(mode = "numeric", length = I)
  pvalcol <- vector(mode = "numeric", length = J)
  pvalcolcircle <- vector(mode = "numeric", length = J)
  for (i in 1:I){
    # Approximate row P-values from Lebart et al.’s (1984) confidence
    # circles
    pvalrowcircle[i] <- 1 - pchisq(n * Imass[i] * (f[i, 1]^2 + f[i, 2]^2),
                                   df = (I-1)*(J-1))
    # Approximate P-values based on Beh’s (2010) confidence ellipses
    if (M > 2){
      pvalrow[i] <- 1 - pchisq(n * Total.Inertia * ((1/Imass[i] -
                                                       sum(a[i, 3:M]^2))^(-1)) * ((f[i, 1]/dmu[1, 1])^2
                                                                                  + (f[i, 2]/dmu[2, 2])^2), df = (I - 1) * (J - 1))
    } else {
      pvalrow[i] <- 1 - pchisq(n * Total.Inertia * Imass[i]*((f[i, 1]/
                                                                dmu[1, 1])^2 + (f[i, 2]/dmu[2, 2])^2),
                               df = (I - 1) * (J - 1))
    }
  }
  for (j in 1:J){
    # Approximate row P-values based on Lebart et al.’s (1984)
    # confidence circles
    pvalcolcircle[j] <- 1 - pchisq(n * Imass[i] * (g[j, 1]^2 + g[j, 2]^2),
                                   df = (I - 1) * (J - 1))
    # Approximate P-values based on Beh’s (2010) confidence ellipses
    if (M > 2){
      pvalcol[j] <- 1 - pchisq(n * Total.Inertia * ((1/Jmass[j] -
                                                       sum(b[j, 3:M]^2))^(-1)) * ((g[j, 1]/dmu[1, 1])^2
                                                                                  + (g[j, 2]/dmu[2, 2])^2), df = (I - 1)*(J - 1))
    } else {
      pvalcol[j] <- 1 - pchisq(n * Total.Inertia * Jmass[j]*((g[j,1]/
                                                                dmu[1, 1])^2 + (g[j, 2]/dmu[2, 2])^2),
                               df = (I - 1) * (J - 1))
    }
  }
  summ.name <- c("HL Axis 1", "HL Axis 2", "P-value-ellipse",
                 "P-value-circle")
  if (region == 1){
    row.summ <- cbind(radii, radii, pvalrow, pvalrowcircle)
    col.summ <- cbind(radij, radij, pvalcol, pvalcolcircle)
  } else if (region == 2){
    row.summ <- cbind(hlax1.row, hlax2.row, pvalrow, pvalrowcircle)
    col.summ <- cbind(hlax1.col, hlax2.col, pvalcol, pvalcolcircle)
  }
  dimnames(row.summ) <- list(paste(Inames), paste(summ.name))
  dimnames(col.summ) <- list(paste(Jnames), paste(summ.name))
  #############################################################
  # #
  # Superimposing the confidence regions #
  # #
  #############################################################
  if (region == 1){
    # Superimposing the confidence circles
    symbols(f[,a1], f[,a2], circles = radii, add = T, fg = cols[1])
    symbols(g[,a1], g[,a2], circles = radij, add = T, fg = cols[2])
  } else if (region == 2){
    # Superimposing the confidence ellipses
    for (i in 1:I){
      regions(f[i,], xcoord = hlax1.row[i], ycoord = hlax2.row[i],
              col = cols[1])
    }
    for (j in 1:J){
      regions(g[j,], xcoord = hlax1.col[j], ycoord = hlax2.col[j],
              col = cols[2])
    }
  }
  #############################################################
  # #
  # Summary of output #
  # #
  #############################################################
  if (region == 1){
    list(Row.Summary = round(row.summ, digits = 4), Column.Summary =
           round(col.summ, digits = 4),
         Inertia = Principal.Inertia)
  } else if (region == 2){
    list(Eccentricity = round(eccentricity, digits = 4), Row.Summary
         = round(row.summ, digits = 4), Column.Summary =
           round(col.summ, digits = 4), P=round(p, digits = 4),
         rm=round(Imass, digits=4),cm=round(Jmass, digits=4), 
         S= round(x, digits=4), dI, dJ,sva, PCB=round(f,digits=4),
         PCK=round(g,digits=4))
  }
}


##### CORRESPONDENCE ANALYSIS ######
corsp.analysis<-function(matrix)
{
  corsp <- ca(matrix)
  KU <- cacoord(corsp,type = c("principal"), dim = NA)
  KUC <- KU$columns
  KUR <- KU$rows
  distance <- dist(KUC,method="euclidean")
  jarak = as.matrix(distance)
  hasil=list(Koordinat.Utama = KUC, Jarak = jarak)
  print(hasil)
}

write.table(cors46, file="cors46.csv", 
            row.names=T,sep=";")

##### VARIABEL TEMPAT PEMBUANGAN AKHIR TINJA #####
cr2 = ca.regions.exe((data.matrix(cont2)), region = 2); cr2
cors2 = corsp.analysis(cont2)
#==== GABUNG DATA ====
x22 = x2
x22$X2[x22$X2 == "Lubang tanah"] <- "LT"
x22$X2[x22$X2 == "Tangki/ instalasi pengelolaan air limbah"] <- "LT"
#==== TABEL KONTINGENSI ====
contingency22 = caconv(x22, from="rpm", to="freq")
cont22 = as.table(contingency22)
names(dimnames(cont22))<-c("Kecamatan","X2")
chisq.test(cont22)
#==== CONFIDENCE REGIONS ====
cr22 = ca.regions.exe((data.matrix(cont22)), region = 2); cr22


##### VARIABEL TEMPAT PEMBUANGAN LIMBAH CAIR #####
cr3 = ca.regions.exe((data.matrix(cont3)), region = 2); cr3
cors3 = corsp.analysis(cont3)
#==== GABUNG DATA ====
x32 = x3
x32$X3[x32$X3 == "Drainase(Got/ Selokan)"] <- "DS"
x32$X3[x32$X3 == "Sungai/ Saluran Irigasi"] <- "DS"
#==== TABEL KONTINGENSI ====
contingency32 = caconv(x32, from="rpm", to="freq")
cont32 = as.table(contingency32)
names(dimnames(cont32))<-c("Kecamatan","X3")
chisq.test(cont32)
#==== CONFIDENCE REGIONS ====
cr32 = ca.regions.exe((data.matrix(cont32)), region = 2); cr32
cors32 = corsp.analysis(cont32)
#==== GABUNG DATA ====
x33 = x32
x33$X3[x33$X3 == "DS"] <- "DSD"
x33$X3[x33$X3 == "Dalam lubang/ Tanah terbuka"] <- "DSD"
#==== TABEL KONTINGENSI ====
contingency33 = caconv(x33, from="rpm", to="freq")
cont33 = as.table(contingency33)
names(dimnames(cont33))<-c("Kecamatan","X3")
chisq.test(cont33)
#==== CONFIDENCE REGIONS ====
cr33 = ca.regions.exe((data.matrix(cont33)), region = 2); cr33
cors33 = corsp.analysis(cont33)

#==== GABUNG DATA ====
#x34 = x33
#x34$X3[x34$X3 == "Dalam lubang/ Tanah terbuka"] <- "DSLD"
#x34$X3[x34$X3 == "DSL"] <- "DSLD"
#==== TABEL KONTINGENSI ====
#contingency34 = caconv(x34, from="rpm", to="freq")
#cont34 = as.table(contingency34)
#names(dimnames(cont34))<-c("Kecamatan","X3")
#chisq.test(cont34)
#==== CONFIDENCE REGIONS ====
#cr34 = ca.regions.exe((data.matrix(cont34)), region = 2); cr34


##### VARIABEL SUMBER AIR MINUM #####
cr4 = ca.regions.exe((data.matrix(cont4)), region = 2); cr4
cors4 = corsp.analysis(cont4)
#==== GABUNG DATA ====
x42 = x4
x42$X4[x42$X4 == "Air isi ulang"] <- "AS"
x42$X4[x42$X4 == "Sumur bor atau pompa"] <- "AS"
#==== TABEL KONTINGENSI ====
contingency42 = caconv(x42, from="rpm", to="freq")
cont42 = as.table(contingency42)
names(dimnames(cont42))<-c("Kecamatan","X4")
chisq.test(cont42)
#==== CONFIDENCE REGIONS ====
cr42 = ca.regions.exe((data.matrix(cont42)), region = 2); cr42
cors42 = corsp.analysis(cont42)
#==== GABUNG DATA ====
x43 = x42
x43$X4[x43$X4 == "AS"] <- "ASS"
x43$X4[x43$X4 == "Sumur"] <- "ASS"
#==== TABEL KONTINGENSI ====
contingency43 = caconv(x43, from="rpm", to="freq")
cont43 = as.table(contingency43)
names(dimnames(cont43))<-c("Kecamatan","X4")
chisq.test(cont43)
#==== CONFIDENCE REGIONS ====
cr43 = ca.regions.exe((data.matrix(cont43)), region = 2); cr43
cors43 = corsp.analysis(cont43)
#==== GABUNG DATA ====
x44 = x43
x44$X4[x44$X4 == "ASS"] <- "ASSM"
x44$X4[x44$X4 == "Mata air"] <- "ASSM"
#==== TABEL KONTINGENSI ====
contingency44 = caconv(x44, from="rpm", to="freq")
cont44 = as.table(contingency44)
names(dimnames(cont44))<-c("Kecamatan","X4")
chisq.test(cont44)
#==== CONFIDENCE REGIONS ====
cr44 = ca.regions.exe((data.matrix(cont44)), region = 2); cr44
cors44 = corsp.analysis(cont44)
#==== GABUNG DATA ====
x45 = x44
x45$X4[x45$X4 == "ASSM"] <- "ASSML"
x45$X4[x45$X4 == "Ledeng dengan meteran"] <- "ASSML"
#==== TABEL KONTINGENSI ====
contingency45 = caconv(x45, from="rpm", to="freq")
cont45 = as.table(contingency45)
names(dimnames(cont45))<-c("Kecamatan","X4")
chisq.test(cont45)
#==== CONFIDENCE REGIONS ====
cr45 = ca.regions.exe((data.matrix(cont45)), region = 2); cr45
cors45 = corsp.analysis(cont45)
#==== GABUNG DATA ====
x46 = x45
x46$X4[x45$X4 == "ASSML"] <- "ASSMLL"
x46$X4[x45$X4 == "Ledeng tanpa meteran"] <- "ASSMLL"
#==== TABEL KONTINGENSI ====
contingency46 = caconv(x46, from="rpm", to="freq")
cont46 = as.table(contingency46)
names(dimnames(cont46))<-c("Kecamatan","X4")
chisq.test(cont46)
#==== CONFIDENCE REGIONS ====
cr46 = ca.regions.exe((data.matrix(cont46)), region = 2); cr46


##### VARIABEL SUMBER AIR MANDI #####
cr5 = ca.regions.exe((data.matrix(cont5)), region = 2); cr5
cors5 = corsp.analysis(cont5)
#==== GABUNG DATA ====
x52 = x5
x52$X5[x52$X5 == "Ledeng dengan meteran"] <- "LS"
x52$X5[x52$X5 == "Sumur bor atau pompa"] <- "LS"
#==== TABEL KONTINGENSI ====
contingency52 = caconv(x52, from="rpm", to="freq")
cont52 = as.table(contingency52)
names(dimnames(cont52))<-c("Kecamatan","X5")
chisq.test(cont52)
#==== CONFIDENCE REGIONS ====
cr52 = ca.regions.exe((data.matrix(cont52)), region = 2); cr52
cors52 = corsp.analysis(cont52)
#==== GABUNG DATA ====
x53 = x52
x53$X5[x53$X5 == "Sungai/ danau/ kolam/ waduk/ situ/ embung/ bendungan"] <- "SS"
x53$X5[x53$X5 == "Sumur"] <- "SS"
#==== TABEL KONTINGENSI ====
contingency53 = caconv(x53, from="rpm", to="freq")
cont53 = as.table(contingency53)
names(dimnames(cont53))<-c("Kecamatan","X5")
chisq.test(cont53)
#==== CONFIDENCE REGIONS ====
cr53 = ca.regions.exe((data.matrix(cont53)), region = 2); cr53
cors53 = corsp.analysis(cont53)
#==== GABUNG DATA ====
x54 = x53
x54$X5[x54$X5 == "LS"] <- "LSSS"
x54$X5[x54$X5 == "SS"] <- "LSSS"
#==== TABEL KONTINGENSI ====
contingency54 = caconv(x54, from="rpm", to="freq")
cont54 = as.table(contingency54)
names(dimnames(cont54))<-c("Kecamatan","X5")
chisq.test(cont54)
#==== CONFIDENCE REGIONS ====
cr54 = ca.regions.exe((data.matrix(cont54)), region = 2); cr54


##### VARIABEL FUNGSI KAWASAN HUTAN #####
cr7 = ca.regions.exe((data.matrix(cont7)), region = 2); cr7
cors7 = corsp.analysis(cont7)
#==== GABUNG DATA ====
x72 = x7
x72$X7[x72$X7 == "Produksi"] <- "PD"
x72$X7[x72$X7 == "Diluar Kawasan Hutan"] <- "PD"
#==== TABEL KONTINGENSI ====
contingency72 = caconv(x72, from="rpm", to="freq")
cont72 = as.table(contingency72)
names(dimnames(cont72))<-c("Kecamatan","X7")
chisq.test(cont72)
#==== CONFIDENCE REGIONS ====
cr72 = ca.regions.exe((data.matrix(cont72)), region = 2); cr72
cors72 = corsp.analysis(cont72)
#==== GABUNG DATA ====
x73 = x72
x73$X7[x73$X7 == "Konservasi"] <- "PDK"
x73$X7[x73$X7 == "PD"] <- "PDK"
#==== TABEL KONTINGENSI ====
contingency73 = caconv(x73, from="rpm", to="freq")
cont73 = as.table(contingency73)
names(dimnames(cont73))<-c("Kecamatan","X7")
chisq.test(cont73)
#==== CONFIDENCE REGIONS ====
cr73 = ca.regions.exe((data.matrix(cont73)), region = 2); cr73
cors73 = corsp.analysis(cont73)
#==== GABUNG DATA ====
x74 = x73
x74$X7[x74$X7 == "Bukan fungsi kawasan hutan"] <- "PDKB"
x74$X7[x74$X7 == "PDK"] <- "PDKB"
#==== TABEL KONTINGENSI ====
contingency74 = caconv(x74, from="rpm", to="freq")
cont74 = as.table(contingency74)
names(dimnames(cont74))<-c("Kecamatan","X7")
chisq.test(cont74)
#==== CONFIDENCE REGIONS ====
cr74 = ca.regions.exe((data.matrix(cont74)), region = 2); cr74


##### MULTIPLE CORRESPONDENCE ANALYSIS #####
df = as.data.frame(c(x22,(x33[2]),(x46[2]),(x54[2]),(x6[2]),(x74[2])))
head(df)

##== JCA ====
#Matriks asosiasi MCA
dataJCA <- mjca(df, nd=2, lambda = "JCA",maxit=100,reti=TRUE,epsilon=0.0001)
datamca=mjca(df, nd=2, lambda = "Burt")
summary(datamca)
summary(dataJCA) ###inersia total
plot(dataJCA)
View(dataJCA$Burt.upd)
View(dataJCA$indmat)
ind=dataJCA$indmat
write.table(ind, file="Matriks Indikator Awal.csv", 
            row.names=T,sep=";")
k=dim(ind)[2]
B=dataJCA27$Burt
write.table(B, file="Matriks Burt Akhir.csv", 
            row.names=T,sep=";")
B2=dataJCA27$Burt.upd
write.table(B2, file="Matriks Rekonstruksi Burt Akhir.csv", 
            row.names=T,sep=";")
b=sum(B)
P=B/b
write.table(P, file="Matriks Korespondensi Awal.csv", 
            row.names=T,sep=";")
View(S)
c=apply(P,2,sum)
write.table(c, file="Vektor Proporsi Kolom Awal.csv", 
            row.names=T,sep=";")
Dc=diag(c)
eP <- c %*% t(c)
diagC=diag(c)
sn <- solve(sqrt(diagC)) %*% (P - eP) %*% solve(sqrt(diagC))
S <- (P - eP) / sqrt(eP)
write.table(S, file="Matriks Residual Akhir.csv", 
            row.names=T,sep=";")
View(sn)

#eigen_f = (datamca$sv)^2
diagv = diag(eigen_f)
View(diage)
dec=eigen(S)
L=k-7
eig_value=dec$values[1:L]
write.table(eig_value, file="Vektor Eigen Awal.csv", 
            row.names=T,sep=";")
diage=diag(eig_value)
V=dec$vectors[,1:L]
write.table(V, file="Matriks V Awal.csv", 
            row.names=T,sep=";")
SR = V%*%diage%*%(t(V))

#__Koordinat Standar__#
KS_r<-cacoord(dataJCA,type = c("standard"))
KSC_r<-KS_r$columns
write.table(KSC_r, file="Matriks Koordinat Standar Awal.csv", 
            row.names=T,sep=";")
View(KSC_r)
plot(dataJCA)

Kor = dataJCA$colpcoord[1:31,]
View(Kor)
write.table(Kor, file="Matriks Koordinat Utama Awal.csv", 
            row.names=T,sep=";")
dataJCA$colcoord[1:31,]
distanceK <- dist(Kor,method="euclidean")
jarakK = as.matrix(distanceK)
write.table(jarakK, file="Euclidean pertama.csv", 
            row.names=T,sep=";")
I = diag(31)
jarakK_ = jarakK + I*100
rownames(jarakK_) = unique(df$Kecamatan)
colnames(jarakK_) = unique(df$Kecamatan)
min(jarakK_)
View(jarakK_)
s = sort(jarakK_)
quantile(s[1:930], probs=0.47)

#==== GABUNG DATA ====
df2 = df


df2$Kecamatan[df2$Kecamatan == "Ciparay"] <- "Ciparay dan Majalaya"
df2$Kecamatan[df2$Kecamatan == "Majalaya"]<- "Ciparay dan Majalaya"
dataJCA2 <- mjca(df2, lambda = "JCA",reti=TRUE, epsilon=0.0001)
summary(dataJCA2) ###inersia total


Kor2 = dataJCA2$colpcoord[1:30,]
distanceK2 <- dist(Kor2,method="euclidean")
jarakK2 = as.matrix(distanceK2)
I = diag(30)
jarakK2_ = jarakK2 + I*100
rownames(jarakK2_) = unique(df2$Kecamatan)
colnames(jarakK2_) = unique(df2$Kecamatan)
sort(jarakK2_)
View(jarakK2_)

df3 = df2

df3$Kecamatan[df3$Kecamatan == "Rancaekek"] <- "Rancaekek dan Solokan Jeruk"
df3$Kecamatan[df3$Kecamatan == "Solokanjeruk"]<- "Rancaekek dan Solokan Jeruk"
dataJCA3 <- mjca(df3, lambda = "JCA")
summary(dataJCA3) ###inersia total
plot(dataJCA3)

Kor3 = dataJCA3$colpcoord[1:29,]
distanceK3 <- dist(Kor3,method="euclidean")
jarakK3 = as.matrix(distanceK3)
I = diag(29)
jarakK3_ = jarakK3 + I*100
rownames(jarakK3_) = unique(df3$Kecamatan)
colnames(jarakK3_) = unique(df3$Kecamatan)
sort(jarakK3_)
View(jarakK3_)


df4 = df3

df4$Kecamatan[df4$Kecamatan == "Cikancung"] <- "Cikancung dan Nagreg"
df4$Kecamatan[df4$Kecamatan == "Nagreg"]<- "Cikancung dan Nagreg"
dataJCA4 <- mjca(df4, lambda = "JCA")
summary(dataJCA4) ###inersia total
plot(dataJCA4)

Kor4 = dataJCA4$colpcoord[1:28,]
distanceK4 <- dist(Kor4,method="euclidean")
jarakK4 = as.matrix(distanceK4)
I = diag(28)
jarakK4_ = jarakK4 + I*100
rownames(jarakK4_) = unique(df4$Kecamatan)
colnames(jarakK4_) = unique(df4$Kecamatan)
sort(jarakK4_)
View(jarakK4_)


df5 = df4
df5$Kecamatan[df5$Kecamatan == "Cikancung dan Nagreg"] <- "Cikancung, Nagreg, Cicalengka"
df5$Kecamatan[df5$Kecamatan == "Cicalengka"]<- "Cikancung, Nagreg, Cicalengka"
dataJCA5 <- mjca(df5, lambda = "JCA")
summary(dataJCA5) ###inersia total
plot(dataJCA5)

Kor5 = dataJCA5$colpcoord[1:27,]
distanceK5 <- dist(Kor5,method="euclidean")
jarakK5 = as.matrix(distanceK5)
I = diag(27)
jarakK5_ = jarakK5 + I*100
rownames(jarakK5_) = unique(df5$Kecamatan)
colnames(jarakK5_) = unique(df5$Kecamatan)
sort(jarakK5_)
View(jarakK5_)


df6 = df5
df6$Kecamatan[df6$Kecamatan == "Bojongsoang"] <- "Bojongsoang dan Dayeuhkolot"
df6$Kecamatan[df6$Kecamatan == "Dayeuhkolot"]<- "Bojongsoang dan Dayeuhkolot"
dataJCA6 <- mjca(df6, nd = 2, lambda = "JCA")
summary(dataJCA6) ###inersia total
plot(dataJCA6)

Kor6 = dataJCA6$colpcoord[1:26,]
distanceK6 <- dist(Kor6,method="euclidean")
jarakK6 = as.matrix(distanceK6)
I = diag(26)
jarakK6_ = jarakK6 + I*100
rownames(jarakK6_) = unique(df6$Kecamatan)
colnames(jarakK6_) = unique(df6$Kecamatan)
sort(jarakK6_)
View(jarakK6_)


df7 = df6
df7$Kecamatan[df7$Kecamatan == "Ciparay dan Majalaya"] <- "Ciparay, Majalaya, Rancaekek, dan Solokan Jeruk"
df7$Kecamatan[df7$Kecamatan == "Rancaekek dan Solokan Jeruk"]<- "Ciparay, Majalaya, Rancaekek, dan Solokan Jeruk"
dataJCA7 <- mjca(df7, nd = 2, lambda = "JCA")
summary(dataJCA7) ###inersia total
plot(dataJCA7)

Kor7 = dataJCA7$colpcoord[1:25,]
distanceK7 <- dist(Kor7,method="euclidean")
jarakK7 = as.matrix(distanceK7)
I = diag(25)
jarakK7_ = jarakK7 + I*100
rownames(jarakK7_) = unique(df7$Kecamatan)
colnames(jarakK7_) = unique(df7$Kecamatan)
sort(jarakK7_)
View(jarakK7_)


df8 = df7
df8$Kecamatan[df8$Kecamatan == "Ciparay, Majalaya, Rancaekek, dan Solokan Jeruk"] <- "Ciparay, Majalaya, Rancaekek, Solokan Jeruk, dan Cileunyi"
df8$Kecamatan[df8$Kecamatan == "Cileunyi"]<- "Ciparay, Majalaya, Rancaekek, Solokan Jeruk, dan Cileunyi"
dataJCA8 <- mjca(df8, nd = 2, lambda = "JCA")
summary(dataJCA8) ###inersia total
plot(dataJCA8)

Kor8 = dataJCA8$colpcoord[1:24,]
distanceK8 <- dist(Kor8,method="euclidean")
jarakK8 = as.matrix(distanceK8)
I = diag(24)
jarakK8_ = jarakK8 + I*100
rownames(jarakK8_) = unique(df8$Kecamatan)
colnames(jarakK8_) = unique(df8$Kecamatan)
sort(jarakK8_)
View(jarakK8_)



df9 = df8  ##0.345
df9$Kecamatan[df9$Kecamatan == "Kutawaringin"] <- "Kutawaringin dan Soreang"
df9$Kecamatan[df9$Kecamatan == "Soreang"] <- "Kutawaringin dan Soreang"
#df9$Kecamatan[df9$Kecamatan == "Ciparay, Majalaya, Rancaekek, Solokan Jeruk, dan Cileunyi"] <- "Ciparay, Majalaya, Rancaekek, Solokan Jeruk, Cileunyi, dan Pacet"
#df9$Kecamatan[df9$Kecamatan == "Pacet"] <- "Ciparay, Majalaya, Rancaekek, Solokan Jeruk, Cileunyi, dan Pacet"
dataJCA9 <- mjca(df9, nd = 2, lambda = "JCA")
summary(dataJCA9) ###inersia total
plot(dataJCA9)

Kor9 = dataJCA9$colpcoord[1:23,]
distanceK9 <- dist(Kor9,method="euclidean")
jarakK9 = as.matrix(distanceK9)
I = diag(23)
jarakK9_ = jarakK9 + I*100
rownames(jarakK9_) = unique(df9$Kecamatan)
colnames(jarakK9_) = unique(df9$Kecamatan)
sort(jarakK9_)
View(jarakK9_)



df10 = df9
df10$Kecamatan[df10$Kecamatan == "Margaasih"] <- "Kutawaringin, Soreang, dan Margaasih"
df10$Kecamatan[df10$Kecamatan == "Kutawaringin dan Soreang"] <- "Kutawaringin, Soreang, dan Margaasih"
dataJCA10 <- mjca(df10, nd = 2, lambda = "JCA")
summary(dataJCA10) ###inersia total
plot(dataJCA10)

Kor10 = dataJCA10$colpcoord[1:22,]
distanceK10 <- dist(Kor10,method="euclidean")
jarakK10 = as.matrix(distanceK10)
I = diag(22)
jarakK10_ = jarakK10 + I*100
rownames(jarakK10_) = unique(df10$Kecamatan)
colnames(jarakK10_) = unique(df10$Kecamatan)
sort(jarakK10_)
View(jarakK10_)



df11 = df10  #0.357
df11$Kecamatan[df11$Kecamatan == "Ibun"] <- "Ibun dan Pacet"
df11$Kecamatan[df11$Kecamatan == "Pacet"]<- "Ibun dan Pacet"
dataJCA11 <- mjca(df11, nd = 2, lambda = "JCA")
summary(dataJCA11) ###inersia total
plot(dataJCA11)

Kor11 = dataJCA11$colpcoord[1:21,]
distanceK11 <- dist(Kor11,method="euclidean")
jarakK11 = as.matrix(distanceK11)
I = diag(21)
jarakK11_ = jarakK11 + I*100
rownames(jarakK11_) = unique(df11$Kecamatan)
colnames(jarakK11_) = unique(df11$Kecamatan)
sort(jarakK11_)
View(jarakK11_)


df12 = df11
df12$Kecamatan[df12$Kecamatan == "Katapang"] <- "Katapang, Kutawaringin, Soreang, dan Margaasih"
df12$Kecamatan[df12$Kecamatan == "Kutawaringin, Soreang, dan Margaasih"]<- "Katapang, Kutawaringin, Soreang, dan Margaasih"
dataJCA12 <- mjca(df12, nd = 2, lambda = "JCA")
summary(dataJCA12) ###inersia total
plot(dataJCA12)

Kor12 = dataJCA12$colpcoord[1:20,]
distanceK12 <- dist(Kor12,method="euclidean")
jarakK12 = as.matrix(distanceK12)
I = diag(20)
jarakK12_ = jarakK12 + I*100
rownames(jarakK12_) = unique(df12$Kecamatan)
colnames(jarakK12_) = unique(df12$Kecamatan)
sort(jarakK12_)
View(jarakK12_)



df13 = df12
df13$Kecamatan[df13$Kecamatan == "Arjasari"] <- "Arjasari dan Banjaran"
df13$Kecamatan[df13$Kecamatan == "Banjaran"]<- "Arjasari dan Banjaran"
dataJCA13 <- mjca(df13, nd = 2, lambda = "JCA")
summary(dataJCA13) ###inersia total
plot(dataJCA13)

Kor13 = dataJCA13$colpcoord[1:19,]
distanceK13 <- dist(Kor13,method="euclidean")
jarakK13 = as.matrix(distanceK13)
I = diag(19)
jarakK13_ = jarakK13 + I*100
rownames(jarakK13_) = unique(df13$Kecamatan)
colnames(jarakK13_) = unique(df13$Kecamatan)
sort(jarakK13_)
View(jarakK13_)


df14 = df13
df14$Kecamatan[df14$Kecamatan == "Ibun dan Pacet"] <- "Ibun, Pacet, dan Kertasari"
df14$Kecamatan[df14$Kecamatan == "Kertasari"]<- "Ibun, Pacet, dan Kertasari"
dataJCA14 <- mjca(df14, nd = 2, lambda = "JCA")
summary(dataJCA14) ###inersia total
plot(dataJCA14)

Kor14 = dataJCA14$colpcoord[1:18,]
distanceK14 <- dist(Kor14,method="euclidean")
jarakK14 = as.matrix(distanceK14)
I = diag(18)
jarakK14_ = jarakK14 + I*100
rownames(jarakK14_) = unique(df14$Kecamatan)
colnames(jarakK14_) = unique(df14$Kecamatan)
sort(jarakK14_)
View(jarakK14_)


df15 = df14  ##0.454
df15$Kecamatan[df15$Kecamatan == "Cikancung, Nagreg, Cicalengka"] <- "Cikancung, Nagreg, Cicalengka, dan Paseh"
df15$Kecamatan[df15$Kecamatan == "Paseh"]<- "Cikancung, Nagreg, Cicalengka, dan Paseh"
dataJCA15 <- mjca(df15, nd = 2, lambda = "JCA")
summary(dataJCA15) ###inersia total
plot(dataJCA15)

Kor15 = dataJCA15$colpcoord[1:17,]
distanceK15 <- dist(Kor15,method="euclidean")
jarakK15 = as.matrix(distanceK15)
I = diag(17)
jarakK15_ = jarakK15 + I*100
rownames(jarakK15_) = unique(df15$Kecamatan)
colnames(jarakK15_) = unique(df15$Kecamatan)
sort(jarakK15_)
View(jarakK15_)



df16 = df15  ##0.46
df16$Kecamatan[df16$Kecamatan == "Pangalengan"] <- "Ibun, Pacet, Kertasari, dan Pangalengan"
df16$Kecamatan[df16$Kecamatan == "Ibun, Pacet, dan Kertasari"]<- "Ibun, Pacet, Kertasari, dan Pangalengan"
dataJCA16 <- mjca(df16, nd = 2, lambda = "JCA")
summary(dataJCA16) ###inersia total
plot(dataJCA16)

Kor16 = dataJCA16$colpcoord[1:16,]
distanceK16 <- dist(Kor16,method="euclidean")
jarakK16 = as.matrix(distanceK16)
I = diag(16)
jarakK16_ = jarakK16 + I*100
rownames(jarakK16_) = unique(df16$Kecamatan)
colnames(jarakK16_) = unique(df16$Kecamatan)
sort(jarakK16_)
View(jarakK16_)



df17 = df16 #0.494
df17$Kecamatan[df17$Kecamatan == "Ciwidey"] <- "Katapang, Kutawaringin, Soreang, Margaasih, dan Ciwidey"
df17$Kecamatan[df17$Kecamatan == "Katapang, Kutawaringin, Soreang, dan Margaasih"]<- "Katapang, Kutawaringin, Soreang, Margaasih, dan Ciwidey"
dataJCA17 <- mjca(df17, nd = 2, lambda = "JCA")
summary(dataJCA17) ###inersia total
plot(dataJCA17)


Kor17 = dataJCA17$colpcoord[1:15,]
distanceK17 <- dist(Kor17,method="euclidean")
jarakK17 = as.matrix(distanceK17)
I = diag(15)
jarakK17_ = jarakK17 + I*100
rownames(jarakK17_) = unique(df17$Kecamatan)
colnames(jarakK17_) = unique(df17$Kecamatan)
sort(jarakK17_)
View(jarakK17_)


df18 = df17
df18$Kecamatan[df18$Kecamatan == "Arjasari dan Banjaran"] <- "Arjasari, Banjaran, Ibun, Pacet, Kertasari, dan Pangalengan"
df18$Kecamatan[df18$Kecamatan == "Ibun, Pacet, Kertasari, dan Pangalengan"]<- "Arjasari, Banjaran, Ibun, Pacet, Kertasari, dan Pangalengan"
dataJCA18 <- mjca(df18, nd = 2, lambda = "JCA")
summary(dataJCA18) ###inersia total
plot(dataJCA18)

Kor18 = dataJCA18$colpcoord[1:14,]
distanceK18 <- dist(Kor18,method="euclidean")
jarakK18 = as.matrix(distanceK18)
I = diag(14)
jarakK18_ = jarakK18 + I*100
rownames(jarakK18_) = unique(df18$Kecamatan)
colnames(jarakK18_) = unique(df18$Kecamatan)
sort(jarakK18_)
View(jarakK18_)


df19 = df18  #0.45
df19$Kecamatan[df19$Kecamatan == "Arjasari, Banjaran, Ibun, Pacet, Kertasari, dan Pangalengan"] <- "Arjasari, Banjaran, Ibun, Pacet, Kertasari, Pangalengan, dan Cimaung"
df19$Kecamatan[df19$Kecamatan == "Cimaung"]<- "Arjasari, Banjaran, Ibun, Pacet, Kertasari, Pangalengan, dan Cimaung"
dataJCA19 <- mjca(df19, nd = 2, lambda = "JCA")
summary(dataJCA19) ###inersia total
plot(dataJCA19)

Kor19 = dataJCA19$colpcoord[1:13,]
distanceK19 <- dist(Kor19,method="euclidean")
jarakK19 = as.matrix(distanceK19)
I = diag(13)
jarakK19_ = jarakK19 + I*100
rownames(jarakK19_) = unique(df19$Kecamatan)
colnames(jarakK19_) = unique(df19$Kecamatan)
sort(jarakK19_)
View(jarakK19_)



df20 = df19  ##0.46
df20$Kecamatan[df20$Kecamatan == "Katapang, Kutawaringin, Soreang, Margaasih, dan Ciwidey"] <- "Katapang, Kutawaringin, Soreang, Margaasih, Ciwidey, dan Baleendah"
df20$Kecamatan[df20$Kecamatan == "Baleendah"]<- "Katapang, Kutawaringin, Soreang, Margaasih, Ciwidey, dan Baleendah"
dataJCA20 <- mjca(df20, nd = 2, lambda = "JCA",maxit=100)
summary(dataJCA20) ###inersia total
plot(dataJCA20)

Kor20 = dataJCA20$colpcoord[1:12,]
distanceK20 <- dist(Kor20,method="euclidean")
jarakK20 = as.matrix(distanceK20)
I = diag(12)
jarakK20_ = jarakK20 + I*100
rownames(jarakK20_) = unique(df20$Kecamatan)
colnames(jarakK20_) = unique(df20$Kecamatan)
sort(jarakK20_)
View(jarakK20_)


df21 = df20
df21$Kecamatan[df21$Kecamatan == "Katapang, Kutawaringin, Soreang, Margaasih, Ciwidey, dan Baleendah"] <- "Katapang, Kutawaringin, Soreang, Margaasih, Ciwidey, Baleendah, dan Pameungpeuk"
df21$Kecamatan[df21$Kecamatan == "Pameungpeuk"]<- "Katapang, Kutawaringin, Soreang, Margaasih, Ciwidey, Baleendah, dan Pameungpeuk"
dataJCA21 <- mjca(df21, nd = 2, lambda = "JCA",maxit=100)
summary(dataJCA21) ###inersia total
plot(dataJCA21)

Kor21 = dataJCA21$colpcoord[1:11,]
distanceK21 <- dist(Kor21,method="euclidean")
jarakK21 = as.matrix(distanceK21)
I = diag(11)
jarakK21_ = jarakK21 + I*100
rownames(jarakK21_) = unique(df21$Kecamatan)
colnames(jarakK21_) = unique(df21$Kecamatan)
sort(jarakK21_)
View(jarakK21_)


df22 = df21
df22$Kecamatan[df22$Kecamatan == "Arjasari, Banjaran, Ibun, Pacet, Kertasari, Pangalengan, dan Cimaung"] <- "Arjasari, Banjaran, Ibun, Pacet, Kertasari, Pangalengan, Cimaung,dan Cangkuang"
df22$Kecamatan[df22$Kecamatan == "Cangkuang"]<- "Arjasari, Banjaran, Ibun, Pacet, Kertasari, Pangalengan, Cimaung,dan Cangkuang"
dataJCA22 <- mjca(df22, nd = 2, lambda = "JCA",maxit=100)
summary(dataJCA22) ###inersia total
plot(dataJCA22)

Kor22 = dataJCA22$colpcoord[1:10,]
distanceK22 <- dist(Kor22,method="euclidean")
jarakK22 = as.matrix(distanceK22)
I = diag(10)
jarakK22_ = jarakK22 + I*100
rownames(jarakK22_) = unique(df22$Kecamatan)
colnames(jarakK22_) = unique(df22$Kecamatan)
sort(jarakK22_)
View(jarakK22_)


df23 = df22 #0.57
df23$Kecamatan[df23$Kecamatan == "Arjasari, Banjaran, Ibun, Pacet, Kertasari, Pangalengan, Cimaung,dan Cangkuang"] <- "Arjasari, Banjaran, Ibun, Pacet, Kertasari, Pangalengan, Cimaung, Cangkuang, dan Pasirjambu"
df23$Kecamatan[df23$Kecamatan == "Pasirjambu"]<- "Arjasari, Banjaran, Ibun, Pacet, Kertasari, Pangalengan, Cimaung, Cangkuang, dan Pasirjambu"
dataJCA23 <- mjca(df23, nd = 2, lambda = "JCA",maxit=100)
summary(dataJCA23) ###inersia total
plot(dataJCA23)

Kor23 = dataJCA23$colpcoord[1:9,]
distanceK23 <- dist(Kor23,method="euclidean")
jarakK23 = as.matrix(distanceK23)
I = diag(9)
jarakK23_ = jarakK23 + I*100
rownames(jarakK23_) = unique(df23$Kecamatan)
colnames(jarakK23_) = unique(df23$Kecamatan)
sort(jarakK23_)
View(jarakK23_)


df24 = df23 ##0.586
df24$Kecamatan[df24$Kecamatan == "Katapang, Kutawaringin, Soreang, Margaasih, Ciwidey, Baleendah, dan Pameungpeuk"] <- "Katapang, Kutawaringin, Soreang, Margaasih, Ciwidey, Baleendah, Pameungpeuk, Bojongsoang, dan Dayeuhkolot"
df24$Kecamatan[df24$Kecamatan == "Bojongsoang dan Dayeuhkolot"] <- "Katapang, Kutawaringin, Soreang, Margaasih, Ciwidey, Baleendah, Pameungpeuk, Bojongsoang, dan Dayeuhkolot"
dataJCA24 <- mjca(df24, nd = 2, lambda = "JCA",maxit=100)
summary(dataJCA24) ###inersia total
plot(dataJCA24)

Kor24 = dataJCA24$colpcoord[1:8,]
distanceK24 <- dist(Kor24,method="euclidean")
jarakK24 = as.matrix(distanceK24)
I = diag(8)
jarakK24_ = jarakK24 + I*100
rownames(jarakK24_) = unique(df24$Kecamatan)
colnames(jarakK24_) = unique(df24$Kecamatan)
sort(jarakK24_)
View(jarakK24_)


df25 = df24  ##0.6960222
df25$Kecamatan[df25$Kecamatan == "Cilengkrang"] <- "Cilengkrang dan Cimenyan"
df25$Kecamatan[df25$Kecamatan == "Cimenyan"]<- "Cilengkrang dan Cimenyan"
dataJCA25 <- mjca(df25, nd = 2, lambda = "JCA",maxit=100)
summary(dataJCA25) ###inersia total
plot(dataJCA25)

Kor25 = dataJCA25$colpcoord[1:7,]
distanceK25 <- dist(Kor25,method="euclidean")
jarakK25 = as.matrix(distanceK25)
I = diag(7)
jarakK25_ = jarakK25 + I*100
rownames(jarakK25_) = unique(df25$Kecamatan)
colnames(jarakK25_) = unique(df25$Kecamatan)
sort(jarakK25_)
View(jarakK25_)


df26 = df25 
df26$Kecamatan[df26$Kecamatan == "Cilengkrang dan Cimenyan"] <- "Ciparay, Majalaya, Rancaekek, Solokan Jeruk, Cileunyi, Cilengkrang, dan Cimenyan"
df26$Kecamatan[df26$Kecamatan == "Ciparay, Majalaya, Rancaekek, Solokan Jeruk, dan Cileunyi"]<- "Ciparay, Majalaya, Rancaekek, Solokan Jeruk, Cileunyi, Cilengkrang, dan Cimenyan"
dataJCA26 <- mjca(df26, nd = 2, lambda = "JCA", maxit=100)
summary(dataJCA26) ###inersia total
plot(dataJCA26)

Kor26 = dataJCA26$colpcoord[1:6,]
distanceK26 <- dist(Kor26,method="euclidean")
jarakK26 = as.matrix(distanceK26)
I = diag(6)
jarakK26_ = jarakK26 + I*100
rownames(jarakK26_) = unique(df26$Kecamatan)
colnames(jarakK26_) = unique(df26$Kecamatan)
sort(jarakK26_)
View(jarakK26_)


df27 = df26
df27$Kecamatan[df27$Kecamatan == "Cikancung, Nagreg, Cicalengka, dan Paseh"] <- "Cikancung, Nagreg, Cicalengka, Paseh, Ciparay, Majalaya, Rancaekek, Solokan Jeruk, Cileunyi, Cilengkrang, dan Cimenyan"
df27$Kecamatan[df27$Kecamatan == "Ciparay, Majalaya, Rancaekek, Solokan Jeruk, Cileunyi, Cilengkrang, dan Cimenyan"]<- "Cikancung, Nagreg, Cicalengka, Paseh, Ciparay, Majalaya, Rancaekek, Solokan Jeruk, Cileunyi, Cilengkrang, dan Cimenyan"
dataJCA27 <- mjca(df27,graph=TRUE, nd = 2, lambda = "JCA",maxit=100,epsilon=0.0001)
summary(dataJCA27) ###inersia total
plot(dataJCA27)

df_final = df27
df_final$Kecamatan[df_final$Kecamatan == "Cikancung, Nagreg, Cicalengka, Paseh, Ciparay, Majalaya, Rancaekek, Solokan Jeruk, Cileunyi, Cilengkrang, dan Cimenyan"]<- "3"
df_final$Kecamatan[df_final$Kecamatan == "Katapang, Kutawaringin, Soreang, Margaasih, Ciwidey, Baleendah, Pameungpeuk, Bojongsoang, dan Dayeuhkolot"]<- "2"
df_final$Kecamatan[df_final$Kecamatan == "Arjasari, Banjaran, Ibun, Pacet, Kertasari, Pangalengan, Cimaung, Cangkuang, dan Pasirjambu"]<- "1"
df_final$Kecamatan[df_final$Kecamatan == "Margahayu"]<- "4"
df_final$Kecamatan[df_final$Kecamatan == "Rancabali"]<- "5"


colnames(df_final)[colnames(df_final) == "Kecamatan"] <- "K"
#df10$Kecamatan[df10$Kecamatan == "Margaasih"] <- "Kutawaringin, Soreang, dan Margaasih"
df_final$X2 = as.numeric(df_final$X2)
df_final$X3 = as.numeric(df_final$X3)
df_final$X4 = as.numeric(df_final$X4)
df_final$X5 = as.numeric(df_final$X5)
df_final$X6 = as.numeric(df_final$X6)
df_final$X7 = as.numeric(df_final$X7)

df_final$X2 = as.factor(df_final$X2)
df_final$X3 = as.factor(df_final$X3)
df_final$X4 = as.factor(df_final$X4)
df_final$X5 = as.factor(df_final$X5)
df_final$X6 = as.factor(df_final$X6)
df_final$X7 = as.factor(df_final$X7)

df_final$X2[df_final$X2 == "2"]<- "2,4"
df_final$X3[df_final$X3 == "1"]<- "1,2,5"
df_final$X4[df_final$X4 == "1"]<- "1,3,4,5,6,7"
df_final$X5[df_final$X5 == "1"]<- "1,5,6"
df_final$X7[df_final$X7 == "1"]<- "1,2"


str(df_final)

JCA_final <- mjca(df_final, nd =2, lambda = "JCA",maxit=100, reti=TRUE,epsilon=0.0001)
summary(JCA_final)
sv_f = JCA_final$sv
#__Koordinat Standar__#
KS_f<-cacoord(JCA_final,type = c("standard"))
KSC_f<-KS_f$columns
head(KSC_f)
#__Koordinat Utama__#
KU_f<-cacoord(JCA_final,type = c("principal"), dim = NA)
KUC_f<-KU_f$columns
rownames(KUC_f)<-rownames(JCA_final$Burt)
write.table(KUC_f, file="Matriks Koordinat Utama Akhir.csv", 
            row.names=T,sep=";")

####====Plot dari JCA=====####
Keterangan=c(1,1,1,1,1,
             2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
KUC_new=cbind(KUC_f,Keterangan)
ggplot(as.data.frame(KUC_new), aes(x = Dim1, y = Dim2, label = rownames(KUC_new),shape=factor(Keterangan))) +
  geom_point(aes(colour = factor(Keterangan)), size = 4)+
  geom_text(hjust = 0, vjust = 1, size = 6.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)) +
  theme(axis.text.x = element_text(vjust = 1, hjust = 1, size = 12)) +
  theme(axis.text.y = element_text(vjust = 1, hjust = 1, size = 12))

#KUC ZOOM
KUC1=KUC_f[1:3,]
KUC2=KUC_f[6:9,]
KUC3=KUC_f[c(11, 13), ]
KUC4=KUC_f[15:22,]

KUC_p = rbind(KUC1,KUC2,KUC3,KUC4)
Ket=c(1,1,1,2,2,
      2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
KUC_p = cbind(KUC_p,Ket)

ggplot(as.data.frame(KUC_p), aes(x = Dim1, y = Dim2, label = rownames(KUC_p),shape=factor(Ket))) +
  geom_point(aes(colour = factor(Ket)), size = 4)+
  geom_text(hjust = 1, vjust = 1, size = 6.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)) +
  theme(axis.text.x = element_text(vjust = 1, hjust = 1, size = 12)) +
  theme(axis.text.y = element_text(vjust = 1, hjust = 1, size = 12))+
  labs(shape = "Keterangan", colour = "Keterangan")

#MATRIKS JARAK EUCLIDEAN untuk melihat kebergantungan
#K_f<-KUC_f2[1:31,1:2]
distance_f=dist(KUC_f,method="euclidean")
df_distance=as.matrix(distance_f)
View(df_distance)


####====Data Kuantitatif=====####
data_mt2 = read.csv(file.choose(),header=TRUE,sep=";",dec=".") 

#Standarisasi
data_mt2=scale(data_mt2[2:4])
data_mt2=as.data.frame(data_mt2)
str(data_mt2)

#Korelasi Kosinus antara KUC_Kec dan data metrik
KUC_K=cacoord(JCA_final, type=c("principal"), dim=NA)
summary(JCA_final)
KUC_C=KUC_K$columns
rownames(KUC_C)<-rownames(JCA_final$Burt)
KUC_Kec <- KUC_C[1:5, 1:2]
korelasi = function(KUC_Kec, data_mt2){
  if (nrow(KUC_Kec)==nrow(data_mt2)) {
    col = ncol(KUC_Kec)
    fz = matrix(0,ncol(KUC_Kec),ncol(data_mt2))
    for (i in 1:col) {
      fz[i,] = colSums(KUC_Kec[,i]*data_mt2)
    }
    flsqr = colSums(KUC_Kec^2)
    zsqr = colSums(data_mt2^2)
    rho = matrix(0,ncol(KUC_Kec),ncol(data_mt2))
    for (i in 1:ncol(KUC_Kec)) {
      for (j in 1:ncol(data_mt2)) {
        rho[i,j] = fz[i,j]/(sqrt(flsqr[i])*sqrt(zsqr[j]))
      }
    }
    return(rho)
  }
  if (nrow(KUC_Kec)!=nrow(data_mt22)) {
    print()
  }
}
r=korelasi(KUC_Kec, data_mt2)
r
cos1=acos(r)
cos1*(180/pi)

#Koordinat vektor kuantitatif (psi)
sv_f = JCA_final$sv
sv_f = sv_f[1:2]
ev = diag((sv_f)^2) ###eigen value
psi = (ev^(0.08))%*%r ###delta=1/2


#ev = diag((JCA_final$sv)^2) ###eigen value
#psi = (ev^(1/2))%*%d ###alpha=1/2

psi=t(psi)
rownames(psi)=c("Z1","Z2","Z3")
colnames(psi)=c("Dim1", "Dim2")

###====Plot Hybrid Joint Correspondence===####
df_cos=rbind(KUC_f[,1:2],psi)
Ket=c(1,1,1,1,1,
      2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3)
df_cos = cbind(df_cos,Ket)
ggplot(as.data.frame(df_cos), aes(x = Dim1, y = Dim2, label = rownames(df_cos),shape=factor(Ket))) +
  geom_point(aes(colour = factor(Ket)), size = 5)+
  geom_text(hjust = 0, vjust = 1, size = 6.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)) +
  theme(axis.text.x = element_text(vjust = 1, hjust = 1, size = 12)) +
  theme(axis.text.y = element_text(vjust = 1, hjust = 1, size = 12))+
  labs(shape = "Keterangan", colour = "Keterangan")

#KUC ZOOM
KUC_f1=df_cos[1:3,]
KUC_f2=df_cos[6:9,]
KUC_f3=df_cos[c(11, 13), ]
KUC_f4=df_cos[15:25,]

KUC_final = rbind(KUC_f1,KUC_f2,KUC_f3,KUC_f4)
ggplot(as.data.frame(KUC_final), aes(x = Dim1, y = Dim2, label = rownames(KUC_final),shape=factor(Ket))) +
  geom_point(aes(colour = factor(Ket)), size = 4)+
  geom_text(hjust = 1, vjust = 1, size = 6.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)) +
  theme(axis.text.x = element_text(vjust = 1, hjust = 1, size = 12)) +
  theme(axis.text.y = element_text(vjust = 1, hjust = 1, size = 12))+
  labs(shape = "Keterangan", colour = "Keterangan")


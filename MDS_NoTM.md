---
title: "Multidimensional scaling plot"
Author: "Akanksha Pandey"
output: 
  html_document:
    keep_md: true
---

R script to create Multidimesional Scaling plot for different amino acid excahnge matices

#Step1: Load data 
Load your data into R. Here, I have a 38 X 38 matrix of euclidean distances between different amino acid exchange rate matrices



```r
# Classical MDS
# N rows (objects) x p columns (variables)
# each row identified by a unique row name

row_lable <- c('BUR','BUR_choanozoa', 'COIL' ,'COIL_BUR' ,'COIL_BUR_choanozoa' ,'COIL_EXP' ,'COIL_EXP_choanozoa' ,'COIL_choanozoa' ,'DAYHOFF-dcmut' ,'EXP' ,'EXP_choanozoa' ,'FLU' ,'HELIX' ,'HELIX_BUR' ,'HELIX_BUR_choanozoa' ,'HELIX_EXP' ,'HELIX_EXP_choanozoa' ,'HELIX_choanozoa' ,'HIVb' ,'HIVw' ,'JONES-dcmut' ,'JTTtm' ,'LG' ,'SHEET' ,'SHEET_BUR' ,'SHEET_BUR_choanozoa' ,'SHEET_EXP' ,'SHEET_EXP_choanozoa' ,'SHEET_choanozoa' ,'VT' ,'WAG' ,'blosum', 'cpREV' ,'mtART' ,'mtMAM' ,'mtREV24', 'mtZOA' ,'rtREV')

mydata <- read.csv("final_aa_mat.csv",header = TRUE,sep =  "\t", as.is = TRUE ,row.names = row_lable )#read the csv file for disimmilarity matrix 
head(mydata)
```

```
##                            BUR BUR_choanozoa       COIL   COIL_BUR
## BUR                0.000000000   0.006189167 0.09856246 0.04554241
## BUR_choanozoa      0.006189167   0.000000000 0.09937737 0.04557887
## COIL               0.098562459   0.099377366 0.00000000 0.08688022
## COIL_BUR           0.045542409   0.045578867 0.08688022 0.00000000
## COIL_BUR_choanozoa 0.047482323   0.046693555 0.09306460 0.01018766
## COIL_EXP           0.142843160   0.143927234 0.05410296 0.13046229
##                    COIL_BUR_choanozoa   COIL_EXP COIL_EXP_choanozoa
## BUR                        0.04748232 0.14284316         0.15240218
## BUR_choanozoa              0.04669355 0.14392723         0.15308167
## COIL                       0.09306460 0.05410296         0.06146644
## COIL_BUR                   0.01018766 0.13046229         0.13848450
## COIL_BUR_choanozoa         0.00000000 0.13615854         0.14360225
## COIL_EXP                   0.13615854 0.00000000         0.01817214
##                    COIL_choanozoa DAYHOFF.dcmut        EXP EXP_choanozoa
## BUR                    0.09960932    0.07152740 0.15244475    0.15872867
## BUR_choanozoa          0.09998060    0.07280753 0.15355321    0.15951763
## COIL                   0.01010200    0.10306360 0.07012344    0.07283461
## COIL_BUR               0.08726632    0.08986567 0.14210726    0.14691091
## COIL_BUR_choanozoa     0.09279669    0.09471535 0.14764472    0.15207665
## COIL_EXP               0.05736957    0.13969339 0.02751393    0.03015792
##                           FLU      HELIX  HELIX_BUR HELIX_BUR_choanozoa
## BUR                0.09663043 0.06222952 0.03071986          0.03183886
## BUR_choanozoa      0.09727629 0.06349953 0.03158564          0.03151670
## COIL               0.10273474 0.07266973 0.10622259          0.10615163
## COIL_BUR           0.09290712 0.08360689 0.07134083          0.07178313
## COIL_BUR_choanozoa 0.09526943 0.08885910 0.07364748          0.07374841
## COIL_EXP           0.13797503 0.10997865 0.14829879          0.14844504
##                     HELIX_EXP HELIX_EXP_choanozoa HELIX_choanozoa
## BUR                0.17660827          0.17612052      0.06035392
## BUR_choanozoa      0.17786536          0.17722961      0.06093515
## COIL               0.11254625          0.10875276      0.07418683
## COIL_BUR           0.17062163          0.16926655      0.08287343
## COIL_BUR_choanozoa 0.17550010          0.17406268      0.08764002
## COIL_EXP           0.08095629          0.07751003      0.11284123
##                          HIVb      HIVw JONES.dcmut     JTTtm         LG
## BUR                0.09633334 0.1204678  0.06285436 0.2091813 0.05519211
## BUR_choanozoa      0.09680167 0.1215920  0.06415401 0.2088819 0.05738113
## COIL               0.11725777 0.1485661  0.08983750 0.1480237 0.06286402
## COIL_BUR           0.11184162 0.1296017  0.08105831 0.1952305 0.06818686
## COIL_BUR_choanozoa 0.11400637 0.1313965  0.08554805 0.1981620 0.07475314
## COIL_EXP           0.14442399 0.1793648  0.12766232 0.1494845 0.10259669
##                         SHEET  SHEET_BUR SHEET_BUR_choanozoa  SHEET_EXP
## BUR                0.04185122 0.03713236          0.03761826 0.13837545
## BUR_choanozoa      0.04467498 0.03928677          0.03884415 0.14036170
## COIL               0.10100616 0.11829726          0.11887701 0.08218506
## COIL_BUR           0.07240119 0.06598757          0.06553296 0.13476025
## COIL_BUR_choanozoa 0.07710934 0.06797938          0.06704000 0.14123880
## COIL_EXP           0.14244199 0.16055104          0.16150485 0.07034198
##                    SHEET_EXP_choanozoa SHEET_choanozoa         VT
## BUR                         0.13716998      0.04036736 0.08624685
## BUR_choanozoa               0.13894833      0.04257762 0.08678508
## COIL                        0.07609479      0.10104017 0.09614316
## COIL_BUR                    0.13236056      0.07141378 0.09222971
## COIL_BUR_choanozoa          0.13890893      0.07584017 0.09591571
## COIL_EXP                    0.06580529      0.14323513 0.13023794
##                           WAG     blosum      cpREV      mtART     mtMAM
## BUR                0.05542871 0.07256707 0.05682412 0.10379861 0.2188076
## BUR_choanozoa      0.05779849 0.07529712 0.05762473 0.10370017 0.2197905
## COIL               0.08373717 0.08576344 0.10242055 0.06712125 0.2280640
## COIL_BUR           0.06869065 0.07886278 0.07094196 0.10371357 0.2227538
## COIL_BUR_choanozoa 0.07479012 0.08519925 0.07373529 0.10919833 0.2257955
## COIL_EXP           0.12113133 0.11759381 0.13810498 0.09536096 0.2420698
##                       mtREV24      mtZOA     rtREV
## BUR                0.08755890 0.08722499 0.1687189
## BUR_choanozoa      0.08768762 0.08728175 0.1708134
## COIL               0.09409796 0.07014140 0.1846437
## COIL_BUR           0.10009196 0.09258686 0.1762293
## COIL_BUR_choanozoa 0.10485032 0.09784613 0.1799498
## COIL_EXP           0.12737815 0.10482439 0.1980667
```
#Step2: Fit cmdscale function
Here I used cmdscal scale function for dimension reduction in 2 coordinates


```r
fit <- cmdscale(mydata, k=2, eig = FALSE)
fit.frame <- data.frame(fit)
head(fit.frame)
```

```
##                             X1           X2
## BUR                -0.04957964 -0.023269778
## BUR_choanozoa      -0.05004301 -0.024910302
## COIL                0.04165137 -0.012738117
## COIL_BUR           -0.03043953 -0.024002428
## COIL_BUR_choanozoa -0.03405126 -0.027141174
## COIL_EXP            0.08686531  0.002525567
```
###Step3: Plot coordinates from cmdscale function using plot function


```r
# plot solution 
coordinate1 <- fit.frame[,1]
coordinate2 <- fit.frame[,2]

rainbow1<-rainbow(38)
plot(coordinate1, coordinate2, xlab="Coordinate 1", ylab="Coordinate 2", main="MDS for Amino acid exchange rates",	type="p", pch=19,col=rainbow1, cex.main=0.7, las=1, bty="n" )
text(coordinate1, coordinate2, pos=3,labels = row_lable, cex=0.4,col = "Black")
```

![](MDS_NoTM_files/figure-html/unnamed-chunk-3-1.png)<!-- -->


#Step4: Plot coordinates from cmdscale using ggplot
Here I used ggplot for better visulatization. I also colored the datapoints based on groups of exchange rate matrices


```r
library(ggplot2)
library(tidyverse)
```

```
## ── Attaching packages ────────────────────────────────── tidyverse 1.2.1 ──
```

```
## ✔ tibble  1.4.2     ✔ purrr   0.2.4
## ✔ tidyr   0.8.0     ✔ dplyr   0.7.4
## ✔ readr   1.1.1     ✔ stringr 1.2.0
## ✔ tibble  1.4.2     ✔ forcats 0.3.0
```

```
## ── Conflicts ───────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::filter() masks stats::filter()
## ✖ dplyr::lag()    masks stats::lag()
```

```r
library(RColorBrewer)
fit.frame$bio <- c(row_lable)
head(mutate(fit.frame, grepl("EXP", bio)))
```

```
##            X1           X2                bio grepl("EXP", bio)
## 1 -0.04957964 -0.023269778                BUR             FALSE
## 2 -0.05004301 -0.024910302      BUR_choanozoa             FALSE
## 3  0.04165137 -0.012738117               COIL             FALSE
## 4 -0.03043953 -0.024002428           COIL_BUR             FALSE
## 5 -0.03405126 -0.027141174 COIL_BUR_choanozoa             FALSE
## 6  0.08686531  0.002525567           COIL_EXP              TRUE
```

```r
fit.frame$col_bio <- ifelse(grepl("EXP",fit.frame$bio), 1, ifelse(grepl("BUR", fit.frame$bio), 2, ifelse(grepl("HELIX", fit.frame$bio), 3, ifelse(grepl("SHEET", fit.frame$bio), 3, ifelse(grepl("COIL", fit.frame$bio), 3, 4)))))
fit.frame$col_bio <- as.factor(fit.frame$col_bio)
plot1 <- ggplot(fit.frame, aes(coordinate1,coordinate2,label = c(row_lable), color = col_bio)) + geom_point() +geom_text(check_overlap = TRUE, size=2, vjust =1.4, hjust = 1) +  scale_colour_brewer("Amino acid matrix categories", palette="Dark2",type = "diverging", direction = 1, labels = c('Exposed','Buried','Secondary Structure','Standard Models')) + theme_minimal()
print(plot1)
```

![](MDS_NoTM_files/figure-html/unnamed-chunk-4-1.png)<!-- -->


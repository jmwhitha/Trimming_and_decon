#Compute effect size of statistically significant differences
#Reference: Trimming and decontamination of metagenomic data can significantly impact assembly and binning metrics, phylogenomic and functional analysis 
#Jason M. Whitham 2021

#Data for raw (assembly metrics)
raw_reads <-c('9117.5.MEGAHIT.assembly', '10158.8.MEGAHIT.assembly', '11263.1.MEGAHIT.assembly', '11306.3.MEGAHIT.assembly', '11306.1.MEGAHIT.assembly', '11260.6.MEGAHIT.assembly', '11260.5.MEGAHIT.assembly', '9053.2.MEGAHIT.assembly', '9672.8.MEGAHIT.assembly', '9053.4.MEGAHIT.assembly', '9053.3.MEGAHIT.assembly', '9053.5.MEGAHIT.assembly', '9108.2.MEGAHIT.assembly', '10186.3.MEGAHIT.assembly', '7333.1.MEGAHIT.assembly', '9117.8.MEGAHIT.assembly', '10158.6.MEGAHIT.assembly', '9117.7.MEGAHIT.assembly', '9117.6.MEGAHIT.assembly', '10186.4.MEGAHIT.assembly', '9108.1.MEGAHIT.assembly', '9117.4.MEGAHIT.assembly', '9041.8.MEGAHIT.assembly')
raw_total_contigs <-c(1183038.0, 362044.0, 1331533.0, 1042261.0, 506569.0, 1061986.0, 1015988.0, 1215308.0, 1096293.0, 1058433.0, 1184448.0, 721833.0, 1209055.0, 1044708.0, 763500.0, 795124.0, 1021022.0, 1134727.0, 1111363.0, 1138489.0, 518629.0, 1067385.0, 718988.0)
raw_10k_contigs <-c(7559.0, 2026.0, 10913.0, 7077.0, 2490.0, 9668.0, 9498.0, 7610.0, 5822.0, 6798.0, 8794.0, 3215.0, 5807.0, 6714.0, 2571.0, 4479.0, 6864.0, 9220.0, 5534.0, 9160.0, 1654.0, 7524.0, 2652.0)
raw_100k_contigs <-c(46.0, 12.0, 139.0, 77.0, 48.0, 119.0, 139.0, 55.0, 40.0, 46.0, 52.0, 11.0, 13.0, 114.0, 8.0, 21.0, 36.0, 70.0, 31.0, 79.0, 15.0, 41.0, 16.0)
raw_largest_contig <-c(238480.0, 214425.0, 981732.0, 1235329.0, 450787.0, 449860.0, 1831404.0, 303994.0, 229075.0, 225822.0, 213821.0, 223153.0, 162369.0, 1052999.0, 143149.0, 384383.0, 533581.0, 477964.0, 213302.0, 939775.0, 246044.0, 334225.0, 189163.0)
raw_total_length <-c(2191896544.0, 624530805.0, 2580725406.0, 1944130774.0, 893918652.0, 2024397572.0, 2011571274.0, 2188276796.0, 1952563350.0, 1931161472.0, 2172094081.0, 1237466876.0, 2150531073.0, 1912053702.0, 1294604033.0, 1402423204.0, 1856329992.0, 2171464829.0, 1956565172.0, 2148003718.0, 865456175.0, 1993953446.0, 1216967811.0)
raw_GC_percent <-c(62.76, 62.14, 62.11, 61.77, 63.06, 59.81, 59.06, 62.79, 63.6, 62.67, 63.27, 63.55, 64.27, 64.26, 65.07, 63.13, 61.74, 63.3, 64.36, 64.72, 63.48, 63.19, 63.44)

#Data for qc (assembly metrics)
qc_reads <-c('9117.5.QC.MEGAHIT.assembly', '10158.8.QC.MEGAHIT.assembly', '11263.1.QC.MEGAHIT.assembly', '11306.3.QC.MEGAHIT.assembly', '11306.1.QC.MEGAHIT.assembly', '11260.6.QC.MEGAHIT.assembly', '11260.5.QC.MEGAHIT.assembly', '9053.2.QC.MEGAHIT.assembly', '9672.8.QC.MEGAHIT.assembly', '9053.4.QC.MEGAHIT.assembly', '9053.3.QC.MEGAHIT.assembly', '9053.5.QC.MEGAHIT.assembly', '9108.2.QC.MEGAHIT.assembly', '10186.3.QC.MEGAHIT.assembly', '7333.1.QC.MEGAHIT.assembly', '9117.8.QC.MEGAHIT.assembly', '10158.6.QC.MEGAHIT.assembly', '9117.7.QC.MEGAHIT.assembly', '9117.6.QC.MEGAHIT.assembly', '10186.4.QC.MEGAHIT.assembly', '9108.1.QC.MEGAHIT.assembly', '9117.4.QC.MEGAHIT.assembly', '9041.8.QC.MEGAHIT.assembly')
qc_total_contigs <-c(1177673.0, 351143.0, 1290706.0, 1020207.0, 492752.0, 1031110.0, 992335.0, 1206733.0, 1081516.0, 1047224.0, 1175033.0, 714009.0, 1199339.0, 1026124.0, 748848.0, 788426.0, 1000991.0, 1127638.0, 1104718.0, 1118259.0, 514301.0, 1060948.0, 701661.0)
qc_10k_contigs <-c(7566.0, 2011.0, 10589.0, 6974.0, 2501.0, 9310.0, 9431.0, 7594.0, 5850.0, 6738.0, 8710.0, 3172.0, 5691.0, 6688.0, 2509.0, 4428.0, 6851.0, 9186.0, 5432.0, 9182.0, 1646.0, 7553.0, 2594.0)
qc_100k_contigs <-c(37.0, 14.0, 147.0, 84.0, 46.0, 128.0, 139.0, 56.0, 35.0, 41.0, 47.0, 14.0, 15.0, 116.0, 8.0, 23.0, 33.0, 66.0, 29.0, 83.0, 11.0, 37.0, 16.0)
qc_largest_contig <-c(238481.0, 213266.0, 969148.0, 1195236.0, 536441.0, 449999.0, 1831404.0, 316258.0, 266395.0, 225864.0, 310180.0, 223153.0, 162369.0, 1172225.0, 136759.0, 384383.0, 628222.0, 477964.0, 251356.0, 909164.0, 186562.0, 404583.0, 201602.0)
qc_total_length <-c(2182276600.0, 607791927.0, 2504659553.0, 1905685212.0, 874519186.0, 1971770469.0, 1972700920.0, 2172230088.0, 1929639607.0, 1909891884.0, 2154726090.0, 1222799627.0, 2132518680.0, 1881685843.0, 1267885014.0, 1390291730.0, 1824542955.0, 2157438614.0, 1944342577.0, 2113742330.0, 857783461.0, 1982261895.0, 1186079132.0)
qc_GC_percent <-c(62.76, 62.12, 62.1, 61.76, 63.04, 59.78, 59.03, 62.79, 63.61, 62.67, 63.27, 63.56, 64.27, 64.27, 65.08, 63.12, 61.76, 63.3, 64.36, 64.73, 63.48, 63.19, 63.44)

#Data for raw
Raw_Reads <-c('9117.5_raw', '10158.8_raw', '11263.1_raw', '11306.3_raw', '11306.1_raw', '11260.6_raw', '11260.5_raw', '9108.1_raw', '9053.2_raw', '9672.8_raw', '9108.2_raw', '9053.4_raw', '9053.3_raw', '9117.4_raw', '9117.6_raw', '9117.7_raw', '9117.8_raw', '10158.6_raw', '10186.3_raw', '10186.4_raw', '7331.1_raw', '9053.5_raw', '9041.8_raw')
Raw_Bin_Counts <-c(65, 47, 139, 99, 55, 90, 115, 38, 86, 87, 69, 71, 95, 70, 62, 95, 65, 78, 85, 109, 45, 49, 52)
Raw_Mean_Completeness <-c(45.96, 58.83, 51.28, 56.54, 56.55, 63.26, 52.23, 58.47, 54.69, 53.7, 58.18, 60.32, 62.41, 65.52, 52.14, 50.97, 56.26, 57.0, 65.39, 54.89, 53.78, 53.56, 52.81)
Raw_Mean_Contamination <-c(24.03, 69.63, 67.13, 60.4, 53.21, 81.43, 55.06, 35.14, 54.71, 74.62, 54.7, 76.81, 66.86, 78.94, 67.85, 46.5, 92.78, 97.57, 121.14, 103.47, 75.71, 57.06, 65.71)
Raw_Mean_Single_Copy <-c(64.74, 103.91, 93.11, 76.5, 99.52, 112.61, 86.86, 98.89, 93.73, 76.82, 107.73, 91.71, 87.35, 81.28, 79.12, 85.37, 95.93, 83.01, 100.92, 80.27, 93.02, 89.42, 88.78)
Raw_Mean_Multi_Copy <-c(27.66, 17.26, 15.24, 20.58, 15.06, 21.89, 22.63, 19.62, 20.36, 18.31, 17.33, 15.53, 20.47, 27.25, 16.86, 9.08, 12.54, 19.53, 26.08, 17.34, 16.47, 16.8, 17.32)
Raw_Median_Completeness <-c(45.61, 65.41, 53.6, 60.38, 62.54, 74.79, 53.46, 58.27, 59.88, 59.48, 61.08, 69.28, 75.5, 77.49, 57.92, 48.54, 57.68, 60.58, 70.25, 58.65, 62.06, 55.89, 56.47)
Raw_Median_Contamination <-c(3.45, 4.12, 1.72, 3.73, 2.49, 3.44, 2.9, 2.76, 2.63, 1.94, 5.07, 2.45, 4.08, 5.22, 2.76, 1.38, 2.2, 2.22, 4.11, 2.5, 2.41, 2.16, 3.45)
Raw_Median_Single_Copy <-c(25.0, 67.5, 31.0, 36.0, 39.0, 59.0, 33.0, 65.0, 39.0, 40.0, 66.5, 39.0, 33.0, 41.0, 46.0, 38.5, 38.0, 37.5, 70.0, 31.0, 38.0, 30.0, 39.0)
Raw_Median_Multi_Copy <-c(2.0, 9.0, 2.0, 4.0, 5.0, 8.0, 4.0, 6.0, 3.0, 4.0, 10.5, 4.5, 8.0, 10.0, 6.0, 2.0, 3.0, 4.0, 10.0, 4.0, 5.0, 4.0, 5.0)

#Data for qc
QC_Reads <-c('9117.5_qc', '10158.8_qc', '11263.1_qc', '11306.3_qc', '11306.1_qc', '11260.6_qc', '11260.5_qc', '9108.1_qc', '9053.2_qc', '9672.8_qc', '9108.2_qc', '9053.4_qc', '9053.3_qc', '9117.4_qc', '9117.6_qc', '9117.7_qc', '9117.8_qc', '10158.6_qc', '10186.3_qc', '10186.4_qc', '7331.1_qc', '9053.5_qc', '9041.8_qc')
QC_Bin_Counts <-c(67, 49, 132, 100, 53, 87, 121, 37, 95, 84, 76, 68, 94, 68, 65, 103, 62, 85, 93, 106, 48, 50, 51)
QC_Mean_Completeness <-c(43.78, 57.52, 48.0, 57.83, 56.04, 61.34, 52.56, 58.38, 51.41, 50.47, 60.85, 59.81, 63.0, 65.41, 52.03, 56.04, 58.18, 60.34, 61.34, 50.62, 54.41, 54.39, 54.49)
QC_Mean_Contamination <-c(21.53, 59.82, 58.14, 60.23, 50.14, 84.7, 59.04, 39.3, 54.22, 69.66, 54.38, 67.25, 65.32, 80.07, 61.0, 47.16, 82.49, 102.6, 119.25, 100.18, 72.7, 58.44, 64.68)
QC_Mean_Single_Copy <-c(69.0, 108.55, 77.43, 80.48, 96.57, 99.11, 80.71, 101.17, 88.28, 64.81, 111.78, 94.92, 99.15, 94.04, 74.4, 101.43, 95.18, 85.24, 90.4, 81.43, 77.26, 91.34, 94.68)
QC_Mean_Multi_Copy <-c(22.82, 16.88, 15.01, 18.69, 15.69, 22.28, 23.36, 16.21, 17.32, 15.42, 20.59, 22.48, 21.1, 29.04, 15.52, 13.0, 16.99, 22.41, 22.93, 19.65, 21.53, 13.94, 15.6)
QC_Median_Completeness <-c(35.71, 60.92, 51.15, 62.19, 64.85, 74.41, 55.16, 57.41, 53.02, 53.1, 67.68, 68.56, 74.65, 79.59, 58.2, 52.42, 61.69, 72.1, 66.38, 48.38, 61.21, 58.51, 60.49)
QC_Median_Contamination <-c(1.75, 4.17, 0.93, 4.12, 2.56, 3.27, 2.59, 3.04, 2.72, 1.94, 3.3, 3.64, 3.68, 5.15, 3.4, 1.72, 3.2, 3.32, 3.51, 1.72, 3.25, 1.94, 2.21)
QC_Median_Single_Copy <-c(23.0, 77.0, 26.0, 38.5, 48.0, 37.0, 30.5, 77.0, 29.0, 32.0, 84.0, 44.0, 46.0, 51.5, 35.5, 61.0, 41.5, 36.5, 61.0, 36.0, 35.0, 34.0, 39.0)
QC_Median_Multi_Copy <-c(2.0, 9.0, 1.0, 5.5, 5.0, 6.0, 3.0, 7.0, 4.0, 2.0, 7.0, 5.0, 7.0, 9.5, 6.0, 4.0, 6.0, 4.5, 7.0, 2.0, 4.0, 3.0, 4.5)

#Packages for calculating effect size metrics
#Cliff's Delta (for non-parametric paired data)
#Hedge's g (for paired parametric data)
#Hedge's g is similar to Cohen's d but with less bias
install.packages("devtools")  ## if not already installed
devtools::install_github("mtorchiano/effsize")

#Load packages
library(effsize)

#Calculate Cliff's Delta for total contigs
res = cliff.delta(raw_total_contigs,qc_total_contigs,return.dm=TRUE)
print(res)

#Cliff's Delta
#delta estimate: 0.08506616 (negligible)
#95 percent confidence interval:
#     lower      upper 
#-0.2508704  0.4027155 

#Calculate Cliff's Delta for total contigs
res = cliff.delta(raw_total_length,qc_total_length,return.dm=TRUE)
print(res)

#Cliff's Delta
#delta estimate: 0.08884688 (negligible)
#95 percent confidence interval:
#     lower      upper 
#-0.2476729  0.4062365

#Compute Hedge's g for 
qc = qc_10k_contigs
raw = raw_10k_contigs
#cohen.d(qc,raw) #paired data
d = (c(qc,raw))
f = rep(c("qc","raw"),each=23) #change each to # of pairs
#d
#f
#cohen.d(d ~ f) #formula interface
#cohen.d(d,f) #data and factor
cohen.d(d,f,hedges.correction=TRUE)

#Hedges's g
#g estimate: -0.02252096 (negligible)
#95 percent confidence interval:
#    lower     upper 
#-0.606651  0.561609 

#Compute Hedge's g for 
qc = QC_Mean_Completeness
raw = Raw_Mean_Completeness
#cohen.d(qc,raw) #paired data
d = (c(qc,raw))
f = rep(c("qc","raw"),each=23) #change each to # of pairs
#d
#f
#cohen.d(d ~ f) #formula interface
#cohen.d(d,f) #data and factor
cohen.d(d,f,hedges.correction=TRUE)

#Hedges's g
#g estimate: -0.02252096 (negligible)
#95 percent confidence interval:
#    lower     upper 
#-0.606651  0.561609 
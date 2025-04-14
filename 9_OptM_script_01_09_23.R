##### OptM ##########

install.packages("Sizer")
install.packages("OptM")
library ("OptM")

setwd ("/Users/Katherin_Otalora/Desktop/Katherin_PhD/Treemix/treemix-1.13/25_JUL_23")

optM( "/Users/Katherin_Otalora/Desktop/Katherin_PhD/Treemix/treemix-1.13/25_JUL_23", orientagraph = F, tsv = "/Users/Katherin_Otalora/Desktop/Katherin_PhD/Treemix/treemix-1.13/25_JUL_23/Result_OptM_allsamples_ML_Tree_m+16_B500_global_se_cabrera ", method = "Evanno", ignore = NULL, thresh = 0.05, m=5))

data <- optM( "/Users/Katherin_Otalora/Desktop/Katherin_PhD/Treemix/treemix-1.13/25_JUL_23", orientagraph = F, tsv = "/Users/Katherin_Otalora/Desktop/Katherin_PhD/Treemix/treemix-1.13/25_JUL_23/Result_OptM_allsamples_ML_Tree_m+16_B500_global_se_cabrera ", method = "Evanno", ignore = NULL, thresh = 0.05, m=5)

write.table(data,"OptM_Plilfordi_K500_m5_rep10_23_JUL_23", row.names=F, quote=F)

plot_optM(data, method = "Evanno", plot = TRUE, pdf = "OptM_Plilfordi_K500_m5_rep10_23_JUL_23.pdf")

################################# 28 AGO 23 ROOT DRAGONERA ###############################

library ("OptM")

setwd ("/Users/katherin_otalora/Documents/Treemix/OUTPUT_LOOP_CORRECTO/25_JUL_23_K500_without_root")

optM( "setwd ("/Users/katherin_otalora/Documents/Treemix/OUTPUT_LOOP_CORRECTO/25_JUL_23_K500_without_root", orientagraph = F, tsv = "/Users/Katherin_Otalora/Desktop/Katherin_PhD/Treemix/treemix-1.13/28_AGO_23/Result_OptM_allsamples_ML_Tree_m+16_B500_global_se_cabrera ", method = "Evanno", ignore = NULL, thresh = 0.05, m=5))

data <- optM( "/Users/Katherin_Otalora/Desktop/Katherin_PhD/Treemix/treemix-1.13/25_JUL_23", orientagraph = F, tsv = "/Users/Katherin_Otalora/Desktop/Katherin_PhD/Treemix/treemix-1.13/28_AGO_23/Result_OptM_allsamples_ML_Tree_m+16_B500_global_se_cabrera ", method = "Evanno", ignore = NULL, thresh = 0.05, m=5)

write.table(data,"OptM_Plilfordi_K500_m5_rep10_28_AGO_23", row.names=F, quote=F)

plot_optM(data, method = "Evanno", plot = TRUE, pdf = "OptM_Plilfordi_K500_m5_rep10_28_AGO_23.pdf")
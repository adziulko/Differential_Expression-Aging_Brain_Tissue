head -32611 GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct | tail -n 1 > VWF_tpm.txt
head -47075 GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct | tail -n 1 > MBP_tpm.txt
head -21360 GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct | tail -n 1 > SEMA3A_tpm.txt
head -18248 GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct | tail -n 1 > TREM2_tpm.txt
head -45808 GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct | tail -n 1 > RBFOX3_tpm.txt
head -3402 GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct | tail -n 1 > ADAMTS4_tpm.txt
head -26695 GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct | tail -n 1 > STOM_tpm.txt
head -40226 GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct | tail -n 1 > CRABP1_tpm.txt
head -29848 GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct | tail -n 1 > DKK3_tpm.txt
head -16597 GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct | tail -n 1 > CSF1R_tpm.txt
head -9064 GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct | tail -n 1 > VHL_tpm.txt
head -37575 GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct | tail -n 1 > HIF1A_tpm.txt
head -4548 GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct | tail -n 1 > PARP1_tpm.txt
head -53044 GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct | tail -n 1 > KU70_tpm.txt
head -24964 GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct | tail -n 1 > CYC1_tpm.txt
head -11647 GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct | tail -n 1 > EIF4A2_tpm.txt

head -2731 GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct | tail -n 1 > TXNIP_tpm.txt
head -17570 GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct | tail -n 1 > HIST1H4C_tpm.txt

head -46075 GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct | tail -n 1 > DLGAP1_tpm.txt
head -24880 GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct | tail -n 1 > ARC_tpm.txt

awk 1 ID_column.txt ADAMTS4_tpm.txt CP_tpm.txt CRABP1_tpm.txt CSF1R_tpm.txt CYC1_tpm.txt DKK3_tpm.txt EIF4A2_tpm.txt HIF1A_tpm.txt KU70_tpm.txt MBP_tpm.txt PARP1_tpm.txt RBFOX3_tpm.txt SEMA3A_tpm.txt SGPP1_tpm.txt STOM_tpm.txt TREM2_tpm.txt VHL_tpm.txt VWF_tpm.txt > GOI_tpm_file.txt
awk 1 GOI_tpm_file.txt TXNIP_tpm.txt HIST1H4C_tpm.txt DLGAP1_tpm.txt ARC_tpm.txt > GOI_tpm_file_2.txt

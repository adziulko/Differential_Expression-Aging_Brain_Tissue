# Differential_Expression-Aging_Brain_Tissue
Using RNA-seq data from the open source Genome-Tissue Expression project at the Broad Institute, we are interested in measuring differential expression of genes in brain tissues amongst aging individuals.
## Getting Started

- A de-identified, open access version of the sample annotations available in dbGaP. We use this file to extract the brain tissue types of individuals. (11M)
```
$ wget https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt
```
- The above file was further simplified by extracting the desired tissue types (all brain tissue and nerve tissue). The resulting file (ID_Brain_Nerve_Tissue.txt) is used in the main python code. The result of the below code will be used for --fnb.
```
$ awk '/Brain/ || /Nerve/' GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt | awk -F'\t' '{print $1, $7}' OFS='\t'  > ID_Brain_Nerve_Tissue.txt
```

- A de-identified, open access version of the subject phenotypes available in dbGaP. We use this file to extract the ages of individuals. (20K)
```
$ wget https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt
```

- Gene TPMs. We use this file to extract the gene transcript per million (TPM) for an individual.
```
$ wget https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz
```
- The above file was subset. To subset to the genes you want, make a file containing the list of ~56,000 RNA transcripts:
```
awk '{print $2 }' GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct > gene_list.txt
```
- With the above list, extract the line number of your desired gene and run the following command in terminal (example for VWF gene):
```
head -32611 GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct | tail -n 1 > VWF_tpm.txt
```
- Once you have your gene.txt files made, join them together with the following command in terminal:
```
awk 1 VWF_tpm.txt MBP_tpm.txt SEMA3a_tpm.txt TREM2_tpm.txt > GOI_tpm_file.txt
```
- You can also use GOI_tpm_subset.sh bash script to get started with creating a subset gene list

## File Instruction

To run code to make violin-plot of TPM for each gene in different brain tissue:
```
$ python gtex_main.py --fnt GOI_tpm_file_2.txt --fnb ID_Brain_Nerve_Tissue.txt --fna GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt --tg 'TXNIP' 'CP' 'DLGAP1' 'VW' --fn 'TXNIP_CP_DLGAP1_VWF.pdf'
```

--fnt = RNA Seq gene tmp file (subset version shown above)

--fnb = Sample attributes (tissue types) file (simplified version shown above)

--fna = Sample phenotypes (ages) file (txt)

--tg = Input list of target genes (genes used in plotting. Best visulization at 4 genes)

--fn = Name of output plot file

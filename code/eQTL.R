library(MatrixEQTL);
useModel = modelLINEAR;
#5 FILE preparation
#1 SNP_genetype_file
SNP_file_name = paste("genetype.txt",sep="")
#2 Gene_file
expression_file_name=paste("norm_fpkm.txt",sep="")
#3 SNP_loc_file
snps_location_file_name = paste("snp_loc.txt",sep="")
#4 Gene_loc_file
gene_location_file_name = paste("gene_loc.txt",sep="")
#5 Covariates_file
covariates_file_name = paste("COV.txt",sep="")

#output_file_name = tempfile()
output_file_name_cis = tempfile();
output_file_name_tra = tempfile();

pvOutputThreshold_cis = 0.01;
pvOutputThreshold_tra = 0.01;

#pvOutputThreshold = 1e-20
errorCovariance = numeric();

snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

cisDist = 10000;

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 10000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 10000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
cvrt$LoadFile(covariates_file_name);
}

me = Matrix_eQTL_main(
snps = snps, 
gene = gene, 
cvrt = cvrt,
output_file_name     = output_file_name_tra,
pvOutputThreshold     = pvOutputThreshold_tra,
useModel = useModel, 
errorCovariance = errorCovariance, 
verbose = TRUE, 
output_file_name.cis = output_file_name_cis,
pvOutputThreshold.cis = pvOutputThreshold_cis,
snpspos = snpspos, 
genepos = genepos,
cisDist = cisDist,
pvalue.hist = "qqplot",
min.pv.by.genesnp = TRUE,
noFDRsaveMemory = FALSE);

unlink(output_file_name_tra);
unlink(output_file_name_cis);

#noFDRsaveMemory = FALSE);
write.table(me$cis$eqtls,file="cis-eQTL",sep="\t",quote=FALSE)
write.table(me$tran$eqtls,file="trans-eQTL",sep="\t",quote=FALSE)
write.table(me$cis$min.pv.gene, "cis.min.pv.gene.txt")

#save (me,file="PAVall.RData")
#verbose = TRUE,
plot(me)
dev.off()
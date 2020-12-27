filepath = "C:/Users/abraks/Documents/CS 548/ARCH5/human_matrix.h5"

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.10")

BiocManager::install("rhdf5")

library("rhdf5")

destination_file = filepath
extracted_expression_file = "example_expression_matrix.tsv"

# Check if gene expression file was already downloaded, if not in current directory download file form repository
if(!file.exists(destination_file)){
    print("Downloading compressed gene expression matrix.")
    url = "https://s3.amazonaws.com/mssm-seq-matrix/human_matrix.h5"
    download.file(url, destination_file, quiet = FALSE, mode = "wb")
} else{
    print("Local file already exists.")
}

# Selected samples to be extracted
samp = c("GSM1224927","GSM1066120","GSM1224923","GSM1224929","GSM1224924","GSM1066118","GSM1066119","GSM1224925","GSM1224930","GSM1872071","GSM2282084","GSM1872064","GSM1872067","GSM1704845")

# Retrieve information from compressed data
samples = h5read(destination_file, "meta/Sample_geo_accession")
# Identify columns to be extracted
sample_locations = which(samples %in% samp)

tissue = h5read(destination_file, "meta/Sample_source_name_ch1")
genes = h5read(destination_file, "meta/genes")
series = h5read(destination_file, "meta/Sample_series_id")

# extract gene expression from compressed data
expression = h5read(destination_file, "data/expression", index=list(1:length(genes), sample_locations))
H5close()
rownames(expression) = genes
colnames(expression) = samples[sample_locations]

# Print file
write.table(expression, file=extracted_expression_file, sep="\t", quote=FALSE)
print(paste0("Expression file was created at ", getwd(), "/", extracted_expression_file))

gene_names = h5read(filepath, "meta/gene_name")
genes = h5read(filepath, "meta/genes")
gene_names[1:10]


sample_meta = h5read(filepath, "meta/Sample_title")


h5ls(filepath)

expression = h5read(filepath, "data/expression", index = list(1:100,1:20))
expression

h5closeAll()



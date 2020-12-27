library(RColorBrewer)
library(gplots)

data_path = "C:/Users/abraks/Documents/CS 548/ARCH5/samples/"


#INPUTS:
#  file_list: list of csv filenames containing expression data. Will concatenate these matrices together and compute distances
preprocess_samples <- function(file_list)
{
  sample_source = c()
  # Merge all of the files into one dataframe
  df_expression = read.csv(paste0(data_path,file_list[1]),sep="\t",row.names=1,header=TRUE)
  sample_source = append(sample_source,rep(strsplit(file_list[1],"_")[[1]][1],ncol(df_expression)))
  for(filename in file_list[2:length(file_list)])
  {
    temp = read.csv(paste0(data_path,filename),sep="\t",row.names=1,header=TRUE)
    sample_source = append(sample_source,rep(strsplit(filename,"_")[[1]][1],ncol(temp)))
    df_expression = cbind(df_expression,temp)
  
  }
  return(list(df_expression,sample_source))
}
normal <- function(x){return((x-mean(x))/sqrt(sum((x-mean(x))^2)))}

#normalizes data before filtering
generate_heatmap_norm_before <- function(file_list,gene.list,samples){
  sample.list = preprocess_samples(file_list)
  df = sample.list[[1]]
  genes <- row.names(sample.list[[1]])
  df_expression = sapply(df,normal)
  row.names(df_expression) = genes
  sample_labels = sample.list[[2]]
  small.eset = df_expression[gene.list,samples]
  sample_labels = sample_labels[colnames(df_expression) %in% samples]
  hmcol <- colorRampPalette(brewer.pal(10,"RdBu"))(256)
  
  cell <- sample_labels 
  csc <- rep(hmcol[50],ncol(small.eset))
  csc[cell=='PANC1'] <- hmcol[200]
  csc[cell=='ALPHA'] <- hmcol[125]
  colnames(small.eset) <- cell
  #small.eset = as.matrix(sapply(small.eset, as.numeric))
  
  heatmap.2(small.eset,scale="row",col=hmcol, ColSideColors=csc,trace="none")
}

#normalizes data after filtering
generate_heatmap_norm_after <- function(file_list,gene.list,samples){
  sample.list = preprocess_samples(file_list)
  df = sample.list[[1]]
  genes <- row.names(sample.list[[1]])
  #df_expression = sapply(df,normal)
  row.names(df_expression) = genes
  sample_labels = sample.list[[2]]
  small.eset = df[gene.list,samples]
  
  sample_labels = sample_labels[colnames(df_expression) %in% samples]
  
  small.eset = sapply(small.eset,normal)
  row.names(small.eset) = gene.list
  
  hmcol <- colorRampPalette(brewer.pal(10,"RdBu"))(256)
  
  cell <- sample_labels 
  csc <- rep(hmcol[50],ncol(small.eset))
  csc[cell=='PANC1'] <- hmcol[200]
  csc[cell=='ALPHA'] <- hmcol[125]
  colnames(small.eset) <- cell
  #small.eset = as.matrix(sapply(small.eset, as.numeric))
  
  heatmap.2(small.eset,scale="row",col=hmcol, ColSideColors=csc,trace="none")
}

#log transforms
generate_heatmap_raw <- function(file_list,gene.list,samples){
  sample.list = preprocess_samples(file_list)
  df = sample.list[[1]]
  genes <- row.names(sample.list[[1]])
  #df_expression = sapply(df,normal)
  row.names(df_expression) = genes
  sample_labels = sample.list[[2]]
  small.eset = df[gene.list,samples]
  
  sample_labels = sample_labels[colnames(df_expression) %in% samples]
  small.eset = sapply(small.eset,as.numeric)
  row.names(small.eset) = gene.list
  
  hmcol <- colorRampPalette(brewer.pal(10,"RdBu"))(256)
  
  cell <- sample_labels 
  csc <- rep(hmcol[50],ncol(small.eset))
  csc[cell=='PANC1'] <- hmcol[200]
  csc[cell=='ALPHA'] <- hmcol[125]
  colnames(small.eset) <- cell
  #small.eset = as.matrix(sapply(small.eset, as.numeric))
  
  heatmap.2(small.eset,scale="row",col=hmcol, ColSideColors=csc,trace="none")
}


file_list = c("Ovary_expression.csv","SKOV3_expression.csv")
#row.names(df_expression) = genes
gene.list <- c('AC004824.1', 'AC006427.2', 'AC010740.1', 'AC093084.1',
               'AC130360.7', 'ADAM30', 'AF196972.4', 'AK9', 'AL355480.2', 'AP2M1',
               'ARPC1A', 'ATG16L2', 'ATP5F1', 'BCL2L15', 'BCRP1', 'C12ORF66',
               'CDR1', 'CEP128', 'CGA', 'CH17-13I23.3', 'CH507-236L23.4',
               'CNN2P10', 'CTB-179I1.1', 'CTC-503J8.2', 'CTD-2014N11.1',
               'CTD-2501O3.2', 'DIMT1', 'DLST', 'DMXL1', 'DNAJB1', 'EARS2',
               'EIF4A1P9', 'FXYD1', 'GOLGA6L4', 'GPATCH8', 'GZMAP1', 'HEMK1',
               'HMGB1P29', 'HS3ST6', 'IGHV2OR16-5', 'ITFG1', 'IYD', 'KDM4E',
               'MCRS1', 'MT-CYB', 'NOL7', 'NUTF2P4', 'OR10A3', 'OR2L13',
               'OR51A9P', 'OR9K2', 'PROSER2', 'PRSS53', 'REG4', 'RNPEP',
               'RP1-20N2.7', 'RP11-404F10.6', 'RP11-411B10.4', 'RP11-548K12.6',
               'RP11-989E6.13', 'RP3-509I19.1', 'RP3-511B24.5', 'RP4-593A12.1',
               'RPL23AP48', 'RPL23AP54', 'RPL34P27', 'RPL36AL', 'SEC62', 'SETD9',
               'SH3GL1P3', 'SLC10A6', 'SMC1B', 'STX7', 'SUGP1', 'TBX1', 'TGIF2',
               'TPPP2', 'TUBB8P2', 'TUBBP10', 'USP17L28', 'VKORC1L1', 'WFDC6',
               'WHRN', 'XCL2', 'ZC3H18', 'ZNF354C', 'ZNF449', 'ZNF705G')

samples = c('GSM742947', 'GSM979870', 'GSM1196045', 'GSM979869', 'GSM1505605', 'GSM1548004', 'GSM1548001', 'GSM1548003', 'GSM1505580', 'GSM1547996', 'GSM1548002', 'GSM1548000', 'GSM1547998', 'GSM1547997', 'GSM1547999', 'GSM2343138', 'GSM1101662', 'GSM2344226', 'GSM1010948', 'GSM2453417', 'GSM2343109', 'GSM3415785', 'GSM3415786', 'GSM3415852', 'GSM3415878', 'GSM3613418', 'GSM3613419', 'GSM3613420', 'GSM3617822', 'GSM3319032', 'GSM3319033', 'GSM3319034', 'GSM3319035', 'GSM3319036', 'GSM3319037', 'GSM3319038', 'GSM3319039', 'GSM3319040', 'GSM3319041', 'GSM3319042', 'GSM3319043', 'GSM3319044', 'GSM3319045', 'GSM3319046', 'GSM3319047', 'GSM3557967', 'GSM3557969', 'GSM3557973', 'GSM3944441', 'GSM3944442', 'GSM3944443', 'GSM3944444', 'GSM3944445', 'GSM3944446', 'GSM3944447', 'GSM3944448', 'GSM3944449', 'GSM3944450', 'GSM3944451', 'GSM3944452', 'GSM3944453', 'GSM3944454', 'GSM3944455', 'GSM3944456', 'GSM3944457', 'GSM3944458', 'GSM3944459', 'GSM3944460', 'GSM3944461', 'GSM3944462', 'GSM3944463', 'GSM3944464', 'GSM3089935', 'GSM3089936', 'GSM3089937', 'GSM3089938', 'GSM3089946', 'GSM4073729', 'GSM4073731', 'GSM4073733', 'GSM4073735', 'GSM4073737', 'GSM4073739', 'GSM4073741', 'GSM4073743', 'GSM4073745', 'GSM4073747', 'GSM4073749', 'GSM4073752', 'GSM4073754', 'GSM4073756', 'GSM4073758', 'GSM4073760', 'GSM4073762', 'GSM4073764', 'GSM4073766', 'GSM4110159', 'GSM4110160', 'GSM4110161', 'GSM1348967', 'GSM1348965', 'GSM1348966', 'GSM1348954', 'GSM1348964', 'GSM1348955', 'GSM1348963', 'GSM1348962', 'GSM1348968', 'GSM2794657', 'GSM2794658', 'GSM2794659', 'GSM2794660', 'GSM2794661', 'GSM2794662', 'GSM2631741', 'GSM2631742', 'GSM2631743', 'GSM2631744', 'GSM2631745', 'GSM2631746', 'GSM2631747', 'GSM2631748', 'GSM3613426', 'GSM3613427', 'GSM3613428', 'GSM3613429', 'GSM3613430', 'GSM3613431', 'GSM2855461', 'GSM2855462', 'GSM2855463', 'GSM2855464')


generate_heatmap_norm_before(file_list,gene.list,samples)

generate_heatmap_norm_after(file_list,gene.list,samples)

generate_heatmap_raw(file_list,gene.list,samples)



#to filter or not to filter. Make your decision by uncommenting the correct line(s)

#small.eset = sample.list[[1]]



heater = read.csv("CS 548/ARCH5/haha.csv",row.names = 1)
hmcol <- colorRampPalette(brewer.pal(10,"RdBu"))(256)

cell <- sample_labels 
csc <- rep(hmcol[50],ncol(heater))
csc[cell=='PANC1'] <- hmcol[200]
csc[cell=='ALPHA'] <- hmcol[125]
colnames(small.eset) <- cell
#small.eset = as.matrix(sapply(small.eset, as.numeric))

heatmap.2(small.eset,scale="row",col=hmcol,trace="none")

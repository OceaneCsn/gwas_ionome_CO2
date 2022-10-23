library(data.table)


load(file = "rdata/gwas_all_traits_emma_top50_astle_maf0.04_mad5.rdata")


############ export manhattan plots for IGV

get_igv_files <- function(gwas){
  traits <- unique(gwas$signSnp$Y$trait)
  for(t in traits){
    
    DT <- as.data.table(gwas$GWAResult$Y)
    dt <- DT[trait == t,]
    gwas_bed <- dt[,c("chr", "pos", "snp", "pValue")]
    gwas_bed$chr <- paste0("Chr", gwas_bed$chr)
    GWAS <- gwas_bed
    colnames(GWAS) <- c("chr", "bp", "snp", "p")
    write.table(GWAS[order(GWAS$chr, GWAS$bp),], 
                row.names = F, quote = F, sep = '\t', 
                file = paste0("results/IGV_input_", t, ".gwas"))
  }
}

get_igv_files(gwas)
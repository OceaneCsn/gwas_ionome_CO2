library(ramwas)
library(data.table)


load(file = "rdata/gwas_all_traits_emma_top50_astle_maf0.04_mad5.rdata")


# manhattan saved for all traits in gwas object
draw_manhattan <- function(gwas, trait = NULL, plot = TRUE){
  DT <- as.data.table(gwas$GWAResult$Y)
  
  if(is.null(trait)){
    for(t in unique(gwas$GWAResult$Y$trait)){
      
      dt <- DT[trait == t]
      man <- manPlotPrepare(pvalues = dt$pValue, chr = dt$chr,
                            pos = dt$pos)
      
      if(plot){
        png(filename = paste0("results/Manhattan_", t, "_trait_maf0.04_toua.png"), 
            width = 1000, height = 700)
        manPlotFast(
          man,
          ylim = NULL,
          colorSet = c("#445A67","#57838D","#B4C9C7", "#F3BFB3", "#CCADB2"),
          yaxmax = NULL,
          lwd = 3,
          axistep = 2,
          cex = 0.75)
        title(paste("Manhattan plot :", t, "trait"))
        dev.off()
      }
      else{
        manPlotFast(
          man,
          ylim = NULL,
          colorSet = c("#445A67","#57838D","#B4C9C7", "#F3BFB3", "#CCADB2"),
          yaxmax = NULL,
          lwd = 3,
          axistep = 2,
          cex = 0.75)
        title(paste("Manhattan plot :", t, "trait"))
      }
      
    }
  }
  else{
    dt <- DT[trait == t]
    man <- manPlotPrepare(pvalues = dt$pValue, chr = dt$chr,
                          pos = dt$pos)
    manPlotFast(
      man,
      ylim = NULL,
      colorSet = c("#445A67","#57838D","#B4C9C7", "#F3BFB3", "#CCADB2"),
      yaxmax = NULL,
      lwd = 3,
      axistep = 2,
      cex = 0.75)
    title(paste("Manhattan plot :", t, "trait"))
  }
}

draw_manhattan(gwas)





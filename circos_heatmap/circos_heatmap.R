###### circlize  heatmap
####### log2 scale

library(ape)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)

luc <- read.table("110RBP_data.txt",sep = "\t",header = T,row.names = 1)

mat = t(log2(luc))
luc_num <- as.numeric(as.matrix(log2(luc)))
col_fun = colorRamp2(c(-2,0,2), c("#505668","white","#C05850"))

factors = rep(letters[1], times = 111)
mat_list = list(a = mat[, factors == "a"])
dend_list = list(a = as.dendrogram(hclust(dist(t(mat_list[["a"]])),method = "ward.D")))
circos.par(cell.padding = c(0, 0,0, 0), gap.degree = 70)
circos.initialize(factors, xlim = cbind(c(0), table(factors)))
circos.track(ylim = c(0, 5.5), bg.border = NA,track.height = 0.5, panel.fun = function(x, y) {
  sector.index = CELL_META$sector.index
  m = mat_list[[sector.index]]
  dend = dend_list[[sector.index]]
  
  m2 = m[, order.dendrogram(dend)]
  col_mat = col_fun(m2)
  nr = nrow(m2)
  nc = ncol(m2)
  print(col_mat)
  print(nr)
  print(nc)
  #
  for(i in 1:nr) {
    circos.rect(1:nc - 1, rep(nr - i, nc), 
                1:nc, rep(nr - i + 1, nc), 
                border = col_mat[i, ], col = col_mat[i, ])
    circos.text(1,i,rownames(col_mat)[nr-i+1], 
                facing = "bending.inside", adj = c(1.1, 0.8), cex = 0.6)
  }
  circos.axis(sector.index = "a" ,h = 4,labels = colnames(col_mat)[1:111],labels.cex = 0.45,minor.ticks = F,
              lwd = 0 ,
              major.at = (1-0.5):(nc-0.5),
              major.tick.percentage	= 0.05,
              labels.facing = "reverse.clockwise")
})
max_height = max(sapply(dend_list, function(x) attr(x, "height")))
circos.track(ylim = c(0, max_height), bg.border = NA, track.height = 0.4, 
             panel.fun = function(x, y) {
               sector.index = get.cell.meta.data("sector.index")
               dend = dend_list[[sector.index]]
               circos.dendrogram(dend, max_height = max_height)
             })

lgd_links = Legend(at = c(-2,0,2), 
                   labels_gp = gpar(fontsize = 7),
                   col_fun = col_fun, 
                   title_position = "topleft", title = "log2(Value)",
                   title_gp = gpar(fontsize = 7))
draw(lgd_links, x = unit(1, "npc") - unit(5, "mm"), y = unit(98, "mm"), just = c("right", "top"))
circos.clear()


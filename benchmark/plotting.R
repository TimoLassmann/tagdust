library(ggplot2)
library(gplots)
library(scales)
require("reshape")

mat = read.table("barread4r.tsv",header =T,sep="\t")



mat$barcodes = paste(mat$barcodes, "Barcodes")
mat$barcodes  = factor(mat$barcodes, levels=c('8 Barcodes','24 Barcodes','48 Barcodes'))

s = subset(mat,mat$Program != "tagdustallarch")
s$Program  = factor(s$Program, levels=c('tagdust','fastx','btrim'))


ggplot(s, aes(x = simerror,y = Precision ,color=Program, fill=factor(Program))) + geom_line()  +facet_grid(scales="free", . ~ barcodes  ) +   scale_x_continuous(breaks=c(0.01,0.015,0.02,0.025,0.03), labels = expression("1","1.5","2","2.5","3")) +  xlab("Simulated Error Rate (%)")

ggsave("barread.pdf",width=6,height=4)

mat = read.table("5barread34r.tsv",header =T,sep="\t")

mat$barcodes = paste(mat$barcodes, "Barcodes")
mat$barcodes  = factor(mat$barcodes, levels=c('8 Barcodes','24 Barcodes','48 Barcodes'))
s = subset(mat,mat$Program != "tagdustallarch")
s$Program  = factor(s$Program, levels=c('tagdust','fastx','btrim'))



ggplot(s, aes(x = simerror,y = Precision ,color=Program, fill=factor(Program))) + geom_line()  +facet_grid(scales="free", . ~ barcodes  ) +   scale_x_continuous(breaks=c(0.01,0.015,0.02,0.025,0.03), labels = expression("1","1.5","2","2.5","3")) +  xlab("Simulated Error Rate (%)")

ggsave("5barread3.pdf",width=6,height=4)




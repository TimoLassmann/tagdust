library(ggplot2)
library(gplots)
library(scales)
require("reshape")

yaxis_label_format <- function(x) {
	lab <- sprintf('%0.2f', x) # Format the strings as HH:MM:SS
}


mat = read.table("barread_4nt_4r.tsv",header =T,sep="\t")



mat$barcodes = paste(mat$barcodes, "Barcodes")
mat$barcodes  = factor(mat$barcodes, levels=c('8 Barcodes','24 Barcodes','48 Barcodes'))

s = subset(mat,mat$Program != "tagdustallarch")
s$Program  = factor(s$Program, levels=c('tagdust','fastx','btrim'))

m = melt(s, id=c("Program","simerror","barcodes"),measure.vars = c("Recall", "Precision"))
ggplot(m, aes(x = simerror,y = value ,color=Program, fill=factor(Program))) + geom_line()  +facet_grid(scales="free", variable ~ barcodes  ) +   scale_x_continuous(breaks=c(0.01,0.015,0.02,0.025,0.03), labels = expression("1","1.5","2","2.5","3")) +  xlab("Simulated Error Rate (%)") + scale_y_continuous(labels=yaxis_label_format)+ ylab("")

ggsave("barread_4nt.pdf",width=8,height=4,dpi = 300)




mat = read.table("barread_6nt_4r.tsv",header =T,sep="\t")

mat$barcodes = paste(mat$barcodes, "Barcodes")
mat$barcodes  = factor(mat$barcodes, levels=c('8 Barcodes','24 Barcodes','48 Barcodes'))

s = subset(mat,mat$Program != "tagdustallarch")
s$Program  = factor(s$Program, levels=c('tagdust','fastx','btrim'))
m = melt(s, id=c("Program","simerror","barcodes"),measure.vars = c("Recall", "Precision"))
ggplot(m, aes(x = simerror,y = value ,color=Program, fill=factor(Program))) + geom_line()  +facet_grid(scales="free", variable ~ barcodes  ) +   scale_x_continuous(breaks=c(0.01,0.015,0.02,0.025,0.03), labels = expression("1","1.5","2","2.5","3")) +  xlab("Simulated Error Rate (%)") + scale_y_continuous(labels=yaxis_label_format)+ ylab("")

ggsave("barread_6nt.pdf",width=8,height=4,dpi = 300)





mat = read.table("5barread3_4nt_4r.tsv",header =T,sep="\t")
mat$barcodes = paste(mat$barcodes, "Barcodes")
mat$barcodes  = factor(mat$barcodes, levels=c('8 Barcodes','24 Barcodes','48 Barcodes'))
s = subset(mat,mat$Program != "tagdustallarch")
s$Program  = factor(s$Program, levels=c('tagdust','fastx','btrim'))
m = melt(s, id=c("Program","simerror","barcodes"),measure.vars = c("Recall", "Precision"))
ggplot(m, aes(x = simerror,y = value ,color=Program, fill=factor(Program))) + geom_line()  +facet_grid(scales="free", variable ~ barcodes  ) +   scale_x_continuous(breaks=c(0.01,0.015,0.02,0.025,0.03), labels = expression("1","1.5","2","2.5","3")) +  xlab("Simulated Error Rate (%)") + scale_y_continuous(labels=yaxis_label_format) + ylab("")

ggsave("5barread3_4nt.pdf",width=8,height=4,dpi = 300)


mat = read.table("5barread3_6nt_4r.tsv",header =T,sep="\t")

mat$barcodes = paste(mat$barcodes, "Barcodes")
mat$barcodes  = factor(mat$barcodes, levels=c('8 Barcodes','24 Barcodes','48 Barcodes'))
s = subset(mat,mat$Program != "tagdustallarch")
s$Program  = factor(s$Program, levels=c('tagdust','fastx','btrim'))
m = melt(s, id=c("Program","simerror","barcodes"),measure.vars = c("Recall", "Precision"))
ggplot(m, aes(x = simerror,y = value ,color=Program, fill=factor(Program))) + geom_line()  +facet_grid(scales="free", variable ~ barcodes  ) +   scale_x_continuous(breaks=c(0.01,0.015,0.02,0.025,0.03), labels = expression("1","1.5","2","2.5","3")) +  xlab("Simulated Error Rate (%)") + scale_y_continuous(labels=yaxis_label_format)+ ylab("")

ggsave("5barread3_6nt.pdf",width=8,height=4,dpi = 300)




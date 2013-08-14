 

library(ggplot2)
library(gplots)

td10 = read.table ("tagdust_benchmark10_norm.csv", header = F)
td10[,8] = sprintf("N:%2d L:%2d" ,td10[,1],td10[,2])
#td10[,8] = paste(td10[,1],td10[,2])
td10[,9] = 10000 * td10[,1] + td10[,2];


td10[,10] ="TagDust"

f1 = read.table("fastx_benchmark_1_norm.csv",header = F)
#f1[,8] = paste(f1[,1],f1[,2])
f1[,8] = sprintf("N:%2d L:%2d" ,f1[,1],f1[,2])
f1[,9] = 10000 * f1[,1] + f1[,2];

f1[,10] ="Fastx"

x = rbind(td10,f1)
x[,11] = 0;

x[,12] = 0;
x[,13] = 0;
colnames(x) = c("NumBar","BarLen","Error Rate","InDelRate","Extracted","Not","MisExtracted","Name","INDEXX","Program","Total","PrintedLabel","RANGE")



for (i in 1:nrow(x)){
	
	x$Total[i] = x$Extracted[i] + x$MisExtracted[i] ;
	
	if(x$Total[i] == 0){
		x$Total[i] = NA;
		x$PrintedLabel[i] =  sprintf("%1.1f%%", 0 )
	}else{
		x$PrintedLabel[i] =  sprintf("%1.1f%%", x$Total[i] /10000)
	}
}



x$RANGE <- cut(x$Total/10000,breaks = c(0,70,80,90.05,91,92,93,94,95,96,97,98,99,101 ),right = FALSE, label = c(1:13) )


p <- ggplot(x, aes(reorder(Name, INDEXX), `Error Rate`, label = PrintedLabel )) + geom_tile(aes(fill = RANGE))
p+scale_fill_manual(breaks=c(1:13),values = union (colorRampPalette(c("white","lightgreen"))(3),colorRampPalette(c("white","firebrick3"))(10)) , na.value = "grey50") + theme(legend.position = "none", axis.ticks = element_blank(),strip.text.y = element_text(size = 15), axis.title.x = element_blank(), axis.text.x = element_text(size = 10 , angle = 310, hjust = 0, colour = "black"))+facet_grid( Program ~ . )+ geom_text(aes(size=10))+ scale_y_continuous(breaks=c(0.0,0.01,0.02,0.03,0.04,0.05), labels = expression("0%","1%","2%","3%","4%","5%"))



ggsave("SimQ10F1NORMExtracted.pdf", width=8, height=5,dpi = 300);


x$PrintedLabel =  sprintf("%1.1f%%", x$MisExtracted /10000)


x$RANGE <- cut(x$MisExtracted /(x$Extracted[i] + x$MisExtracted[i]) *100,breaks = c(0.0,0.1,0.2,0.3,0.4,0.5,1,2,5,10,40),right = FALSE, label = c(1,2,3,4,5,6,7,8,9,10) )


for (i in 1:nrow(x)){
	
	x$Total[i] = x$Extracted[i] + x$MisExtracted[i] ;
	
	if(x$Total[i] == 0){
		x$RANGE[i] = NA;
	}
}





p <- ggplot(x, aes(reorder(Name, INDEXX), `Error Rate`,  label = PrintedLabel )) + geom_tile(aes(fill = RANGE))
p+scale_fill_manual(breaks=c(1,2,3,4,5,6,7,8,9,10),values =  colorRampPalette(c("white","firebrick3"))(10), na.value = "grey50") +  theme(legend.position = "none", axis.ticks = element_blank(),strip.text.y = element_text(size = 15), axis.title.x = element_blank(), axis.text.x = element_text(size = 10 , angle = 310, hjust = 0, colour = "black"))+facet_grid( Program ~ . )+ geom_text(aes(size=10))+ scale_y_continuous(breaks=c(0.0,0.01,0.02,0.03,0.04,0.05), labels = expression("0%","1%","2%","3%","4%","5%"))

ggsave("SimQ10F1NORMErrors.pdf", width=8, height=5,dpi = 300);





td10 = read.table ("tagdust_benchmark10_INDEL.csv", header = F)
td10[,8] = sprintf("N:%2d L:%2d" ,td10[,1],td10[,2])
#td10[,8] = paste(td10[,1],td10[,2])
td10[,9] = 10000 * td10[,1] + td10[,2];


td10[,10] ="TagDust"

f1 = read.table("fastx_benchmark_1_INDEL.csv",header = F)
#f1[,8] = paste(f1[,1],f1[,2])
f1[,8] = sprintf("N:%2d L:%2d" ,f1[,1],f1[,2])
f1[,9] = 10000 * f1[,1] + f1[,2];

f1[,10] ="Fastx"


x = rbind(td10,f1)
x[,11] = 0;

x[,12] = 0;
x[,13] = 0;
colnames(x) = c("NumBar","BarLen","Error Rate","InDelRate","Extracted","Not","MisExtracted","Name","INDEXX","Program","Total","PrintedLabel","RANGE")



for (i in 1:nrow(x)){
	
	x$Total[i] = x$Extracted[i] + x$MisExtracted[i] ;
	
	if(x$Total[i] == 0){
		x$Total[i] = NA;
		x$PrintedLabel[i] =  sprintf("%1.1f%%", 0 )
	}else{
		x$PrintedLabel[i] =  sprintf("%1.1f%%", x$Total[i] /10000)
	}
}



x$RANGE <- cut(x$Total/10000,breaks = c(0,70,80,90.05,91,92,93,94,95,96,97,98,99,101 ),right = FALSE, label = c(1:13) )


p <- ggplot(x, aes(reorder(Name, INDEXX), `Error Rate`, label = PrintedLabel )) + geom_tile(aes(fill = RANGE))
p+scale_fill_manual(breaks=c(1:13),values = union (colorRampPalette(c("white","lightgreen"))(3),colorRampPalette(c("white","firebrick3"))(10)) , na.value = "grey50") + theme(legend.position = "none",strip.text.y = element_text(size = 15), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.text.x = element_text(size = 10 , angle = 310, hjust = 0, colour = "black"))+facet_grid( Program ~ . )+ geom_text(aes(size=10))+ scale_y_continuous(breaks=c(0.0,0.01,0.02,0.03,0.04,0.05), labels = expression("0%","1%","2%","3%","4%","5%"))



ggsave("SimQ10F1INDELExtracted.pdf", width=8, height=5,dpi = 300);


x$PrintedLabel =  sprintf("%1.1f%%", x$MisExtracted /10000)


x$RANGE <- cut(x$MisExtracted /(x$Extracted[i] + x$MisExtracted[i]) *100,breaks = c(0.0,0.1,0.2,0.3,0.4,0.5,1,2,5,10,40),right = FALSE, label = c(1,2,3,4,5,6,7,8,9,10) )


for (i in 1:nrow(x)){
	
	x$Total[i] = x$Extracted[i] + x$MisExtracted[i] ;
	
	if(x$Total[i] == 0){
		x$RANGE[i] = NA;
	}
}





p <- ggplot(x, aes(reorder(Name, INDEXX), `Error Rate`,  label = PrintedLabel )) + geom_tile(aes(fill = RANGE))
p+scale_fill_manual(breaks=c(1,2,3,4,5,6,7,8,9,10),values =  colorRampPalette(c("white","firebrick3"))(10), na.value = "grey50") +  theme(legend.position = "none", axis.ticks = element_blank(),strip.text.y = element_text(size = 15), axis.title.x = element_blank(), axis.text.x = element_text(size = 10 , angle = 310, hjust = 0, colour = "black"))+facet_grid( Program ~ . )+ geom_text(aes(size=10)) + scale_y_continuous(breaks=c(0.0,0.01,0.02,0.03,0.04,0.05), labels = expression("0%","1%","2%","3%","4%","5%"))
ggsave("SimQ10F1INDELErrors.pdf", width=8, height=5,dpi = 300);











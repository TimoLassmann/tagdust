library(ggplot2)
library(gplots)
library(scales)

mat = read.table("../data/NsummaryQ20.csv",header =F )


s = subset(mat,mat[,9] == "agggaggacgatgcgg" & mat[,11] == 0 & mat[,16] == 0.1)
names(s)  = c("date","Program","good","rand","wrong","error")
names(s)[15] = "errorRate"

ggplot(s, aes(x = errorRate,y = good+rand+wrong , fill=factor(Program))) + geom_bar(position="dodge",stat="identity") +scale_x_continuous(labels = percent_format(),breaks=seq(0, 0.05, 0.01))  + scale_y_continuous(limits=c(80000,91000),oob = rescale_none) + ylab("Number of Extracted Reads") + xlab("Simulated Error Rate") + theme(legend.title = element_blank() ,legend.justification=c(0,0), legend.position=c(0.1,0.1),legend.background = element_rect(colour = "black"),legend.text = element_text(size = 12 ), axis.ticks = element_blank(),strip.text.y = element_text(size = 15), axis.text.x = element_text(size = 10, hjust = 0, colour = "black"))

ggsave("Q20ReadTrimming1extracted.pdf" ,width =16, height =8, units = c("cm"), dpi = 300);


ggplot(s, aes(x = errorRate,y = error , fill=factor(Program))) + geom_bar(position="dodge",stat="identity") +scale_x_continuous(labels = percent_format(),breaks=seq(0, 0.05, 0.01))  + scale_y_continuous(breaks=seq(0, 0.07, 0.01),limits=c(0,0.075),labels = percent_format(),oob = rescale_none) + ylab("Error Rate in Extracted Reads") + xlab("Simulated Error Rate") + theme(legend.title = element_blank() ,legend.justification=c(0,0), legend.position=c(0.1,0.5),legend.background = element_rect(colour = "black"),legend.text = element_text(size = 12 ), axis.ticks = element_blank(),strip.text.y = element_text(size = 15), axis.text.x = element_text(size = 10, hjust = 0, colour = "black"))
ggsave("Q20ReadTrimming1errorrateinreads.pdf" , width =16, height =8, units = c("cm"), dpi = 300);



mat = read.table("../data/NsummaryQ10.csv",header =F )


s = subset(mat,mat[,9] == "agggaggacgatgcgg" & mat[,11] == 0 & mat[,16] == 0.1)
names(s)  = c("date","Program","good","rand","wrong","error")
names(s)[15] = "errorRate"

ggplot(s, aes(x = errorRate,y = good+rand+wrong , fill=factor(Program))) + geom_bar(position="dodge",stat="identity") +scale_x_continuous(labels = percent_format(),breaks=seq(0, 0.05, 0.01))  + scale_y_continuous(limits=c(80000,91000),oob = rescale_none) + ylab("Number of Extracted Reads") + xlab("Simulated Error Rate") + theme(legend.title = element_blank() ,legend.justification=c(0,0), legend.position=c(0.1,0.1),legend.background = element_rect(colour = "black"),legend.text = element_text(size = 12 ), axis.ticks = element_blank(),strip.text.y = element_text(size = 15), axis.text.x = element_text(size = 10, hjust = 0, colour = "black"))

ggsave("Q10ReadTrimming1extracted.pdf"  ,width =16, height =8, units = c("cm"), dpi = 300);


ggplot(s, aes(x = errorRate,y = error , fill=factor(Program))) + geom_bar(position="dodge",stat="identity") +scale_x_continuous(labels = percent_format(),breaks=seq(0, 0.05, 0.01))  + scale_y_continuous(breaks=seq(0, 0.07, 0.01),limits=c(0,0.075),labels = percent_format(),oob = rescale_none) + ylab("Error Rate in Extracted Reads") + xlab("Simulated Error Rate") + theme(legend.title = element_blank() ,legend.justification=c(0,0), legend.position=c(0.1,0.5),legend.background = element_rect(colour = "black"),legend.text = element_text(size = 12 ), axis.ticks = element_blank(),strip.text.y = element_text(size = 15), axis.text.x = element_text(size = 10, hjust = 0, colour = "black"))
ggsave("Q10ReadTrimming1errorrateinreads.pdf"  ,width =16, height =8, units = c("cm"), dpi = 300);




mat = read.table("../data/NsummaryQ20.csv",header =F )
s = subset(mat,mat[,9] == "agggaggacgatgcgg" & mat[,11] == 0 & mat[,16] == 0.5)
names(s)  = c("date","Program","good","rand","wrong","error")
names(s)[15] = "errorRate"

ggplot(s, aes(x = errorRate,y = good+rand+wrong , fill=factor(Program))) + geom_bar(position="dodge",stat="identity") +scale_x_continuous(labels = percent_format(),breaks=seq(0, 0.05, 0.01))  + scale_y_continuous(limits=c(80000,91000),oob = rescale_none) + ylab("Number of Extracted Reads") + xlab("Simulated Error Rate") + theme(legend.title = element_blank() ,legend.justification=c(0,0), legend.position=c(0.1,0.1),legend.background = element_rect(colour = "black"),legend.text = element_text(size = 12 ), axis.ticks = element_blank(),strip.text.y = element_text(size = 15), axis.text.x = element_text(size = 10, hjust = 0, colour = "black"))

ggsave("Q20ReadTrimming5extracted.pdf"  ,width =16, height =8, units = c("cm"), dpi = 300);


ggplot(s, aes(x = errorRate,y = error , fill=factor(Program))) + geom_bar(position="dodge",stat="identity") +scale_x_continuous(labels = percent_format(),breaks=seq(0, 0.05, 0.01))  + scale_y_continuous(breaks=seq(0, 0.07, 0.01),limits=c(0,0.075),labels = percent_format(),oob = rescale_none) + ylab("Error Rate in Extracted Reads") + xlab("Simulated Error Rate") + theme(legend.title = element_blank() ,legend.justification=c(0,0), legend.position=c(0.1,0.7),legend.background = element_rect(colour = "black"),legend.text = element_text(size = 12 ), axis.ticks = element_blank(),strip.text.y = element_text(size = 15), axis.text.x = element_text(size = 10, hjust = 0, colour = "black"))
ggsave("Q20ReadTrimming5errorrateinreads.pdf" ,width =16, height =8, units = c("cm"), dpi = 300 );



mat = read.table("../data/NsummaryQ10.csv",header =F )


s = subset(mat,mat[,9] == "agggaggacgatgcgg" & mat[,11] == 0 & mat[,16] == 0.5)
names(s)  = c("date","Program","good","rand","wrong","error")
names(s)[15] = "errorRate"

ggplot(s, aes(x = errorRate,y = good+rand+wrong , fill=factor(Program))) + geom_bar(position="dodge",stat="identity") +scale_x_continuous(labels = percent_format(),breaks=seq(0, 0.05, 0.01))  + scale_y_continuous(limits=c(80000,91000),oob = rescale_none) + ylab("Number of Extracted Reads") + xlab("Simulated Error Rate") + theme(legend.title = element_blank() ,legend.justification=c(0,0), legend.position=c(0.1,0.1),legend.background = element_rect(colour = "black"),legend.text = element_text(size = 12 ), axis.ticks = element_blank(),strip.text.y = element_text(size = 15), axis.text.x = element_text(size = 10, hjust = 0, colour = "black"))

ggsave("Q10ReadTrimming5extracted.pdf" ,width =16, height =8, units = c("cm"), dpi = 300 );


ggplot(s, aes(x = errorRate,y = error , fill=factor(Program))) + geom_bar(position="dodge",stat="identity") +scale_x_continuous(labels = percent_format(),breaks=seq(0, 0.05, 0.01))  + scale_y_continuous(breaks=seq(0, 0.07, 0.01),limits=c(0,0.075),labels = percent_format(),oob = rescale_none) + ylab("Error Rate in Extracted Reads") + xlab("Simulated Error Rate") + theme(legend.title = element_blank() ,legend.justification=c(0,0), legend.position=c(0.1,0.7),legend.background = element_rect(colour = "black"),legend.text = element_text(size = 12 ), axis.ticks = element_blank(),strip.text.y = element_text(size = 15), axis.text.x = element_text(size = 10, hjust = 0, colour = "black"))
ggsave("Q10ReadTrimming5errorrateinreads.pdf" ,width =16, height =8, units = c("cm"), dpi = 300 );









mat = read.table("../data/NsummaryQ20.csv",header =F )
s = subset(mat,is.na(mat[,9]) & mat[,16] == 0.1 )

x = cbind(s[,11],s[,12],s[,15],s[,16] ,s[,3] + s[,4], 100000 - (s[,3] + s[,4]), s[,4]+s[5],sprintf("N:%2d L:%2d",s[,11],s[,12]  ),10000 * s[,11] + s[,12],s[,2] , s[,3]+s[,4] ) 

x[,12] = 0;
x[,13] = 0;

colnames(x) = c("NumBar","BarLen","Error Rate","InDelRate","Extracted","Not","MisExtracted","Name","INDEXX","Program","Total","PrintedLabel","RANGE")


for (i in 1:nrow(x)){
	
	
	if(x$Total[i] == 0){
		x$Total[i] = NA;
		x$PrintedLabel[i] =  sprintf("%1.1f%%", 0 )
	}else{
		x$PrintedLabel[i] =  sprintf("%1.1f%%", x$Total[i] /1000)
	}
}



x$RANGE <- cut(x$Total/1000,breaks = c(0,50,80,90.05,91,92,93,94,95,96,97,98,99,101 ),right = FALSE, label = c(1:13) )


p <- ggplot(x, aes(Name, `Error Rate`, label = PrintedLabel )) + geom_tile(aes(fill = RANGE))
p+scale_fill_manual(breaks=c(1:13),values = union (colorRampPalette(c("white","lightgreen"))(3),colorRampPalette(c("white","firebrick3"))(10)) , na.value = "grey50") + theme(legend.position = "none", axis.ticks = element_blank(),strip.text.y = element_text(size = 15), axis.title.x = element_blank(), axis.text.x = element_text(size = 10 , angle = 310, hjust = 0, colour = "black"))+facet_grid( Program ~ . )+ geom_text(aes(size=10))+ scale_y_continuous(breaks=c(0.0,0.01,0.02,0.03,0.04,0.05), labels = expression("0%","1%","2%","3%","4%","5%"))
ggsave("Q20onlyBarcode1totalextracted.pdf" );

x$PrintedLabel =  sprintf("%1.1f%%", x$MisExtracted /1000)


x$RANGE <- cut(x$MisExtracted /(x$Extracted[i] + x$MisExtracted[i]) *100,breaks = c(0.0,0.1,0.2,0.3,0.4,0.5,1,2,5,10,40),right = FALSE, label = c(1,2,3,4,5,6,7,8,9,10) )


for (i in 1:nrow(x)){
	
	x$Total[i] = x$Extracted[i] + x$MisExtracted[i] ;
	
	if(x$Total[i] == 0){
		x$RANGE[i] = NA;
	}
}





p <- ggplot(x, aes(Name, `Error Rate`,  label = PrintedLabel )) + geom_tile(aes(fill = RANGE))
p+scale_fill_manual(breaks=c(1,2,3,4,5,6,7,8,9,10),values =  colorRampPalette(c("white","firebrick3"))(10), na.value = "grey50") +  theme(legend.position = "none", axis.ticks = element_blank(),strip.text.y = element_text(size = 15), axis.title.x = element_blank(), axis.text.x = element_text(size = 10 , angle = 310, hjust = 0, colour = "black"))+facet_grid( Program ~ . )+ geom_text(aes(size=10))+ scale_y_continuous(breaks=c(0.0,0.01,0.02,0.03,0.04,0.05), labels = expression("0%","1%","2%","3%","4%","5%"))

ggsave("Q20onlyBarcode1misextracted.pdf" );


mat = read.table("../data/NsummaryQ10.csv",header =F )
s = subset(mat,is.na(mat[,9]) & mat[,16] == 0.1 )

x = cbind(s[,11],s[,12],s[,15],s[,16] ,s[,3] + s[,4], 100000 - (s[,3] + s[,4]), s[,4]+s[5],sprintf("N:%2d L:%2d",s[,11],s[,12]  ),10000 * s[,11] + s[,12],s[,2] , s[,3]+s[,4] )

x[,12] = 0;
x[,13] = 0;

colnames(x) = c("NumBar","BarLen","Error Rate","InDelRate","Extracted","Not","MisExtracted","Name","INDEXX","Program","Total","PrintedLabel","RANGE")


for (i in 1:nrow(x)){
	
	
	if(x$Total[i] == 0){
		x$Total[i] = NA;
		x$PrintedLabel[i] =  sprintf("%1.1f%%", 0 )
	}else{
		x$PrintedLabel[i] =  sprintf("%1.1f%%", x$Total[i] /1000)
	}
}



x$RANGE <- cut(x$Total/1000,breaks = c(0,50,80,90.05,91,92,93,94,95,96,97,98,99,101 ),right = FALSE, label = c(1:13) )


p <- ggplot(x, aes(Name, `Error Rate`, label = PrintedLabel )) + geom_tile(aes(fill = RANGE))
p+scale_fill_manual(breaks=c(1:13),values = union (colorRampPalette(c("white","lightgreen"))(3),colorRampPalette(c("white","firebrick3"))(10)) , na.value = "grey50") + theme(legend.position = "none", axis.ticks = element_blank(),strip.text.y = element_text(size = 15), axis.title.x = element_blank(), axis.text.x = element_text(size = 10 , angle = 310, hjust = 0, colour = "black"))+facet_grid( Program ~ . )+ geom_text(aes(size=10))+ scale_y_continuous(breaks=c(0.0,0.01,0.02,0.03,0.04,0.05), labels = expression("0%","1%","2%","3%","4%","5%"))
ggsave("Q10onlyBarcode1totalextracted.pdf" );

x$PrintedLabel =  sprintf("%1.1f%%", x$MisExtracted /1000)


x$RANGE <- cut(x$MisExtracted /(x$Extracted[i] + x$MisExtracted[i]) *100,breaks = c(0.0,0.1,0.2,0.3,0.4,0.5,1,2,5,10,40),right = FALSE, label = c(1,2,3,4,5,6,7,8,9,10) )


for (i in 1:nrow(x)){
	
	x$Total[i] = x$Extracted[i] + x$MisExtracted[i] ;
	
	if(x$Total[i] == 0){
		x$RANGE[i] = NA;
	}
}





p <- ggplot(x, aes(Name, `Error Rate`,  label = PrintedLabel )) + geom_tile(aes(fill = RANGE))
p+scale_fill_manual(breaks=c(1,2,3,4,5,6,7,8,9,10),values =  colorRampPalette(c("white","firebrick3"))(10), na.value = "grey50") +  theme(legend.position = "none", axis.ticks = element_blank(),strip.text.y = element_text(size = 15), axis.title.x = element_blank(), axis.text.x = element_text(size = 10 , angle = 310, hjust = 0, colour = "black"))+facet_grid( Program ~ . )+ geom_text(aes(size=10))+ scale_y_continuous(breaks=c(0.0,0.01,0.02,0.03,0.04,0.05), labels = expression("0%","1%","2%","3%","4%","5%"))

ggsave("Q10onlyBarcode1misextracted.pdf" );












mat = read.table("../data/NsummaryQ20.csv",header =F )
s = subset(mat,is.na(mat[,9]) & mat[,16] == 0.5 )

x = cbind(s[,11],s[,12],s[,15],s[,16] ,s[,3] + s[,4], 100000 - (s[,3] + s[,4]), s[,4]+s[5],sprintf("N:%2d L:%2d",s[,11],s[,12]  ),10000 * s[,11] + s[,12],s[,2] , s[,3]+s[,4] )

x[,12] = 0;
x[,13] = 0;

colnames(x) = c("NumBar","BarLen","Error Rate","InDelRate","Extracted","Not","MisExtracted","Name","INDEXX","Program","Total","PrintedLabel","RANGE")


for (i in 1:nrow(x)){
	
	
	if(x$Total[i] == 0){
		x$Total[i] = NA;
		x$PrintedLabel[i] =  sprintf("%1.1f%%", 0 )
	}else{
		x$PrintedLabel[i] =  sprintf("%1.1f%%", x$Total[i] /1000)
	}
}



x$RANGE <- cut(x$Total/1000,breaks = c(0,50,80,90.05,91,92,93,94,95,96,97,98,99,101 ),right = FALSE, label = c(1:13) )


p <- ggplot(x, aes(Name, `Error Rate`, label = PrintedLabel )) + geom_tile(aes(fill = RANGE))
p+scale_fill_manual(breaks=c(1:13),values = union (colorRampPalette(c("white","lightgreen"))(3),colorRampPalette(c("white","firebrick3"))(10)) , na.value = "grey50") + theme(legend.position = "none", axis.ticks = element_blank(),strip.text.y = element_text(size = 15), axis.title.x = element_blank(), axis.text.x = element_text(size = 10 , angle = 310, hjust = 0, colour = "black"))+facet_grid( Program ~ . )+ geom_text(aes(size=10))+ scale_y_continuous(breaks=c(0.0,0.01,0.02,0.03,0.04,0.05), labels = expression("0%","1%","2%","3%","4%","5%"))
ggsave("Q20onlyBarcode5totalextracted.pdf" );

x$PrintedLabel =  sprintf("%1.1f%%", x$MisExtracted /1000)


x$RANGE <- cut(x$MisExtracted /(x$Extracted[i] + x$MisExtracted[i]) *100,breaks = c(0.0,0.1,0.2,0.3,0.4,0.5,1,2,5,10,40),right = FALSE, label = c(1,2,3,4,5,6,7,8,9,10) )


for (i in 1:nrow(x)){
	
	x$Total[i] = x$Extracted[i] + x$MisExtracted[i] ;
	
	if(x$Total[i] == 0){
		x$RANGE[i] = NA;
	}
}





p <- ggplot(x, aes(Name, `Error Rate`,  label = PrintedLabel )) + geom_tile(aes(fill = RANGE))
p+scale_fill_manual(breaks=c(1,2,3,4,5,6,7,8,9,10),values =  colorRampPalette(c("white","firebrick3"))(10), na.value = "grey50") +  theme(legend.position = "none", axis.ticks = element_blank(),strip.text.y = element_text(size = 15), axis.title.x = element_blank(), axis.text.x = element_text(size = 10 , angle = 310, hjust = 0, colour = "black"))+facet_grid( Program ~ . )+ geom_text(aes(size=10))+ scale_y_continuous(breaks=c(0.0,0.01,0.02,0.03,0.04,0.05), labels = expression("0%","1%","2%","3%","4%","5%"))

ggsave("Q20onlyBarcode5misextracted.pdf" );


mat = read.table("../data/NsummaryQ10.csv",header =F )
s = subset(mat,is.na(mat[,9]) & mat[,16] == 0.5 )

x = cbind(s[,11],s[,12],s[,15],s[,16] ,s[,3] + s[,4], 100000 - (s[,3] + s[,4]), s[,4]+s[5],sprintf("N:%2d L:%2d",s[,11],s[,12]  ),10000 * s[,11] + s[,12],s[,2] , s[,3]+s[,4] )

x[,12] = 0;
x[,13] = 0;

colnames(x) = c("NumBar","BarLen","Error Rate","InDelRate","Extracted","Not","MisExtracted","Name","INDEXX","Program","Total","PrintedLabel","RANGE")


for (i in 1:nrow(x)){
	
	
	if(x$Total[i] == 0){
		x$Total[i] = NA;
		x$PrintedLabel[i] =  sprintf("%1.1f%%", 0 )
	}else{
		x$PrintedLabel[i] =  sprintf("%1.1f%%", x$Total[i] /1000)
	}
}



x$RANGE <- cut(x$Total/1000,breaks = c(0,50,80,90.05,91,92,93,94,95,96,97,98,99,101 ),right = FALSE, label = c(1:13) )


p <- ggplot(x, aes(Name, `Error Rate`, label = PrintedLabel )) + geom_tile(aes(fill = RANGE))
p+scale_fill_manual(breaks=c(1:13),values = union (colorRampPalette(c("white","lightgreen"))(3),colorRampPalette(c("white","firebrick3"))(10)) , na.value = "grey50") + theme(legend.position = "none", axis.ticks = element_blank(),strip.text.y = element_text(size = 15), axis.title.x = element_blank(), axis.text.x = element_text(size = 10 , angle = 310, hjust = 0, colour = "black"))+facet_grid( Program ~ . )+ geom_text(aes(size=10))+ scale_y_continuous(breaks=c(0.0,0.01,0.02,0.03,0.04,0.05), labels = expression("0%","1%","2%","3%","4%","5%"))
ggsave("Q10onlyBarcode5totalextracted.pdf" );

x$PrintedLabel =  sprintf("%1.1f%%", x$MisExtracted /1000)


x$RANGE <- cut(x$MisExtracted /(x$Extracted[i] + x$MisExtracted[i]) *100,breaks = c(0.0,0.1,0.2,0.3,0.4,0.5,1,2,5,10,40),right = FALSE, label = c(1,2,3,4,5,6,7,8,9,10) )


for (i in 1:nrow(x)){
	
	x$Total[i] = x$Extracted[i] + x$MisExtracted[i] ;
	
	if(x$Total[i] == 0){
		x$RANGE[i] = NA;
	}
}





p <- ggplot(x, aes(Name, `Error Rate`,  label = PrintedLabel )) + geom_tile(aes(fill = RANGE))
p+scale_fill_manual(breaks=c(1,2,3,4,5,6,7,8,9,10),values =  colorRampPalette(c("white","firebrick3"))(10), na.value = "grey50") +  theme(legend.position = "none", axis.ticks = element_blank(),strip.text.y = element_text(size = 15), axis.title.x = element_blank(), axis.text.x = element_text(size = 10 , angle = 310, hjust = 0, colour = "black"))+facet_grid( Program ~ . )+ geom_text(aes(size=10))+ scale_y_continuous(breaks=c(0.0,0.01,0.02,0.03,0.04,0.05), labels = expression("0%","1%","2%","3%","4%","5%"))

ggsave("Q10onlyBarcode5misextracted.pdf" );








mat = read.table("../data/NsummaryQ20.csv",header =F )

s = subset(mat,!(is.na(mat[,9])) & mat[,16] == 0.1 & mat[,11] > 0 & mat[,2]!= "CUTADAPT" )

s[,2] = sub("CUTADAPTFASTX","CutFastx",s[,2])



x = cbind(s[,11],s[,12],s[,15],s[,16] ,s[,3] + s[,4], 100000 - (s[,3] + s[,4]), s[,4]+s[5],sprintf("N:%2d L:%2d",s[,11],s[,12]  ),10000 * s[,11] + s[,12],s[,2] , s[,3]+s[,4] )

x[,12] = 0;
x[,13] = 0;

colnames(x) = c("NumBar","BarLen","Error Rate","InDelRate","Extracted","Not","MisExtracted","Name","INDEXX","Program","Total","PrintedLabel","RANGE")


for (i in 1:nrow(x)){
	
	
	if(x$Total[i] == 0){
		x$Total[i] = NA;
		x$PrintedLabel[i] =  sprintf("%1.1f%%", 0 )
	}else{
		x$PrintedLabel[i] =  sprintf("%1.1f%%", x$Total[i] /1000)
	}
}



x$RANGE <- cut(x$Total/1000,breaks = c(0,50,80,90.05,91,92,93,94,95,96,97,98,99,101 ),right = FALSE, label = c(1:13) )


p <- ggplot(x, aes(Name, `Error Rate`, label = PrintedLabel )) + geom_tile(aes(fill = RANGE))
p+scale_fill_manual(breaks=c(1:13),values = union (colorRampPalette(c("white","lightgreen"))(3),colorRampPalette(c("white","firebrick3"))(10)) , na.value = "grey50") + theme(legend.position = "none", axis.ticks = element_blank(),strip.text.y = element_text(size = 15), axis.title.x = element_blank(), axis.text.x = element_text(size = 10 , angle = 310, hjust = 0, colour = "black"))+facet_grid( Program ~ . )+ geom_text(aes(size=10))+ scale_y_continuous(breaks=c(0.0,0.01,0.02,0.03,0.04,0.05), labels = expression("0%","1%","2%","3%","4%","5%"))

ggsave("Q20BarcodeandAdapter1totalextracted.pdf" );

x$PrintedLabel =  sprintf("%1.1f%%", x$MisExtracted /1000)


x$RANGE <- cut(x$MisExtracted /(x$Extracted[i] + x$MisExtracted[i]) *100,breaks = c(0.0,0.1,0.2,0.3,0.4,0.5,1,2,5,10,40),right = FALSE, label = c(1,2,3,4,5,6,7,8,9,10) )


for (i in 1:nrow(x)){
	
	x$Total[i] = x$Extracted[i] + x$MisExtracted[i] ;
	
	if(x$Total[i] == 0){
		x$RANGE[i] = NA;
	}
}





p <- ggplot(x, aes(Name, `Error Rate`,  label = PrintedLabel )) + geom_tile(aes(fill = RANGE))
p+scale_fill_manual(breaks=c(1,2,3,4,5,6,7,8,9,10),values =  colorRampPalette(c("white","firebrick3"))(10), na.value = "grey50") +  theme(legend.position = "none", axis.ticks = element_blank(),strip.text.y = element_text(size = 15), axis.title.x = element_blank(), axis.text.x = element_text(size = 10 , angle = 310, hjust = 0, colour = "black"))+facet_grid( Program ~ . )+ geom_text(aes(size=10))+ scale_y_continuous(breaks=c(0.0,0.01,0.02,0.03,0.04,0.05), labels = expression("0%","1%","2%","3%","4%","5%"))


ggsave("Q20BarcodeandAdapter1misextracted.pdf" );


mat = read.table("../data/NsummaryQ10.csv",header =F )

s = subset(mat,!(is.na(mat[,9])) & mat[,16] == 0.1 & mat[,11] > 0 & mat[,2]!= "CUTADAPT" )

s[,2] = sub("CUTADAPTFASTX","CutFastx",s[,2])



x = cbind(s[,11],s[,12],s[,15],s[,16] ,s[,3] + s[,4], 100000 - (s[,3] + s[,4]), s[,4]+s[5],sprintf("N:%2d L:%2d",s[,11],s[,12]  ),10000 * s[,11] + s[,12],s[,2] , s[,3]+s[,4] )

x[,12] = 0;
x[,13] = 0;

colnames(x) = c("NumBar","BarLen","Error Rate","InDelRate","Extracted","Not","MisExtracted","Name","INDEXX","Program","Total","PrintedLabel","RANGE")


for (i in 1:nrow(x)){
	
	
	if(x$Total[i] == 0){
		x$Total[i] = NA;
		x$PrintedLabel[i] =  sprintf("%1.1f%%", 0 )
	}else{
		x$PrintedLabel[i] =  sprintf("%1.1f%%", x$Total[i] /1000)
	}
}



x$RANGE <- cut(x$Total/1000,breaks = c(0,50,80,90.05,91,92,93,94,95,96,97,98,99,101 ),right = FALSE, label = c(1:13) )


p <- ggplot(x, aes(Name, `Error Rate`, label = PrintedLabel )) + geom_tile(aes(fill = RANGE))
p+scale_fill_manual(breaks=c(1:13),values = union (colorRampPalette(c("white","lightgreen"))(3),colorRampPalette(c("white","firebrick3"))(10)) , na.value = "grey50") + theme(legend.position = "none", axis.ticks = element_blank(),strip.text.y = element_text(size = 15), axis.title.x = element_blank(), axis.text.x = element_text(size = 10 , angle = 310, hjust = 0, colour = "black"))+facet_grid( Program ~ . )+ geom_text(aes(size=10))+ scale_y_continuous(breaks=c(0.0,0.01,0.02,0.03,0.04,0.05), labels = expression("0%","1%","2%","3%","4%","5%"))

ggsave("Q10BarcodeandAdapter1totalextracted.pdf" );

x$PrintedLabel =  sprintf("%1.1f%%", x$MisExtracted /1000)


x$RANGE <- cut(x$MisExtracted /(x$Extracted[i] + x$MisExtracted[i]) *100,breaks = c(0.0,0.1,0.2,0.3,0.4,0.5,1,2,5,10,40),right = FALSE, label = c(1,2,3,4,5,6,7,8,9,10) )


for (i in 1:nrow(x)){
	
	x$Total[i] = x$Extracted[i] + x$MisExtracted[i] ;
	
	if(x$Total[i] == 0){
		x$RANGE[i] = NA;
	}
}





p <- ggplot(x, aes(Name, `Error Rate`,  label = PrintedLabel )) + geom_tile(aes(fill = RANGE))
p+scale_fill_manual(breaks=c(1,2,3,4,5,6,7,8,9,10),values =  colorRampPalette(c("white","firebrick3"))(10), na.value = "grey50") +  theme(legend.position = "none", axis.ticks = element_blank(),strip.text.y = element_text(size = 15), axis.title.x = element_blank(), axis.text.x = element_text(size = 10 , angle = 310, hjust = 0, colour = "black"))+facet_grid( Program ~ . )+ geom_text(aes(size=10))+ scale_y_continuous(breaks=c(0.0,0.01,0.02,0.03,0.04,0.05), labels = expression("0%","1%","2%","3%","4%","5%"))


ggsave("Q10BarcodeandAdapter1misextracted.pdf" );















mat = read.table("../data/NsummaryQ20.csv",header =F )

s = subset(mat,!(is.na(mat[,9])) & mat[,16] == 0.5 & mat[,11] > 0 & mat[,2]!= "CUTADAPT" )

s[,2] = sub("CUTADAPTFASTX","CutFastx",s[,2])



x = cbind(s[,11],s[,12],s[,15],s[,16] ,s[,3] + s[,4], 100000 - (s[,3] + s[,4]), s[,4]+s[5],sprintf("N:%2d L:%2d",s[,11],s[,12]  ),10000 * s[,11] + s[,12],s[,2] , s[,3]+s[,4] )

x[,12] = 0;
x[,13] = 0;

colnames(x) = c("NumBar","BarLen","Error Rate","InDelRate","Extracted","Not","MisExtracted","Name","INDEXX","Program","Total","PrintedLabel","RANGE")


for (i in 1:nrow(x)){
	
	
	if(x$Total[i] == 0){
		x$Total[i] = NA;
		x$PrintedLabel[i] =  sprintf("%1.1f%%", 0 )
	}else{
		x$PrintedLabel[i] =  sprintf("%1.1f%%", x$Total[i] /1000)
	}
}



x$RANGE <- cut(x$Total/1000,breaks = c(0,50,80,90.05,91,92,93,94,95,96,97,98,99,101 ),right = FALSE, label = c(1:13) )


p <- ggplot(x, aes(Name, `Error Rate`, label = PrintedLabel )) + geom_tile(aes(fill = RANGE))
p+scale_fill_manual(breaks=c(1:13),values = union (colorRampPalette(c("white","lightgreen"))(3),colorRampPalette(c("white","firebrick3"))(10)) , na.value = "grey50") + theme(legend.position = "none", axis.ticks = element_blank(),strip.text.y = element_text(size = 15), axis.title.x = element_blank(), axis.text.x = element_text(size = 10 , angle = 310, hjust = 0, colour = "black"))+facet_grid( Program ~ . )+ geom_text(aes(size=10))+ scale_y_continuous(breaks=c(0.0,0.01,0.02,0.03,0.04,0.05), labels = expression("0%","1%","2%","3%","4%","5%"))

ggsave("Q20BarcodeandAdapter5totalextracted.pdf" );

x$PrintedLabel =  sprintf("%1.1f%%", x$MisExtracted /1000)


x$RANGE <- cut(x$MisExtracted /(x$Extracted[i] + x$MisExtracted[i]) *100,breaks = c(0.0,0.1,0.2,0.3,0.4,0.5,1,2,5,10,40),right = FALSE, label = c(1,2,3,4,5,6,7,8,9,10) )


for (i in 1:nrow(x)){
	
	x$Total[i] = x$Extracted[i] + x$MisExtracted[i] ;
	
	if(x$Total[i] == 0){
		x$RANGE[i] = NA;
	}
}





p <- ggplot(x, aes(Name, `Error Rate`,  label = PrintedLabel )) + geom_tile(aes(fill = RANGE))
p+scale_fill_manual(breaks=c(1,2,3,4,5,6,7,8,9,10),values =  colorRampPalette(c("white","firebrick3"))(10), na.value = "grey50") +  theme(legend.position = "none", axis.ticks = element_blank(),strip.text.y = element_text(size = 15), axis.title.x = element_blank(), axis.text.x = element_text(size = 10 , angle = 310, hjust = 0, colour = "black"))+facet_grid( Program ~ . )+ geom_text(aes(size=10))+ scale_y_continuous(breaks=c(0.0,0.01,0.02,0.03,0.04,0.05), labels = expression("0%","1%","2%","3%","4%","5%"))


ggsave("Q20BarcodeandAdapter5misextracted.pdf" );


mat = read.table("../data/NsummaryQ10.csv",header =F )

s = subset(mat,!(is.na(mat[,9])) & mat[,16] == 0.5 & mat[,11] > 0 & mat[,2]!= "CUTADAPT" )

s[,2] = sub("CUTADAPTFASTX","CutFastx",s[,2])



x = cbind(s[,11],s[,12],s[,15],s[,16] ,s[,3] + s[,4], 100000 - (s[,3] + s[,4]), s[,4]+s[5],sprintf("N:%2d L:%2d",s[,11],s[,12]  ),10000 * s[,11] + s[,12],s[,2] , s[,3]+s[,4] )

x[,12] = 0;
x[,13] = 0;

colnames(x) = c("NumBar","BarLen","Error Rate","InDelRate","Extracted","Not","MisExtracted","Name","INDEXX","Program","Total","PrintedLabel","RANGE")


for (i in 1:nrow(x)){
	
	
	if(x$Total[i] == 0){
		x$Total[i] = NA;
		x$PrintedLabel[i] =  sprintf("%1.1f%%", 0 )
	}else{
		x$PrintedLabel[i] =  sprintf("%1.1f%%", x$Total[i] /1000)
	}
}



x$RANGE <- cut(x$Total/1000,breaks = c(0,50,80,90.05,91,92,93,94,95,96,97,98,99,101 ),right = FALSE, label = c(1:13) )


p <- ggplot(x, aes(Name, `Error Rate`, label = PrintedLabel )) + geom_tile(aes(fill = RANGE))
p+scale_fill_manual(breaks=c(1:13),values = union (colorRampPalette(c("white","lightgreen"))(3),colorRampPalette(c("white","firebrick3"))(10)) , na.value = "grey50") + theme(legend.position = "none", axis.ticks = element_blank(),strip.text.y = element_text(size = 15), axis.title.x = element_blank(), axis.text.x = element_text(size = 10 , angle = 310, hjust = 0, colour = "black"))+facet_grid( Program ~ . )+ geom_text(aes(size=10))+ scale_y_continuous(breaks=c(0.0,0.01,0.02,0.03,0.04,0.05), labels = expression("0%","1%","2%","3%","4%","5%"))

ggsave("Q10BarcodeandAdapter5totalextracted.pdf" );

x$PrintedLabel =  sprintf("%1.1f%%", x$MisExtracted /1000)


x$RANGE <- cut(x$MisExtracted /(x$Extracted[i] + x$MisExtracted[i]) *100,breaks = c(0.0,0.1,0.2,0.3,0.4,0.5,1,2,5,10,40),right = FALSE, label = c(1,2,3,4,5,6,7,8,9,10) )


for (i in 1:nrow(x)){
	
	x$Total[i] = x$Extracted[i] + x$MisExtracted[i] ;
	
	if(x$Total[i] == 0){
		x$RANGE[i] = NA;
	}
}





p <- ggplot(x, aes(Name, `Error Rate`,  label = PrintedLabel )) + geom_tile(aes(fill = RANGE))
p+scale_fill_manual(breaks=c(1,2,3,4,5,6,7,8,9,10),values =  colorRampPalette(c("white","firebrick3"))(10), na.value = "grey50") +  theme(legend.position = "none", axis.ticks = element_blank(),strip.text.y = element_text(size = 15), axis.title.x = element_blank(), axis.text.x = element_text(size = 10 , angle = 310, hjust = 0, colour = "black"))+facet_grid( Program ~ . )+ geom_text(aes(size=10))+ scale_y_continuous(breaks=c(0.0,0.01,0.02,0.03,0.04,0.05), labels = expression("0%","1%","2%","3%","4%","5%"))


ggsave("Q10BarcodeandAdapter5misextracted.pdf" );



















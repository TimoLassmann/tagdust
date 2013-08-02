 

library(ggplot2)

td10 = read.table ("tagdust_benchmark10.csv", header = F)
td20 = read.table ("tagdust_benchmark20.csv", header = F)

fx1 = read.table("fastx_benchmark_1.csv",header = F)

fx0 = read.table("fastx_benchmark_0.csv",header = F)


t10 = subset(td10,td10[,1] == 8 )
t20 = subset(td20,td20[,1] == 8 )
f1 = subset(fx1,fx1[,1] == 8 )
f0 = subset(fx0,fx0[,1] == 8 )

t10[,8] = "Tagdust10";
t20[,8] = "Tagdust20";
f1[,8] = "FastX1"
f0[,8] = "FastX0"


x =  rbind(t10,t20,f1,f0)

x = x[ order(x[,2]), ]

colnames(x)[8] = "Program"

colnames(x)[2] = "BarcodeLength"
colnames(x)[3] = "ErrorRate"
colnames(x)[5] = "Correctly Extracted"
colnames(x)[6] = "Not Extracted"
colnames(x)[7] = "Incorrectly Assigned"


x[,2] = paste("Barcode Lenght" ,x[,2])

x$BarcodeLength <- factor(x$BarcodeLength, levels = c("Barcode Lenght 4", "Barcode Lenght 6", "Barcode Lenght 8", "Barcode Lenght 10"));


pdf("8Barcodes10Indels.pdf",  height = 10,width = 12,bg="white",pointsize = 36)

p <- ggplot(x, aes((`Correctly Extracted` +`Incorrectly Assigned`) / 10000, `Incorrectly Assigned`/ (`Correctly Extracted` +`Incorrectly Assigned`) * 100, size= ErrorRate , label=Program,colour = Program))
p+geom_point(aes(colour = Program) ) +geom_text(size=0)+  facet_wrap( ~ BarcodeLength,scales = "free") + scale_size(name = "Error Rate",range = c(1, 5))  +scale_y_continuous (name = "Percentage of misassigned reads")+scale_x_continuous  (name = "Percentage of extracted reads",limits=c(75,100))+ labs(title = "8 Barcodes 1% Indels") + geom_vline(xintercept = 90,colour="red", linetype = "longdash") + geom_path(alpha = 0.1);

dev.off();

remove(x)


p <- ggplot(x, aes(x=ErrorRate  , y=`Incorrectly Assigned`/ (`Correctly Extracted` +`Incorrectly Assigned`) * 100  , group=Program,color = Program))
p + geom_line()+  facet_wrap( ~ BarcodeLength,scales = "free")

p <- ggplot(x, aes((`Correctly Extracted` , `Incorrectly Assigned`, size= ErrorRate , label=Program,colour = Program))
p+geom_point(aes(colour = Program) ) +geom_text(size=0)+  facet_wrap( ~ BarcodeLength,scales = "free") + scale_size(name = "Error Rate",range = c(1, 5))  +scale_y_continuous (name = "Percentage of misassigned reads")+scale_x_continuous  (name = "Percentage of extracted reads",limits=c(75,100))+ labs(title = "8 Barcodes 1% Indels") + geom_vline(xintercept = 90,colour="red", linetype = "longdash") + geom_path(alpha = 0.1);



library(gplots)
require(tikzDevice)

td10 = read.table ("tagdust_benchmark10.csv", header = F)

#Mistakes...

x = cbind(
subset(td10, td10[,2] == 4 & td10[,1] == 8)[7],
subset(td10, td10[,2] == 4 & td10[,1] == 24)[7],
subset(td10, td10[,2] == 4 & td10[,1] == 96)[7],
subset(td10, td10[,2] == 6& td10[,1] == 8)[7],
subset(td10, td10[,2] == 6& td10[,1] == 24)[7],
subset(td10, td10[,2] == 6& td10[,1] == 96)[7],
subset(td10, td10[,2] == 8& td10[,1] == 8)[7],
subset(td10, td10[,2] == 8& td10[,1] == 24)[7],
subset(td10, td10[,2] == 8& td10[,1] == 96)[7],
subset(td10, td10[,2] == 10& td10[,1] == 8)[7],
subset(td10, td10[,2] == 10& td10[,1] == 24)[7],
subset(td10, td10[,2] == 10& td10[,1] == 96)[7])
td10 = read.table("fastx_benchmark_1.csv",header = F)

y = cbind(
subset(td10, td10[,2] == 4 & td10[,1] == 8)[7],
subset(td10, td10[,2] == 4 & td10[,1] == 24)[7],
subset(td10, td10[,2] == 4 & td10[,1] == 96)[7],
subset(td10, td10[,2] == 6& td10[,1] == 8)[7],
subset(td10, td10[,2] == 6& td10[,1] == 24)[7],
subset(td10, td10[,2] == 6& td10[,1] == 96)[7],
subset(td10, td10[,2] == 8& td10[,1] == 8)[7],
subset(td10, td10[,2] == 8& td10[,1] == 24)[7],
subset(td10, td10[,2] == 8& td10[,1] == 96)[7],
subset(td10, td10[,2] == 10& td10[,1] == 8)[7],
subset(td10, td10[,2] == 10& td10[,1] == 24)[7],
subset(td10, td10[,2] == 10& td10[,1] == 96)[7])




z = rbind(x,y)

z = as.matrix(z)

heatmap.2(z/10000,dendrogram="none",Rowv=F,Colv = F,trace="none",cellnote=round(z/10000,digits = 1),breaks = c(0,0.5,1,3,7) ,col= c("lightblue","skyblue","steelblue","dodgerblue4"),notecol="black",key = F)



td10 = read.table ("tagdust_benchmark10.csv", header = F)


td10[,8] = td10[,5] + td10[,7];
 #total extracted...

x = cbind(
subset(td10, td10[,2] == 4 & td10[,1] == 8)[8],
subset(td10, td10[,2] == 4 & td10[,1] == 24)[8],
subset(td10, td10[,2] == 4 & td10[,1] == 96)[8],
subset(td10, td10[,2] == 6& td10[,1] == 8)[8],
subset(td10, td10[,2] == 6& td10[,1] == 24)[8],
subset(td10, td10[,2] == 6& td10[,1] == 96)[8],
subset(td10, td10[,2] == 8& td10[,1] == 8)[8],
subset(td10, td10[,2] == 8& td10[,1] == 24)[8],
subset(td10, td10[,2] == 8& td10[,1] == 96)[8],
subset(td10, td10[,2] == 10& td10[,1] == 8)[8],
subset(td10, td10[,2] == 10& td10[,1] == 24)[8],
subset(td10, td10[,2] == 10& td10[,1] == 96)[8])
td10 = read.table("fastx_benchmark_1.csv",header = F)
td10[,8] = td10[,5] + td10[,7];
y = cbind(
subset(td10, td10[,2] == 4 & td10[,1] == 8)[8],
subset(td10, td10[,2] == 4 & td10[,1] == 24)[8],
subset(td10, td10[,2] == 4 & td10[,1] == 96)[8],
subset(td10, td10[,2] == 6& td10[,1] == 8)[8],
subset(td10, td10[,2] == 6& td10[,1] == 24)[8],
subset(td10, td10[,2] == 6& td10[,1] == 96)[8],
subset(td10, td10[,2] == 8& td10[,1] == 8)[8],
subset(td10, td10[,2] == 8& td10[,1] == 24)[8],
subset(td10, td10[,2] == 8& td10[,1] == 96)[8],
subset(td10, td10[,2] == 10& td10[,1] == 8)[8],
subset(td10, td10[,2] == 10& td10[,1] == 24)[8],
subset(td10, td10[,2] == 10& td10[,1] == 96)[8])




z = rbind(x,y)

z = as.matrix(z)

heatmap.2(z/10000,dendrogram="none",Rowv=F,Colv = F,trace="none",cellnote=round(z/10000,digits = 0),breaks = c(70,75,80,85,90) ,col= c("lightblue","skyblue","steelblue","dodgerblue4"),notecol="black",key = F)



rownames(x) = c("2%","3%","4%","5%","6%","7%");
colnames(x) = c("20nt","22nt","24nt","26nt","28nt","30nt","32nt","34nt","36nt","38nt","40nt");






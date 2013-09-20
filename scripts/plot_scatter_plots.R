


files = c("comp_15c_noumi.csv","comp_15c_umi.csv","comp_25c_noumi.csv", "comp_25c_umi.csv", "comp_15c_repeat1_25c_repeat1_noumi.csv","comp_15c_repeat1_25c_repeat1_umi.csv", "comp_15c_repeat2_25c_repeat1_lane2_noumi.csv" , "comp_15c_repeat2_25c_repeat1_lane2_umi.csv");



library(ggplot2)




mat = read.table("15_25_combined_noumi.csv",header = F);

cor = cor(mat[,2],mat[,3]);
frame=ggplot(mat,aes(mat[,2],mat[,3]))


frame + geom_point(colour="black",size=1,shape=20,alpha=0.3) + labs(x="15c ",y="25c",title="15 vs 25 PCR cycles - no UMI") + annotate("text", label =paste("R=", round (cor,6), sep=""), x = min(mat[,2])+10, y = max(mat[,3]) / 10, size = 6, colour = "red") +scale_x_log10() + scale_y_log10();

ggsave("15_25_combined_noumi.pdf")
mat = read.table("15_25_combined_UMI.csv",header = F);

cor = cor(mat[,2],mat[,3]);
frame=ggplot(mat,aes(mat[,2],mat[,3]))


frame + geom_point(colour="black",size=1,shape=20,alpha=0.3) + labs(x="15c ",y="25c",title="15 vs 25 PCR cycles - UMI") + annotate("text", label =paste("R=", round (cor,6), sep=""), x = min(mat[,2])+10, y = max(mat[,3]) / 10, size = 6, colour = "red") +scale_x_log10() + scale_y_log10();

ggsave("15_25_combined_UMI.pdf")





mat = read.table("comp_15c_noumi.csv",header = F);

cor = cor(mat[,2],mat[,3]);
frame=ggplot(mat,aes(mat[,2],mat[,3]))


frame + geom_point(colour="black",size=1,shape=20,alpha=0.3) + labs(x="15c Rep1",y="15c Rep2 ",title="15 PCR cycles - no UMI") + annotate("text", label =paste("R=", round (cor,6), sep=""), x = min(mat[,2])+10, y = max(mat[,3]) / 10, size = 6, colour = "red") +scale_x_log10() + scale_y_log10();

ggsave("comp_15c_noumi.pdf")


mat = read.table("comp_15c_umi.csv",header = F);

cor = cor(mat[,2],mat[,3]);
frame=ggplot(mat,aes(mat[,2],mat[,3]))


frame + geom_point(colour="black",size=1,shape=20,alpha=0.3) + labs(x="15c Rep1",y="15c Rep2 ",title="15 PCR cycles - UMI") + annotate("text", label =paste("R=", round (cor,6), sep=""), x = min(mat[,2])+10, y = max(mat[,3]) / 10, size = 6, colour = "red") +scale_x_log10() + scale_y_log10();

ggsave("comp_15c_umi.pdf")


mat = read.table("comp_25c_noumi.csv",header = F);

cor = cor(mat[,2],mat[,3]);
frame=ggplot(mat,aes(mat[,2],mat[,3]))


frame + geom_point(colour="black",size=1,shape=20,alpha=0.3) + labs(x="25c Rep1",y="25c Rep2 ",title="25 PCR cycles - no UMI") + annotate("text", label =paste("R=", round (cor,6), sep=""), x = min(mat[,2])+10, y = max(mat[,3]) / 10, size = 6, colour = "red") +scale_x_log10() + scale_y_log10();

ggsave("comp_25c_noumi.pdf")


mat = read.table("comp_25c_umi.csv",header = F);

cor = cor(mat[,2],mat[,3]);
frame=ggplot(mat,aes(mat[,2],mat[,3]))

frame + geom_point(colour="black",size=1,shape=20,alpha=0.3) + labs(x="25c Rep1",y="25c Rep2 ",title="25 PCR cycles - UMI") + annotate("text", label =paste("R=", round (cor,6), sep=""), x = min(mat[,2])+10, y = max(mat[,3]) / 10, size = 6, colour = "red") +scale_x_log10() + scale_y_log10();

ggsave("comp_25c_umi.pdf")


mat = read.table("comp_15c_repeat1_25c_repeat1_noumi.csv",header = F);

cor = cor(mat[,2],mat[,3]);
frame=ggplot(mat,aes(mat[,2],mat[,3]))

frame + geom_point(colour="black",size=1,shape=20,alpha=0.3) + labs(x="15c Rep1",y="25c Rep1 ",title="15 vs 25 PCR cycles - noUMI") + annotate("text", label =paste("R=", round (cor,6), sep=""), x = min(mat[,2])+10, y = max(mat[,3]) / 10, size = 6, colour = "red") +scale_x_log10() + scale_y_log10();

ggsave("comp_15c_repeat1_25c_repeat1_noumi.pdf")


mat = read.table("comp_15c_repeat1_25c_repeat1_umi.csv",header = F);

cor = cor(mat[,2],mat[,3]);
frame=ggplot(mat,aes(mat[,2],mat[,3]))

frame + geom_point(colour="black",size=1,shape=20,alpha=0.3) + labs(x="15c Rep1",y="25c Rep1 ",title="15 vs 25 PCR cycles - UMI") + annotate("text", label =paste("R=", round (cor,6), sep=""), x = min(mat[,2])+10, y = max(mat[,3]) / 10, size = 6, colour = "red") +scale_x_log10() + scale_y_log10();

ggsave("comp_15c_repeat1_25c_repeat1_umi.pdf")

mat = read.table("comp_15c_repeat2_25c_repeat1_lane2_noumi.csv",header = F);

cor = cor(mat[,2],mat[,3]);
frame=ggplot(mat,aes(mat[,2],mat[,3]))

frame + geom_point(colour="black",size=1,shape=20,alpha=0.3) + labs(x="15c Rep2",y="25c Rep2 ",title="15 vs 25 PCR cycles - no UMI") + annotate("text", label =paste("R=", round (cor,6), sep=""), x = min(mat[,2])+10, y = max(mat[,3]) / 10, size = 6, colour = "red") +scale_x_log10() + scale_y_log10();

ggsave("comp_15c_repeat2_25c_repeat1_lane2_noumi.pdf")


mat = read.table("comp_15c_repeat2_25c_repeat1_lane2_umi.csv",header = F);

cor = cor(mat[,2],mat[,3]);
frame=ggplot(mat,aes(mat[,2],mat[,3]))

frame + geom_point(colour="black",size=1,shape=20,alpha=0.3) + labs(x="15c Rep2",y="25c Rep2 ",title="15 vs 25 PCR cycles - UMI") + annotate("text", label =paste("R=", round (cor,6), sep=""), x = min(mat[,2])+10, y = max(mat[,3]) / 10, size = 6, colour = "red") +scale_x_log10() + scale_y_log10();

ggsave("comp_15c_repeat2_25c_repeat1_lane2_umi.pdf")


mat = read.table("comp_15c_repeat1_25c_repeat1_lane2_noumi.csv",header = F);

cor = cor(mat[,2],mat[,3]);
frame=ggplot(mat,aes(mat[,2],mat[,3]))

frame + geom_point(colour="black",size=1,shape=20,alpha=0.3) + labs(x="15c Rep1",y="25c Rep2 ",title="15 vs 25 PCR cycles - no UMI") + annotate("text", label =paste("R=", round (cor,6), sep=""), x = min(mat[,2])+10, y = max(mat[,3]) / 10, size = 6, colour = "red") +scale_x_log10() + scale_y_log10();

ggsave("comp_15c_repeat1_25c_repeat1_lane2_noumi.pdf")

mat = read.table("comp_15c_repeat1_25c_repeat1_lane2_umi.csv",header = F);

cor = cor(mat[,2],mat[,3]);
frame=ggplot(mat,aes(mat[,2],mat[,3]))

frame + geom_point(colour="black",size=1,shape=20,alpha=0.3) + labs(x="15c Rep1",y="25c Rep2 ",title="15 vs 25 PCR cycles -  UMI") + annotate("text", label =paste("R=", round (cor,6), sep=""), x = min(mat[,2])+10, y = max(mat[,3]) / 10, size = 6, colour = "red") +scale_x_log10() + scale_y_log10();

ggsave("comp_15c_repeat1_25c_repeat1_lane2_umi.pdf")



mat = read.table("comp_15c_repeat2_25c_repeat1_lane1_noumi.csv",header = F);

cor = cor(mat[,2],mat[,3]);
frame=ggplot(mat,aes(mat[,2],mat[,3]))

frame + geom_point(colour="black",size=1,shape=20,alpha=0.3) + labs(x="15c Rep2",y="25c Rep1 ",title="15 vs 25 PCR cycles -  no UMI") + annotate("text", label =paste("R=", round (cor,6), sep=""), x = min(mat[,2])+10, y = max(mat[,3]) / 10, size = 6, colour = "red") +scale_x_log10() + scale_y_log10();

ggsave("comp_15c_repeat2_25c_repeat1_lane1_noumi.pdf")


mat = read.table("comp_15c_repeat2_25c_repeat1_lane1_umi.csv",header = F);

cor = cor(mat[,2],mat[,3]);
frame=ggplot(mat,aes(mat[,2],mat[,3]))

frame + geom_point(colour="black",size=1,shape=20,alpha=0.3) + labs(x="15c Rep2",y="25c Rep1 ",title="15 vs 25 PCR cycles -  UMI") + annotate("text", label =paste("R=", round (cor,6), sep=""), x = min(mat[,2])+10, y = max(mat[,3]) / 10, size = 6, colour = "red") +scale_x_log10() + scale_y_log10();

ggsave("comp_15c_repeat2_25c_repeat1_lane1_umi.pdf")







mat = read.table("comp_10cells_100cells_noumi.csv",header = F);

cor = cor(mat[,2],mat[,3]);
frame=ggplot(mat,aes(mat[,2],mat[,3]))

frame + geom_point(colour="black",size=1,shape=20,alpha=0.3) + labs(x="10 cells",y="100 cells",title="10 vs 100 cells - no UMI") + annotate("text", label =paste("R=", round (cor,6), sep=""), x = min(mat[,2])+10, y = max(mat[,3]) / 10, size = 6, colour = "red") +scale_x_log10() + scale_y_log10();

ggsave("comp_10cells_100cells_noumi.pdf")

mat = read.table("comp_10cells_1000cells_noumi.csv",header = F);

cor = cor(mat[,2],mat[,3]);
frame=ggplot(mat,aes(mat[,2],mat[,3]))

frame + geom_point(colour="black",size=1,shape=20,alpha=0.3) + labs(x="10 cells",y="1000 cells",title="10 vs 1000 cells - no UMI") + annotate("text", label =paste("R=", round (cor,6), sep=""), x = min(mat[,2])+10, y = max(mat[,3]) / 10, size = 6, colour = "red") +scale_x_log10() + scale_y_log10();

ggsave("comp_10cells_1000cells_noumi.pdf")

mat = read.table("comp_100cells_1000cells_noumi.csv",header = F);

cor = cor(mat[,2],mat[,3]);
frame=ggplot(mat,aes(mat[,2],mat[,3]))

frame + geom_point(colour="black",size=1,shape=20,alpha=0.3) + labs(x="100 cells",y="1000 cells",title="100 vs 1000 cells - no UMI") + annotate("text", label =paste("R=", round (cor,6), sep=""), x = min(mat[,2])+10, y = max(mat[,3]) / 10, size = 6, colour = "red") +scale_x_log10() + scale_y_log10();

ggsave("comp_100cells_1000cells_noumi.pdf")







mat = read.table("comp_10cells_100cells_umi.csv",header = F);

cor = cor(mat[,2],mat[,3]);
frame=ggplot(mat,aes(mat[,2],mat[,3]))

frame + geom_point(colour="black",size=1,shape=20,alpha=0.3) + labs(x="10 cells",y="100 cells",title="10 vs 100 cells - UMI") + annotate("text", label =paste("R=", round (cor,6), sep=""), x = min(mat[,2])+10, y = max(mat[,3]) / 10, size = 6, colour = "red") +scale_x_log10() + scale_y_log10();

ggsave("comp_10cells_100cells_umi.pdf")

mat = read.table("comp_10cells_1000cells_umi.csv",header = F);

cor = cor(mat[,2],mat[,3]);
frame=ggplot(mat,aes(mat[,2],mat[,3]))

frame + geom_point(colour="black",size=1,shape=20,alpha=0.3) + labs(x="10 cells",y="1000 cells",title="10 vs 1000 cells - UMI") + annotate("text", label =paste("R=", round (cor,6), sep=""), x = min(mat[,2])+10, y = max(mat[,3]) / 10, size = 6, colour = "red") +scale_x_log10() + scale_y_log10();

ggsave("comp_10cells_1000cells_umi.pdf")

mat = read.table("comp_100cells_1000cells_umi.csv",header = F);

cor = cor(mat[,2],mat[,3]);
frame=ggplot(mat,aes(mat[,2],mat[,3]))

frame + geom_point(colour="black",size=1,shape=20,alpha=0.3) + labs(x="100 cells",y="1000 cells",title="100 vs 1000 cells - UMI") + annotate("text", label =paste("R=", round (cor,6), sep=""), x = min(mat[,2])+10, y = max(mat[,3]) / 10, size = 6, colour = "red") +scale_x_log10() + scale_y_log10();

ggsave("comp_100cells_1000cells_umi.pdf")
















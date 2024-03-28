####wilcox test####
aa <- data.frame()
for (i in 1:ncol(CK)) {
  bb=data.frame(CK=mean(CK[,i]),Treat=mean(Treat[,i]),p=wilcox.test(CK[,i],Treat[,i],paired =F)[[3]])
  aa <- rbind(aa,bb)
}
row.names(aa) <- species$Group.1
aa <- subset(aa,aa$p<0.05)

###exponential decay regression####
paired_data %>%subset(.,pd!=0) %>% 
  ggplot(aes(x=all_shan,y=pd,color=type))+
  geom_point()+theme_zg()+
geom_smooth(method = "lm",formula = y~exp(-x),size=.5,alpha=.2)+
  scale_color_manual(values = cols)+
  stat_poly_eq(formula = y ~exp(-x),aes(label = paste(..rr.label.., ..p.value.label..,sep = "~`,`~")),size=2.5)+
  labs(y="Relative abundance",x="Bacterial diversity",title = "Total pathogens")
  

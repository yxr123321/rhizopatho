meta.df <- read.csv("data/RR.csv",row.names = 1)
ggplot(meta.df,mapping = aes(x=factor,y=estimate,fill=group))+
  geom_hline(yintercept=0,linetype = "dashed",size=0.2)+
  geom_errorbar(position=position_dodge(-0.8),aes(ymin = ci.lb, ymax = ci.ub), width=0.2,size=0.3)+
geom_point(position=position_dodge(-0.8), size=3, stroke = 0.3,shape=21)+scale_fill_manual(values = CS_cols)+
  geom_text(aes(x = factor, y = ci.ub+0.015, label = size),
            position = position_dodge(width = -0.8),vjust = 0.4, hjust=0.4, size = 2, check_overlap = FALSE)+
  geom_text(aes(x = factor, y = ci.ub+0.04, label = star),
            position = position_dodge(width = -0.8),vjust = 0.4, hjust=0.4, size = 2, check_overlap = FALSE)+
  scale_x_discrete(limits=rev(c("Total","Dryland","Wetland","Gramineae","Leguminosae","Solanaceae")))+
  labs(x = " ", y = "Response ratio",colour = 'black')+
  theme(legend.title = element_blank(),
        legend.position=c(0.2,0.94),
        #legend.direction = "horizontal",
        legend.key = element_rect(fill = "white",size = 2),
        legend.key.width = unit(0.5,"lines"),
        legend.key.height= unit(0.8,"lines"),
        legend.background = element_blank(),
        legend.text=element_text(size=6),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        axis.title=element_text(size=9),
        #axis.title.x = element_blank(),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(colour = 'black', size = 8),
        axis.text.x = element_text(colour = 'black', size = 8),
        axis.line = element_line(colour = 'black',size=0.4),
        axis.line.y = element_blank(),
        axis.ticks = element_line(colour = 'black',size=0.4),
        axis.ticks.y = element_blank())+
  guides(fill = "none")+
  coord_flip()+ggtitle("Bacterial diversity")+theme(plot.title = element_text(hjust = 0.5,size = 9,face = "bold"))#+
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())

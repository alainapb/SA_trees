library(tidyverse)
library(patchwork)

#Rescale functions
mm_rescale <- function(data, comparison){
  (data - min(data))/diff(range(data))*diff(range(comparison)) + min(comparison)
}

mm_rescale_getpoint <- function(data, comparison, point){
  my_perc <- (point - min(comparison))/diff(range(comparison))
  diff(range(data))*my_perc + min(data)
}

widen <- function(data, span){
  mean(range(data)) + c(-0.5, 0.5) *diff(range(data))*span
}

myplot <- function(mydf, xaxis, xlab = "", 
                   ylab = "Standard Effect Size of MPD", ylab2 = "Density",
                   myrange, xbreaks = NULL,
                   mybreaks = -3:0, secbreaks = NULL){
    ggplot()+
    geom_histogram(data = bgr_data2, 
                   mapping=aes(x=get(xaxis), y = ..density..), fill = "grey80")+
    geom_line(mydf, mapping = aes(x= get(xaxis), 
                                  y = mm_rescale(yhat, myrange)), size = 1)+
    geom_hline(yintercept = mm_rescale_getpoint(myrange, 
                                                mydf$yhat, 0),
               linetype="dashed")+
    theme_classic()+
    scale_x_continuous(name = xlab,
                       breaks = if(!is.null(xbreaks)){xbreaks} else{waiver()}) +
    scale_y_continuous(name =  ylab,
                       breaks = mm_rescale_getpoint(myrange, 
                                                    mydf$yhat, mybreaks),
                       labels = mybreaks,
                       sec.axis = sec_axis(trans =~.,
                                           name = ylab2,
                                           breaks = if(!is.null(secbreaks)){secbreaks} else{waiver()},
                                           labels = if(!is.null(secbreaks)){secbreaks} else{waiver()}
                                           )) +
    coord_cartesian(ylim = widen(mm_rescale_getpoint(myrange, mydf$yhat, mybreaks), 1.3)) +
    theme(text = element_text(size=12),
          axis.title = element_text(size =12),
          axis.text.y = element_text(color="black", size = 10),
          axis.text.x = element_text(color="black", size = 10))
}

#load best fitting model data
ls <-readRDS("leafshape_brt.rds")
th <-readRDS("treeheight_brt.rds")
la <-readRDS("lat_brt.rds")
lo <-readRDS("long_brt.rds")
sr <-readRDS("speciesrichness_brt.rds")
np <-readRDS("npp_brt.rds")

#load distribution data
bgr_data2<-readRDS("bgr_data2.rds")

leaf_pdp <- myplot(ls, "log_leafshape", "log(Leaf Shape)", ylab = "", ylab2 = " ", seq(0, 2, 0.2))
height_pdp <- myplot(th, "log_treeheight", "log(Tree Height)", ylab = "", ylab2 = "", seq(.6, 1.1, 0.2))
lat_pdp <- myplot(la, "lat", "Latitude", ylab = "", ylab2 = "", seq(.02, .1, .01))
long_pdp <- myplot(lo, "long", "Longitude", ylab = "", ylab2 = "", seq(0.03 , .05, 0.001))
sr_pdp <- myplot(sr, "ntaxa", "Species Richness", ylab = "", ylab2 = "", 
                 xbreaks = seq(0, 500, 250), seq(-0.001, .007, .001))
npp_pdp <- myplot(np, "NPP", "Net Primary Productivity", ylab = "", ylab2 = "", seq(.0063, .0077, .00001))


#title
title <- ggplot(data.frame(l = "Standard Effect Size of MPD", x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), angle = 90) + 
  theme_void() +
  coord_cartesian(clip = "off")

leaf_pdp$labels$y <- " "

comb<- (leaf_pdp + height_pdp + sr_pdp) / (lat_pdp + long_pdp + npp_pdp )
l_combo<-comb + plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A')
title + l_combo + plot_layout(widths = c(1, 25))

ggsave("figure2.png", height = 6, width = 12)

## 
# 08/24/2018 3D figures
##

# This is to calculate the distance matrix using the package "vegan" and the metaMDS which tests for the best solution with random starts, notice I add the "t" before MRcounts to transpose the matrix
#You can change the distance to a few other options, look into the options (maybe euclidean, raup, manhattan?)
metaMDS_amr_bray <- metaMDS(t(MRcounts(AMR_analytic_data[[2]])),k=3,symmetric=TRUE,plot=TRUE,binary = T,trymax = 1000)
#not recommended at this point, but permutation test 
anosim(scores(metaMDS_amr_bray),pData(AMR_analytic_data[[2]])$Treatment)

# Some ordination plots #########
### This is optional and takes more code, but generally people like using "ggplot2" because there are SOOO many options for editing your figures 

sample.scores <- as.data.frame(scores(metaMDS_amr_bray, "sites"))  #Using the scores function from vegan to extract the sample scores and convert to a data.frame
sample.scores$Group<- pData(AMR_analytic_data[[2]])$Group


#sample.scores$sample <- rownames(sample.scores)  # create a column of sample, from the rownames of sample.scores
#sample.scores$AMU <- pData(amr)$AMU
#sample.scores$DDD <- pData(amr)$sumAll_nAUDD
head(sample.scores)  #look at the data

# 2d plot with first 2 NMDS axes
ggplot() + 
  #geom_text_repel(data=sample.scores,aes(x=NMDS1,y=NMDS2,label=sample),alpha=0.5) +  # add the sample labels
  geom_point(data=sample.scores,aes(x=NMDS1,y=NMDS2),size=2, colour= as.numeric(sample.scores$Group)) +
  coord_equal() +
  theme_bw() +
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=20), # remove x-axis labels
        axis.title.y = element_text(size=20), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())


# A version of a 3D plot
ggplot(sample.scores,aes(x=NMDS1,y=NMDS2,z=NMDS3, color= sample.scores$Group) )+
  theme_void() +
  axes_3D() +
  stat_3D() 

## 3D plotly
p <- plot_ly(sample.scores, x = ~NMDS1, y = ~NMDS2, z = ~NMDS3, color = ~Group, colors = c('#BF382A', '#0C4B8E','black','yellow')) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'NMDS1'),
                      yaxis = list(title = 'NMDS2'),
                      zaxis = list(title = 'NMDS3')))
p
chart_link = api_create(p, filename="NMDS3d-isolate-AMRby-media")
chart_link

### 3D hull plot
# this line of code makes the plot bigger. Otherwise too tiny on screen.
# experiment with the sizing
r3dDefaults$windowRect <- c(0,50, 800, 800)

# to make the interactive plot
hullPlot(sample.scores, "Group")

# to save a still
rgl.snapshot(filename = "3d_Packaging_NMDS.png", fmt = "png")

# to save a gif
# need package "magick" installed
# NB dir is the directory - you'll need to change
# "vids" or create a directory called that in your
# working directory
rot <- spin3d( axis= c( 0, 0, 1 ) )
movie3d( rot, duration= 9, dir = "graphs")


## Biplot using principal components analysis 
isolate.pca <- prcomp(t(MRcounts(AMR_analytic_data[[1]])),
                      center = TRUE) #scale. = TRUE 
plot(isolate.pca, type = "l")
g <- ggbiplot(isolate.pca , obs.scale = 1, var.scale = 1, groups = pData(AMR_analytic_data[[1]])$Group, ellipse = TRUE, circle = TRUE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', legend.position = 'top')
print(g)

## Yes - selective media
Ymedia_isolate.pca <- prcomp(t(MRcounts(Ymedia_amr)),
                             center = TRUE) #scale. = TRUE 
plot(Ymedia_isolate.pca, type = "l")
g <- ggbiplot(Ymedia_isolate.pca , obs.scale = 1, var.scale = 1, groups = pData(Ymedia_amr)$Site, ellipse = TRUE, circle = TRUE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', legend.position = 'top')
print(g)

Nmedia_isolate.pca <- prcomp(t(MRcounts(Nmedia_amr)),
                             center = TRUE) #scale. = TRUE 
plot(Nmedia_isolate.pca, type = "l")
g <- ggbiplot(Nmedia_isolate.pca , obs.scale = 1, var.scale = 1, groups = pData(Nmedia_amr)$Site, ellipse = TRUE, circle = TRUE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', legend.position = 'top')
print(g)

#Heatmap plots #######

# This pdf has examples for choosing color palettes
# https://www.nceas.ucsb.edu/~frazier/RSpatialGuides/colorPaletteCheatsheet.pdf 
#heatmapCols = colorRampPalette(brewer.pal(9,"Reds"),bias=1)(50)
#easy hetmap
#plotMRheatmap(obj = amr,n=20, cexRow = 0.4, cexCol = 0.4,trace = "none", norm=FALSE)

# more complicated heatmap that we can edit
col_distance = dist(t(MRcounts(amr)), method = "euclidean")
col_cluster = hclust(col_distance, method = "ward.D")

heatmap.2(MRcounts(amr),
          main= "Heatmap", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=c("grey","blue"),       # use on color palette defined earlier
          dendrogram="row",     # only draw a row dendrogram
          #Colv=as.dendrogram(col_cluster),
          colRow = as.integer(Pen) )           # turn off column clustering

## 3D PCOA plots, ## This part doesn't work yet
open3d()
#### We plot the data points coming from the ord (the ordination object), with x, y, and z axes in black (ax.col="black)
#### and plot the points using type text (type = "t") and tips that describe the treatment levels (text = c(rep("A",5),rep("B",4)))
ordirgl(metaMDS_amr_jaccard,ax.col="black",type="t",text= c(rep("A",5),rep("B",4)))
#### We try and plot an elipsoide per treatment level (groups=Trt) and using the standard error (kind ="se")
#### with confidence region at 0.9 (conf = 0.9), and we shade (type = "shade") with colors "red" and "blue"
#### and make the shades transparent (not solid) with alpha = 0.25
orglellipse(metaMDS_amr_jaccard,groups=Trt, kind = "se", conf = 0.9,type="shade",col=c("red","blue"),alpha=0.25)
#### We also plot spider or ray lines centered on the centroid of the 3D plot using ordspider()
orglspider(metaMDS_amr_jaccard,groups=Trt,col=c("red","blue"))

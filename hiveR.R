require(HiveR)
require(plyr)
require(colorspace)
require(classInt)

d = ggplot2::diamonds
d = d[,c(1:4,7)]
head(d); dim(d)

# separate carat-size data into equal interval groups
brks = classIntervals(d$carat, n=11, style="quantile")$brks[1:11] # also try 'equal' style
d$carat = findInterval(d$carat, brks)

## NODES DATA

nodegroups = list()
for(i in 1:4){
  vals = as.numeric(unique(d[[i]]))
  nodegroup = data.frame(id = 1:length(vals), lab = unique(d[[i]]), vals = vals,
                         radius = 100 * vals/max(vals), axis = i)
  sizes = table(d[[i]])
  nodegroup$size = as.numeric(sizes[ match(nodegroup$lab, names(sizes)) ])
  nodegroup$size = 2 * nodegroup$size / max(nodegroup$size)
  if(i>1) nodegroup$id = nodegroup$id + max(nodegroups[[i-1]]$id)
  nodegroups[[ names(d)[i] ]] = nodegroup
}
nodegroups

nodes = rbind(nodegroups[[1]], nodegroups[[2]], nodegroups[[3]], nodegroups[[4]])
nodes$lab = as.character(nodes$lab)
nodes$axis = as.integer(nodes$axis)
nodes$radius = as.numeric(nodes$radius)
nodes$color = "#ffffff"
head(nodes)

## EDGES DATA

# first update edge data with new node IDs
head(d)
for(i in 1:4) {
  header = paste0(names(nodegroups)[i], 'id')
  d[[header]] = nodegroups[[i]]$id[ match(as.numeric(d[[i]]), nodegroups[[i]]$vals) ]
}
head(d)

# edges between the 4 axes in terms of node IDs
for(i in 6:8){
  edgegroup = data.frame(id1 = d[[i]], id2 = d[[i+1]], price = d[[5]])
  if(i==6) all_edges = edgegroup else all_edges = rbind(all_edges, edgegroup)
}
head(all_edges); dim(all_edges)

# summarise edge data
edges = aggregate(all_edges$price, by=list(all_edges$id1, all_edges$id2), FUN='mean')
names(edges) = c('id1','id2','price')
edges = edges[with(edges, order(id1,id2)),] # reorder

# set edge weights (stroke thickness)
weights = count(all_edges, vars = c('id1', 'id2')) # summary data
weights = weights[with(weights, order(id1,id2)),] # reorder to match egdes
all(weights$id1 == edges$id1, weights$id2 == edges$id2) # check all IDs match up
edges$weight = weights$freq * 0.004
edges$weight = pmax(edges$weight, 0.2) # set min edge weight to still visible
range(weights$freq)
range(edges$weight)

# normalise prices for each group of edges (to utilise full colour range)
p = edges$price
edges$colorvals = 0

for(i in nodegroups[1:3]){
  sel = edges$id1 %in% range(i$id)[1] : range(i$id)[2]
  edges$colorvals[sel] = (p[sel] - min(p[sel])) / (max(p[sel]) - min(p[sel]))
}

edges$color = paste0(hex(HSV(edges$colorvals * 300, 1, 1)), '60') # set alpha
edges = edges[order(edges$weight, decreasing=T),] # draw thin edges last

head(edges)

hpd = list()
hpd$nodes = nodes
hpd$edges = edges
hpd$type = "2D"
hpd$desc = "Diamonds"
hpd$axis.cols = rep('#00000000', 4) # make invisible
hpd$axLabs = c("carats","cut","colour","clarity")
class(hpd) = "HivePlotData"

# Check data correctly formatted
chkHPD(hpd, confirm = TRUE)

plotHive(hpd, axLabs = hpd$axLabs, ch = 0.1)

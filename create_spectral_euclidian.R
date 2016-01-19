library(plyr)
library(RPostgreSQL)
library(RJSONIO)
library(igraph)
library(tm)
library(ggplot2)
library(timeSeries)
library(wavelets)
library(edf)

drv <- dbDriver("PostgreSQL")

foo2 <- read.edf('file.edf')

tmp <- c()

for(i in 1:length(foo2$signal)){
  tmp <- cbind(tmp,foo2$signal[[i]]$data)
}

events<-list()
specs<-list()
events<-list()
specs<-list()

slices <- length(foo2$signal[[1]]$data)/foo2$header.signal[[1]]$samplingrate

for(i in 1:slices){
  events[[paste(i,sep='')]] <- c()
  for(j in 1:length(foo2$signal)){
    ds1 <- tmp[((i*foo2$header.signal[[1]]$samplingrate)-(foo2$header.signal[[1]]$samplingrate-1)):(i*foo2$header.signal[[1]]$samplingrate),j]
    specs[[j]]<-spectrum(ds1, plot=FALSE, method = "ar")
    events[[paste(i,sep='')]]<-cbind(events[[paste(i,sep='')]],specs[[j]]$spec)
  }
}

nodes <- c()
links <- c()
rootsl <- list()

for(j in length(events):1){
  matched<-list()
  interact<-list()
  wt.t1 <- list()
  for(k in length(events):1){
    for(l in 1:length(foo2$signal)){
      if(k<j && is.null(matched[[paste(j,l,sep='_')]])){
        eucl<-dist(rbind(data.frame(events[j])[,l],data.frame(events[k])[,l]), method = "euclidian")
        
        if(eucl<=200){
          #correl<-cor(data.frame(events[j])[,l],data.frame(events[k])[,l])
          print(eucl)
          #bar <- data.frame(data.frame(events[j])[,l]-data.frame(events[k])[,l])
          links<-rbind(links,c(k,j,paste(l,'_',eucl,sep='')))
          #links<-rbind(links,c(k,j,paste(l,'_',paste(which(bar[,1]==0),collapse='-'),'_',paste(data.frame(events[j])[which(bar[,1]==0),l],collapse='-'),sep='')))
          rootsl[[paste(l,'_',eucl,sep='')]] <- c(k,paste(l,'_',eucl,sep=''))
          if(is.null(rootsl[[paste(l,sep='')]]) || rootsl[[paste(l,sep='')]][1]>k) rootsl[[paste(l,sep='')]] <- c(k,paste(l,'_',eucl,sep=''))
          matched[[paste(j,l,sep='_')]]<-1
          if(is.null(interact[[paste(l,sep='')]])) interact[[paste(l,sep='')]]<-l
        }
      }
    }
  }
  sdate <- foo2$header.global$startdate
  stime <- foo2$header.global$starttime
  ndate <- strptime(paste(gsub('\\.','-',sdate), gsub('\\.',':',stime), sep=" "), "%m-%d-%y %H:%M:%S")
  if(ndate$sec==59){
    ndate$sec <- 00
    if(ndate$min==59){
      ndate$min <- 00
      ndate$hour <- ndate$hour+1
    } else{
      ndate$min <- ndate$min+1
    }
  } else{
    ndate$sec <- ndate$sec+1
  }
  nodes <- rbind(nodes, c(as.numeric(j),paste(interact,collapse=', '),gsub(" BST","",ndate)))
}

roots<-NULL
roots<-do.call(rbind.data.frame, rootsl)

colnames(nodes) <- c('node_id','tags','dpub')
colnames(links) <- c('source','target','tag')
colnames(roots) <- c('root_node_id','tag')
nodes<-transform(nodes, node_id = as.numeric(node_id))
nodes <- nodes[order(nodes[,1]),]


conn1 <- dbConnect(drv, host="localhost", port=5432, dbname="cascades_eeg_large_wavelet_dissimilarity", user="postgres", password="postgres")

for(l in 1:nrow(nodes)){
  #print(paste('{"name":"',nodes[l,1],'","tags":"',gsub('-','',gsub('_','',ds2[which(ds2==nodes[l,1]),2])),'"},',sep=''))
  dbSendQuery(conn1, paste("INSERT INTO nodes (uri, dpub, interactions) VALUES('", nodes[l,1], "', to_timestamp('", nodes[l,3] ,"','YYYY-MM-DD HH24:MI:SS'),'", nodes[l,2] ,"')",sep=""))
}
dbDisconnect(conn1)
conn1 <- dbConnect(drv, host="localhost", port=5432, dbname="cascades_eeg_large_wavelet_dissimilarity", user="postgres", password="postgres")
for(k in 1:nrow(links)){
  dbSendQuery(conn1, paste("INSERT INTO links (source_node_uri, target_node_uri, interaction) VALUES('", links[k,1], "','", links[k,2] ,"','", links[k,3] ,"')",sep=""))
}

dbDisconnect(conn1)
conn1 <- dbConnect(drv, host="localhost", port=5432, dbname="cascades_eeg_large_wavelet_dissimilarity", user="postgres", password="postgres")
for(m in 1:nrow(roots)){
  dbSendQuery(conn1, paste("INSERT INTO roots (root_node_uri, interaction) VALUES('", roots[m,1], "','", roots[m,2] ,"')",sep=""))
}

edges <- links[,c(1,2)]
colnames(edges) <- c("id1","id2")

#generate the full graph
g <- graph.data.frame(edges,directed=FALSE)

clus <- clusters(g)$membership
cascades <- list()

for(z in 1:length(clus)){
  cascades[[clus[z]]] <- c(cascades[clus[z]][[1]],get.vertex.attribute(g,"name")[z])
}
dbDisconnect(conn1)
conn1 <- dbConnect(drv, host="localhost", port=5432, dbname="cascades_eeg_large_wavelet_dissimilarity", user="postgres", password="postgres")

for(w in 1:length(cascades)){
  dbSendQuery(conn1, paste("INSERT INTO cascades (path_cnt) VALUES(",length(cascades[[w]]),")",sep=""))
  nextcas <- dbGetQuery(conn1,"SELECT id FROM cascades ORDER BY id DESC LIMIT 1;")
  for(v in 1:length(cascades[[w]])){
    dbSendQuery(conn1, paste("INSERT INTO nodes_cascades (node_uri, cascades_id) VALUES('", cascades[[w]][v] , "','",nextcas[1,1],"')",sep=""))
  }
}
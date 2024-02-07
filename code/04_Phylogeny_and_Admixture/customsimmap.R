customsimmap<-function (tree, colors = NULL, fsize = 1, ftype = "reg", lwd = 2, 
          pts = FALSE, node.numbers = FALSE, mar = NULL, add = FALSE, 
          offset = NULL, direction = "rightwards", type = "phylogram", 
          setEnv = TRUE, part = if (type == "arc") 0.5 else 1, xlim = NULL, 
          ylim = NULL, nodes = "intermediate", tips = NULL, maxY = NULL, 
          hold = TRUE, split.vertical = FALSE, lend = 2, asp = NA, 
          outline = FALSE, plot = TRUE, underscore = FALSE, arc_height = 2) 
{
  if (inherits(tree, "multiPhylo")) {
    par(ask = TRUE)
    for (i in 1:length(tree)) plotSimmap(tree[[i]], colors = colors, 
                                         fsize = fsize, ftype = ftype, lwd = lwd, pts = pts, 
                                         node.numbers = node.numbers, mar, add, offset, direction, 
                                         type, setEnv, part, xlim, ylim, nodes, tips, maxY, 
                                         hold, split.vertical, lend, asp, outline, plot, underscore)
  }
  else {
    if (!inherits(tree, "phylo")) 
      stop("tree should be object of class \"phylo\"")
    if (is.null(tree$maps)) 
      stop("tree should contain mapped states on edges.")
    ftype <- which(c("off", "reg", "b", "i", "bi") == ftype) - 
      1
    if (!ftype) 
      fsize = 0
    if (is.null(colors)) {
      st <- sort(unique(unlist(sapply(tree$maps, names))))
      colors <- palette()[1:length(st)]
      names(colors) <- st
      if (length(st) > 1) {
        cat("no colors provided. using the following legend:\n")
        print(colors)
      }
    }
    if (!underscore) 
      tree$tip.label <- gsub("_", " ", tree$tip.label)
    if (is.null(mar)) 
      mar = rep(0.1, 4)
    if (hold) 
      null <- dev.hold()
    if (type == "phylogram") {
      if (direction %in% c("upwards", "downwards")) {
        if (outline) {
          fg <- par()$fg
          par(fg = "transparent")
          black <- colors
          black[] <- fg
          updownPhylogram(tree, colors = black, fsize, 
                          ftype, lwd = lwd + 2, pts, node.numbers, 
                          mar, add, offset, direction, setEnv, xlim, 
                          ylim, nodes, tips, split.vertical, lend, 
                          asp, plot, underscore)
          par(fg = fg)
        }
        updownPhylogram(tree, colors, fsize, ftype, lwd, 
                        pts, node.numbers, mar, add = if (outline) 
                          TRUE
                        else add, offset, direction, setEnv, xlim, 
                        ylim, nodes, tips, split.vertical, lend, asp, 
                        plot, underscore)
      }
      else {
        if (outline) {
          fg <- par()$fg
          par(fg = "transparent")
          black <- colors
          black[] <- fg
          plotPhylogram(tree, colors = black, fsize, 
                        ftype, lwd = lwd + 2, pts, node.numbers, 
                        mar, add, offset, direction, setEnv, xlim, 
                        ylim, nodes, tips, split.vertical, lend, 
                        asp, plot, underscore)
          par(fg = fg)
        }
        plotPhylogram(tree, colors, fsize, ftype, lwd, 
                      pts, node.numbers, mar, add = if (outline) 
                        TRUE
                      else add, offset, direction, setEnv, xlim, 
                      ylim, nodes, tips, split.vertical, lend, asp, 
                      plot, underscore)
      }
    }
    else if (type == "fan") {
      if (outline) {
        fg <- par()$fg
        par(fg = "transparent")
        black <- colors
        black[] <- fg
        plotFan(tree, colors = black, fsize, ftype, lwd = lwd + 
                  0.75, mar, add, part, setEnv, xlim, ylim, tips, 
                maxY, lend, plot, offset)
        par(fg = fg)
      }
      plotFan(tree, colors, fsize, ftype, lwd, mar, add = if (outline) 
        TRUE
        else add, part, setEnv, xlim, ylim, tips, maxY, lend, 
        plot, offset)
    }
    else if (type == "arc") {
      if (outline) {
        fg <- par()$fg
        par(fg = "transparent")
        black <- colors
        black[] <- fg
        arcPhylogram(tree, colors = black, fsize, ftype, 
                     lwd = lwd + 1, mar, add, part, setEnv, xlim, 
                     ylim, tips, maxY, lend, plot, offset, arc_height)
        par(fg = fg)
      }
      arcPhylogram(tree, colors, fsize, ftype, lwd, mar, 
                   add = if (outline) 
                     TRUE
                   else add, part, setEnv, xlim, ylim, tips, maxY, 
                   lend, plot, offset, arc_height)
    }
    else if (type == "cladogram") {
      if (outline) {
        fg <- par()$fg
        par(fg = "transparent")
        black <- colors
        black[] <- fg
        plotCladogram(tree, colors = black, fsize, ftype, 
                      lwd = lwd + 2, mar, add, offset, direction, 
                      xlim, ylim, nodes, tips, lend, asp, plot)
        par(fg = fg)
      }
      plotCladogram(tree, colors, fsize, ftype, lwd, mar, 
                    add = if (outline) 
                      TRUE
                    else add, offset, direction, xlim, ylim, nodes, 
                    tips, lend, asp, plot)
    }
    if (hold) 
      null <- dev.flush()
  }
}

plotFan<-function(tree,colors,fsize,ftype,lwd,mar,add,part,setEnv,xlim,ylim,tips,maxY,lend,plot,offset){
  if(!plot) cat("plot=FALSE option is not permitted for type=\"fan\". Tree will be plotted.\n")
  if(is.null(offset)) offset<-1
  # reorder
  cw<-reorder(tree)
  pw<-reorder(tree,"pruningwise")
  # count nodes and tips
  n<-Ntip(cw)
  m<-cw$Nnode 
  # get Y coordinates on uncurved space
  Y<-vector(length=m+n)
  if(is.null(tips)) tips<-1:n
  if(part<1.0) Y[cw$edge[cw$edge[,2]<=n,2]]<-0:(n-1)
  else Y[cw$edge[cw$edge[,2]<=n,2]]<-tips
  nodes<-unique(pw$edge[,1])
  for(i in 1:m){
    desc<-pw$edge[which(pw$edge[,1]==nodes[i]),2]
    Y[nodes[i]]<-(min(Y[desc])+max(Y[desc]))/2
  }
  if(is.null(maxY)) maxY<-max(Y)
  Y<-setNames(Y/maxY*2*pi,1:(n+m))
  Y<-part*cbind(Y[as.character(cw$edge[,2])],Y[as.character(cw$edge[,2])])
  R<-nodeHeights(cw)
  # now put into a circular coordinate system
  x<-R*cos(Y)
  y<-R*sin(Y)
  # optimize x & y limits
  par(mar=mar)
  offsetFudge<-1.37 # empirically determined
  OFFSET<-0
  pp<-par("pin")[1]
  sw<-fsize*(max(strwidth(cw$tip.label,units="inches")))+
    offsetFudge*OFFSET*fsize*strwidth("W",units="inches") 
  alp<-optimize(function(a,H,sw,pp) (2*a*1.04*max(H)+2*sw-pp)^2,H=R,sw=sw,pp=pp,
                interval=c(0,1e6))$minimum
  if(part<=0.25) x.lim<-y.lim<-c(0,max(R)+sw/alp)
  else if(part>0.25&&part<=0.5){ 
    x.lim<-c(-max(R)-sw/alp,max(R)+sw/alp)
    y.lim<-c(0,max(R)+sw/alp)
  } else x.lim<-y.lim<-c(-max(R)-sw/alp,max(R)+sw/alp)
  if(is.null(xlim)) xlim<-x.lim
  if(is.null(ylim)) ylim<-y.lim
  # plot tree
  if(!add) plot.new()
  plot.window(xlim=xlim,ylim=ylim,asp=1)
  # plot radial lines (edges)
  ## first, the lines emerging from the root (if there are only two):
  jj<-which(cw$edge[,1]==(Ntip(cw)+1))
  if(length(jj)==2){
    m.left<-cumsum(cw$maps[[jj[1]]])/sum(cw$maps[[jj[1]]])
    xx.left<-c(x[jj[1],1],x[jj[1],1]+(x[jj[1],2]-x[jj[1],1])*m.left)
    yy.left<-c(y[jj[1],1],y[jj[1],1]+(y[jj[1],2]-y[jj[1],1])*m.left)
    m.right<-cumsum(cw$maps[[jj[2]]])/sum(cw$maps[[jj[2]]])
    xx.right<-c(x[jj[2],1],x[jj[2],1]+(x[jj[2],2]-x[jj[2],1])*m.right)
    yy.right<-c(y[jj[2],1],y[jj[2],1]+(y[jj[2],2]-y[jj[2],1])*m.right)
    xx<-c(xx.left[length(xx.left):1],xx.right[2:length(xx.right)])
    yy<-c(yy.left[length(yy.left):1],yy.right[2:length(yy.right)])
    col<-colors[c(names(m.left)[length(m.left):1],names(m.right))]
    segments(xx[2:length(xx)-1],yy[2:length(yy)-1],xx[2:length(xx)],yy[2:length(yy)],
             col=col,lwd=lwd,lend=lend)
  } else jj<-NULL
  for(i in 1:nrow(cw$edge)){
    if(i%in%jj==FALSE){
      maps<-cumsum(cw$maps[[i]])/sum(cw$maps[[i]])
      xx<-c(x[i,1],x[i,1]+(x[i,2]-x[i,1])*maps)
      yy<-c(y[i,1],y[i,1]+(y[i,2]-y[i,1])*maps)
      for(i in 1:(length(xx)-1)) lines(xx[i+0:1],yy[i+0:1],col=colors[names(maps)[i]],
                                       lwd=lwd,lend=lend)
    }
  }
  # plot circular lines
  for(i in 1:m+n){
    r<-R[match(i,cw$edge)]
    a1<-min(Y[which(cw$edge==i)])
    a2<-max(Y[which(cw$edge==i)])
    draw.arc(0,0,r,a1,a2,lwd=lwd,col=colors[names(cw$maps[[match(i,cw$edge[,1])]])[1]])
  }
  # plot labels
  for(i in 1:n){
    ii<-which(cw$edge[,2]==i)
    aa<-Y[ii,2]/(2*pi)*360
    adj<-if(aa>90&&aa<270) c(1,0.25) else c(0,0.25)
    tt<-if(aa>90&&aa<270) paste(cw$tip.label[i],paste(rep(" ",offset),
                                                      collapse=""),sep="") else paste(paste(rep(" ",offset),collapse=""),
                                                                                      cw$tip.label[i],sep="")
    aa<-if(aa>90&&aa<270) 180+aa else aa
    if(ftype) text(x[ii,2],y[ii,2],tt,srt=aa,adj=adj,cex=fsize,font=ftype)
  }
  if(setEnv){
    PP<-list(type="fan",use.edge.length=TRUE,node.pos=1,
             show.tip.label=if(ftype) TRUE else FALSE,show.node.label=FALSE,
             font=ftype,cex=fsize,adj=0,srt=0,no.margin=FALSE,label.offset=offset,
             x.lim=xlim,y.lim=ylim,direction="rightwards",tip.color="black",
             Ntip=Ntip(cw),Nnode=cw$Nnode,edge=tree$edge,
             xx=c(x[sapply(1:n,function(x,y) which(x==y)[1],y=cw$edge[,2]),2],x[1,1],
                  if(m>1) x[sapply(2:m+n,function(x,y) which(x==y)[1],y=cw$edge[,2]),2] else c()),
             yy=c(y[sapply(1:n,function(x,y) which(x==y)[1],y=cw$edge[,2]),2],y[1,1],
                  if(m>1) y[sapply(2:m+n,function(x,y) which(x==y)[1],y=cw$edge[,2]),2] else c()))
    assign("last_plot.phylo",PP,envir=.PlotPhyloEnv)
  }
}
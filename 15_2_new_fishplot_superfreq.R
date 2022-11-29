
recttext <- function(xl, yb, xr, yt, text, rectArgs = NULL, textArgs = NULL) {
  center <- c(mean(c(xl, xr)), mean(c(yb, yt)))
  do.call('rect', c(list(xleft = xl, ybottom = yb, xright = xr, ytop = yt), rectArgs))
  do.call('text', c(list(x = xl*1.01, y = center[2], labels = text, adj=c(0, 0.5)), textArgs))
}

library(superFreq)
library(ggthemes)
library(dplyr)


# change a bit some superFreq functions to plot
addSubclone = function(cT, stories, ylims, dodgyness, colourPool=c(), margin=0.02, preNorm=1) {
  if ( length(colourPool) == 0 ) colourPool = mcri(c('black', 'blue', 'red', 'green', 'orange', 'magenta', 'cyan', 'violet', 'lightblue', 'grey', 'darkblue'))
  subClones = names(cT)
  subStories = stories[names(cT),,drop=F]
  maxSize = ylims[2,]-ylims[1,]
  cloneSum = margin*maxSize+colsums(t(margin*maxSize + t(subStories)))
  cloneSum[cloneSum == 0] = 1   #this happens if all subclones are 0 at a sample. this avoids NaNs.
  norm = pmin(1, maxSize/cloneSum)
  base = ylims[1,]
  usedCols = c()
  x = c()
  y = c()
  for ( subClone in subClones ) {
    subStory = subStories[subClone,]*norm
    range = rbind(base+margin*maxSize*norm, base + margin*maxSize*norm + subStory)

    range[2,] = pmax(base+margin*maxSize*norm, range[2,])
    subrange = rbind(range[1,], range[2]-margin*maxSize*norm)
    pos = addStream(range, col=colourPool[1], dodgyness=dodgyness[subClone])
    print(pos)
    x = c(x, pos$x_label)
    y = c(y, pos$y_label)
    usedCols = c(usedCols, colourPool[1])
    names(usedCols)[length(usedCols)] = subClone
    colourPool = colourPool[-1]
    if ( length(colourPool) == 0 ) colourPool = mcri(c('black', 'blue', 'red', 'green', 'orange', 'magenta', 'cyan', 'violet', 'lightblue', 'grey', 'darkblue'))
    if ( length(cT[[subClone]]) > 0 ) {
      out = addSubclone(cT[[subClone]], stories, range, dodgyness, colourPool, margin, preNorm=norm)
      usedCols = c(usedCols, out$usedCols)
      colourPool = out$colourPool
      x = c(x, out$x)
      y = c(y, out$y)
      print(out$usedCols)
      print(x)
    }
    base = range[2,]
  }

  return(list('colourPool'=colourPool, 'usedCols'=usedCols, 'x'=x, 'y'=y))
}

addStream = function(ylims, col='grey', dodgyness=0) {
  x_label = 0
  y_label = 0
  for ( sample in 2:ncol(ylims) ) {
    x1 = sample-1
    x2 = sample
    y1h = ylims[2, sample-1]
    y1l = ylims[1, sample-1]
    y2h = ylims[2, sample]
    y2l = ylims[1, sample]
    range = c(0,1)
    if ( y1h == y1l & y2h > y2l & !any(ylims[2,]-ylims[1,] != 0 & 1:ncol(ylims) < sample) )  {
      range[1] = 1-sqrt(y2h-y2l)
    }
    if ( y2h == y2l & y1h > y1l & !any(ylims[2,]-ylims[1,] != 0 & 1:ncol(ylims) > sample)) {
      range[2] = sqrt(y1h-y1l)
    }
    pos = addStreamSegment(x1, x2, y1l, y1h, y2l, y2h, range=range, col=col,dodgyness=dodgyness)
    if (x_label == 0) {
      if (max(pos$y) - min(pos$y) > 0.1 ) {
        x_label = min(pos$x)
        y_label = min(pos$y)
      }
    }
  }
  return(list('x_label'=x_label, 'y_label'=y_label))
}

#helper third degree polynomial
third = function(x, x0, a, b, c) a + b*(x-x0) + c*(x-x0)^3

#plots a smooth stream segment.
addStreamSegment = function(x1, x2, y1low, y1high, y2low, y2high, range=c(0,1), col, parts = 100, dodgyness=0) {
  sh = (y1high-y2high)/2
  sl = (y1low-y2low)/2
  z = (x1-x2)/2

  ah = (y1high+y2high)/2
  bh = 3*sh/(2*z)
  ch = -sh/(2*z*z*z)
  al = (y1low+y2low)/2
  bl = 3*sl/(2*z)
  cl = -sl/(2*z*z*z)
  x0 = (x1+x2)/2

  x = seq(from=x1, to=x2, length.out=parts)
  yhigh = third(x, x0, ah, bh, ch)
  ylow = third(x, x0, al, bl, cl)

  if ( range[1] == 0 & range[2] < 1 ) {
    xnorm = (x - x1)/(x2-x1)
    r = range[2]
    yShift = (yhigh-ylow)*(1.5/r^2*xnorm^2 - xnorm^3/r^3)^3*4
    yShift = ifelse(xnorm < r, yShift, 0)
    ylow = ylow + yShift
    yhigh = yhigh - yShift
    x = x[xnorm < r]
    ylow = ylow[xnorm < r]
    yhigh = yhigh[xnorm < r]
  }
  if ( range[1] > 0 & range[2] == 1 ) {
    xnorm = (x2 - x)/(x2-x1)
    r = (1-range[1])
    yShift = (yhigh-ylow)*(1.5/r^2*xnorm^2 - xnorm^3/r^3)^3*4
    yShift = ifelse(xnorm < r, yShift, 0)
    ylow = ylow + yShift
    yhigh = yhigh - yShift
    x = x[xnorm < r]
    ylow = ylow[xnorm < r]
    yhigh = yhigh[xnorm < r]
  }
  polygon(c(x, rev(x)), c(yhigh, rev(pmin(ylow, yhigh))), col=col, border=NA)
  return(list('x'=c(x, rev(x)), 'y'=c(yhigh, rev(pmin(ylow, yhigh)))))
}



plotRiver = function(cloneTree, cloneStories, storyList, allStories, variants, cex=1, genome='hg19', normalise=T, xlim='default', ylim='default', labels=T, setPar=T, sampleOrder='default', excludeClones=c(), markDodgy=T, colourPool = c(), ignoreStoriesBelowSigma=2, smallPlot=F, annotationMethod='VariantAnnotation', plotSampleNames=T, ...) {
  variants$variants = variants$variants[colnames(cloneStories$stories)]
  
  if ( length(colourPool) == 0 )
    colourPool = mcri(c('black', 'blue', 'red', 'green', 'orange', 'magenta', 'cyan', 'violet', 'lightblue', 'grey', 'darkblue'))

  if ( length(sampleOrder) == 1 && sampleOrder == 'default' ) sampleOrder=colnames(cloneStories$stories)
  cloneStories$stories = cloneStories$stories[,sampleOrder,drop=F]
  cloneStories$errors = cloneStories$errors[,sampleOrder,drop=F]

  if ( length(excludeClones) > 0 ) {
    cloneTree = superFreq:::excludeClonesFromTree(cloneTree, excludeClones)
    cloneStories = cloneStories[!(rownames(cloneStories) %in% excludeClones), ]
    storyList = storyList[!(names(storyList) %in% excludeClones)]
  }
  
  dodgyness = superFreq:::getDodgyness(storyList, cloneStories)
  if ( !markDodgy ) dodgyness = 0*dodgyness

  maxLength = 20
  if ( smallPlot ) maxLength = 17
  cloneLabels = lapply(storyList, function(rows) superFreq:::storyToLabel(allStories[rows,], variants, genome=genome,maxLength=maxLength,  mergeCNAs=T, annotationMethod=annotationMethod))
  names(cloneLabels) = names(storyList)
  stories = abs(cloneStories$stories)
  rownames(stories) = rownames(cloneStories)
  purity = sapply(1:ncol(cloneStories$stories), function(i) max(abs(stories[names(cloneTree),i])))
  leadingNormal = F
  if ( any(purity==0) ) {
    leadingNormal = T
    purity[purity==0] = 1
  }
  stories[stories < cloneStories$errors*ignoreStoriesBelowSigma & stories < 0.5] = 0
  if ( normalise ) stories = t(t(stories)/purity)
  if ( !leadingNormal ) {
    stories = cbind(rep(0, nrow(stories)), stories)
    colnames(stories)[1] = 'germline'
    if ( any(rownames(stories)=='germline') ) stories['germline','germline'] = 1
  }
  x = 1:ncol(stories)

  if ( setPar ) {
    par(oma=rep(0, 4))
    par(mar=c(2, 5, 1, 0))
  }
  if ( !labels & xlim[1] == 'default' ) xlim = c(1, max(x))
  if ( labels & xlim[1] == 'default' ) xlim = c(1, max(x)+ceiling(nrow(cloneStories)/2)*1.5)
  if ( smallPlot ) xlim[2] = xlim[2]*1.05
  if ( ylim[1] == 'default' ) ylim = c(-0.02,1)

  plot(1, type='n', xlim=xlim, ylim=ylim, xaxt='n', frame.plot=F,
       ylab='clonality', xlab = '', cex.axis=cex, cex.lab=cex)
  cloneCols = addSubclone(cloneTree, stories, ylims = matrix(rep(c(0,1), ncol(stories)), nrow=2), dodgyness, colourPool=colourPool, margin=0.02)
  prev_y = 0.4
  prev_x = 0
  for (i in nrow(cloneStories):1) {
    clone = rownames(cloneStories)[i]
    xText = max(x) + ceiling(i/2)*1.5 - 0.9

    if ( smallPlot ) xText = max(x) + ceiling(i/2)*1.5 - 1.2
    y0 = i/2 - floor(i/2)

      clone = names(cloneCols$usedCols)[i]
      retained_labels = cloneLabels[[clone]]
      retained_labels = retained_labels[retained_labels$font==2,]
      if (dim(retained_labels)[1] > 0) {
        font = retained_labels$font
        col = cloneCols$usedCols[i]
        uu = stories[clone,] >0
        names(uu) = 1:(length(uu))
        min_x = as.numeric(names(which(uu>0))[1])
        xlab = cloneCols$x[i] + 0.03
        ylab = cloneCols$y[i] + 0.02
        if (abs(xlab - prev_x) >0.3){
          prev_y = 0.4
        }
        prev_x = xlab
        #text(min_x - 0.5 + i * 0.02, 0.7 + i*0.05, retained_labels$label, col=col, adj=0, cex=0.9, font=font)
        recttext(xlab , prev_y, xlab  + 0.037*(length(sampleOrder)**1.9) , (prev_y +   cex/27 * length(retained_labels$label)), paste(sapply(strsplit(as.character(retained_labels$label), " "), `[`, 1), collapse='\n'),
         rectArgs = list(col = col, border = NA),
         textArgs = list(col = 'white', cex = cex))
        segments(xlab+0.002, ylab+.03, x1 = xlab+0.002, y1 = prev_y,
         col = col, lty=1, lwd=3)

        prev_y = prev_y + cex/27 * length(retained_labels$label) + 0.01
    }
  }

  segments(1:ncol(stories), 0.02, 1:ncol(stories), ylim[2], lwd=5, col=rgb(0.7, 0.7, 0.7, 0.3))
  segments(1:ncol(stories), 0.02, 1:ncol(stories), ylim[2], lwd=2, col=rgb(0.3, 0.3, 0.3, 0.3))
  if ( plotSampleNames ) text(1:ncol(stories), -0.02, gsub('germline', 'GL_6', colnames(stories)), srt=20, cex=cex)

  if ( setPar ) {
    par(oma=rep(0, 4))
    par(mar=rep(4, 4))
  }

  return(cloneCols$usedCols)
}


# script
genome='hg19'
annotationMethod='VariantAnnotation'
noneg= superFreq:::noneg
colsums=superFreq:::colsums



for (patient in 1:6) {
  left_color = c("#FFFFFF", pull(ggthemes_data$tableau[['color-palettes']]$regular$Summer)[2:8])
  right_color = c("#FFFFFF", pull(ggthemes_data$tableau[['color-palettes']]$regular[['Nuriel Stone']]))
  if (patient==3) {
    right_color = right_color[-4]
  }
  for (side in c('L', 'R')) {
    load(paste("results/superFreq/R/P", patient, "_", side, "/stories.Rdata", sep=''))
    variants = stories$variants
    normalVariants = stories$normalVariants
    stories = stories$stories
    rm(normalVariants)
    ts = names(stories)
    cloneTree=stories[[ts]]$consistentClusters$cloneTree
    cloneStories=stories[[ts]]$consistentClusters$cloneStories
    storyList=stories[[ts]]$consistentClusters$storyList
    allStories=stories[[ts]]$allConsistent
    patVar = variants
    patVar$variants = patVar$variants[colnames(stories[[ts]]$consistentClusters$cloneStories$stories)]

    if (side=="L") {
      col = left_color
    } else {
      col = right_color
    }
    xlimnolab = 'default'
    xlimlab = 'default'
  substories = abs(cloneStories$stories)
  rownames(substories) = rownames(cloneStories)
  purity = sapply(1:ncol(cloneStories$stories), function(i) max(abs(substories[names(cloneTree),i])))
  leadingNormal = F
  if ( any(purity==0) ) {
    leadingNormal = T
    purity[purity==0] = 1
  }
  normalise = T

  if ( normalise ) substories = t(t(substories)/purity)
  if ( !leadingNormal ) {
    substories = cbind(rep(0, nrow(substories)), substories)
    colnames(substories)[1] = 'germline'
    if ( any(rownames(substories)=='germline') ) substories['germline','germline'] = 1
  }
  x = 1:ncol(substories)

    if ((patient==6) && (side=='R')) {
      sampleOrder = c("PT6A_R", "PT6B_R", "GL_6")
    } else if ((patient==6) && (side=='L')) {
      sampleOrder = c("PT6A_L", "PT6B_L", "GL_6")
    } else if ((patient==4) && (side=='R')) {
      sampleOrder = c("GL_4", "PT4_R")
      x = 1:ncol(substories)
      xlimnolab = c(2, max(x))
      xlimlab = c(2, max(x)+ceiling(nrow(cloneStories)/2)*1.5)
    } else {
      sampleOrder = 'default'
      x = 1:ncol(substories)
      xlimnolab = c(2, max(x)) 
      xlimlab = c(2, max(x)+ceiling(nrow(cloneStories)/2)*1.5)
    }
    factor = 1
    if ((patient==6) || (patient==3) || ((patient==1) && (side=='L')) || ((patient==5) && (side=='R'))) {
      factor = 2
    }

    pdf(paste('results/superFreq/new_fishplot/', patient, "_", side, '_fish.pdf', sep=''), width=8*factor*1.5, height=8)
    plotRiver(cloneTree=cloneTree, cloneStories=cloneStories,
              storyList=storyList, allStories=allStories,
              variants=patVar, cex=8*factor/11*1.5, genome=genome, annotationMethod=annotationMethod,
              labels=F, colourPool=col, sampleOrder=sampleOrder, xlim=xlimnolab)
    dev.off()

    # pdf(paste('results/superFreq/new_fishplot/', patient, "_", side, '_fish_annot.pdf', sep=''), width=15, height=10)
    # plotRiver(cloneTree=cloneTree, cloneStories=cloneStories,
    #           storyList=storyList, allStories=allStories,
    #           variants=patVar, genome=genome, annotationMethod=annotationMethod,
    #           labels=F, colourPool=col, sampleOrder=sampleOrder, xlim=xlimnolab)
    # dev.off()
}}

# for annotation
#https://stackoverflow.com/questions/31371296/how-to-write-text-inside-a-rectangle-in-r


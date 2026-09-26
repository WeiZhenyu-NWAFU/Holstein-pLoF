# Updated labels are illustrative, not a new ranking or confidence tier.
# No biological function is used to change a test statistic or membership.
# Hand-positioned R/grid figure. Read-only statistical input; no fitting/testing.
# Run from the project root: Rscript analysis/section5_representative_final/plot_figure5.R
library(grid)
Sys.setlocale('LC_CTYPE','English_United States.utf8')
argv <- commandArgs(trailingOnly=FALSE)
script <- normalizePath(sub('^--file=', '', argv[grepl('^--file=',argv)][1]), winslash='/')
root <- dirname(script)
out <- if(dir.exists(file.path(dirname(root),'source_tables'))) dirname(root) else file.path(root,'final')
input <- file.path(out,'source_tables/source_plot_data.tsv')
dir.create(file.path(out,'pdf'), recursive=TRUE, showWarnings=FALSE)
before <- tools::md5sum(input)
d <- read.delim(input, check.names=FALSE, stringsAsFactors=FALSE, na.strings=c('nan','','NE'),fileEncoding='UTF-8')
truth <- function(x) !is.na(x) & tolower(as.character(x)) %in% c('true','1')
getpanel <- function(p, plotted=TRUE) d[d$panel==p & (!plotted | truth(d$plotted)),]
A <- getpanel('A'); B <- getpanel('B'); C <- getpanel('C'); D <- getpanel('D')
stopifnot(all(!truth(A$ring_overlay) | truth(A$primary_support)),
          all(!truth(C$ring_overlay) | truth(C$primary_support)))
# Filled colour = support; thin same-colour outer ring = nested subset.
# Coincident variants remain coincident; no jitter or special candidate sizes.
red <- '#B63247'; blue <- '#197793'; grey <- '#BFC3C7'; ink <- '#303238'
cols <- c('Informative absence'=red,'Homozygotes observed'=blue,'Underpowered absence'=grey)
ff <- 'Arial'
txt <- function(x,y,label,size=8.3,just='left',face='plain',col=ink,rot=0,units='native') {
  grid.text(label,x=unit(x,units),y=unit(y,units),just=just,rot=rot,
            gp=gpar(fontfamily=ff,fontsize=size,fontface=face,col=col))
}
seg <- function(x0,y0,x1,y1,col='#85898D',width=.55,lty=1,units='native') {
  grid.segments(unit(x0,units),unit(y0,units),unit(x1,units),unit(y1,units),
                gp=gpar(col=col,lwd=width*96/72,lty=lty,lineend='butt'))
}
mark <- function(x,y,col,open=FALSE,diam=if(open) 3.9 else 2.35,alpha=1,units='native') {
  if(!length(x)) return(invisible(NULL))
  grid.circle(unit(x,units),unit(y,units),r=unit(diam/2,'points'),
              gp=gpar(fill=if(open) NA else adjustcolor(col,alpha.f=alpha),
                      col=if(open) col else NA,lwd=.45*96/72))
}
cloud <- function(z,col,solid_diam=2.35,ring_diam=3.9,bg_diam=1.35) {
  ring <- truth(z$ring_overlay); solid <- truth(z$primary_support) & !ring
  bg <- !truth(z$primary_support) & !ring
  mark(z$x[bg],z$y[bg],grey,diam=bg_diam,alpha=.26)
  mark(z$x[solid],z$y[solid],col,diam=solid_diam,alpha=.82)
  mark(z$x[ring],z$y[ring],col,diam=solid_diam*.78,alpha=.85)
  mark(z$x[ring],z$y[ring],col,open=TRUE,diam=ring_diam)
}
rects <- list(c(.093,.580,.390,.370),c(.590,.580,.390,.370),
              c(.093,.090,.390,.370),c(.590,.090,.390,.370))
ranges <- list(c(-.015,1.01,-35,110),c(-.08,260,-.15,14),
               c(-.025,1.04,-.1,4.15),c(-.075,.105,-.35,13))
ticksx <- list(seq(0,1,.2),c(0,3,20,100,250),seq(0,1,.2),c(-.05,0,.05,.10))
ticksy <- list(seq(-20,100,20),c(0,1,3,5,10),0:4,seq(0,12,2))
sx <- function(x) sqrt(x+1)-1
axis_panel <- function(k) {
  r <- rects[[k]]; lim <- ranges[[k]]
  if(k==2) lim <- sx(lim)
  pushViewport(viewport(x=r[1],y=r[2],width=r[3],height=r[4],just=c('left','bottom'),
                        xscale=lim[1:2],yscale=lim[3:4],clip='off'))
  seg(0,0,1,0,col=ink,width=.65,units='npc');seg(0,0,0,1,col=ink,width=.65,units='npc')
  xt <- ticksx[[k]]; yt <- ticksy[[k]]
  xp <- if(k==2) sx(xt) else xt; yp <- if(k==2) sx(yt) else yt
  labsx <- if(k %in% c(1,3)) sprintf('%.1f',xt) else if(k==4) c('-0.05','0','0.05','0.10') else as.character(xt)
  for(i in seq_along(xp)) {
    grid.segments(unit(xp[i],'native'),unit(0,'npc'),unit(xp[i],'native'),unit(-2.5,'points'),gp=gpar(col=ink,lwd=.6*96/72))
    grid.text(labsx[i],unit(xp[i],'native'),unit(-5,'points'),just=c('centre','top'),gp=gpar(fontfamily=ff,fontsize=8.4,col=ink))
  }
  for(i in seq_along(yp)) {
    grid.segments(unit(0,'npc'),unit(yp[i],'native'),unit(-2.5,'points'),unit(yp[i],'native'),gp=gpar(col=ink,lwd=.6*96/72))
    grid.text(as.character(yt[i]),unit(-5,'points'),unit(yp[i],'native'),just='right',gp=gpar(fontfamily=ff,fontsize=8.4,col=ink))
  }
  xl <- list('pLoF allele frequency','Expected pLoF homozygotes',
             'Extreme-iHS fraction in 100 kb',expression('Estimated selection coefficient, '*italic(s)))
  yl <- list(expression('Signed '*-log[10](italic(P))),'Observed pLoF homozygotes',
             expression('Empirical '*-log[10](italic(P))),'CLUES logLR')
  grid.text(xl[[k]],unit(.5,'npc'),unit(-22,'points'),gp=gpar(fontfamily=ff,fontsize=9.3,col=ink))
  grid.text(yl[[k]],unit(-32,'points'),unit(.5,'npc'),rot=90,gp=gpar(fontfamily=ff,fontsize=9.3,col=ink))
  pushViewport(viewport(clip='on',xscale=lim[1:2],yscale=lim[3:4]))
}
key <- function(x,y,label,col,open=FALSE) {
  if(open) mark(x,y,col,diam=1.8,units='npc')
  mark(x,y,col,open=open,diam=if(open) 4.0 else 2.8,units='npc')
  txt(x+.035,y,label,size=7.6,units='npc')
}
label <- function(z,id,tx,ty,ex,ey,k,name=NULL) {
  q <- z[z$variant==id,]; stopifnot(nrow(q)==1)
  if(!is.null(name)) q$label <- name
  stopifnot(!is.na(q$label), nzchar(q$label))
  xy <- c(q$x,q$y); if(k==2) {xy <- sx(xy);tx<-sx(tx);ty<-sx(ty);ex<-sx(ex);ey<-sx(ey)}
  seg(xy[1],xy[2],ex,ey,width=.45,col='#7B8084')
  txt(tx,ty,q$label,size=if(startsWith(q$label,'ENSBTAG')) 7.1 else 8.3,
      face=if(startsWith(q$label,'ENSBTAG') | startsWith(q$label,'BTA')) 'plain' else 'italic')
}
inset <- function() {
  e <- getpanel('A_inset');q <- e[truth(e$primary_support),]
  pushViewport(viewport(x=.10,y=.82,width=.32,height=.105,just=c('left','bottom'),
                        xscale=c(.88,1.45),yscale=c(-.65,1.65),clip='off'))
  seg(1,-.65,1,1.5,col='#CBCED1',width=.5,lty=2)
  seg(.88,-.65,1.45,-.65,col=ink,width=.5)
  mark(1,1,'#969CA1',diam=3.1);mark(q$x,0,red,diam=3.1)
  seg(q$CI_low,0,q$CI_high,0,col=red,width=.65)
  for(v in c(q$CI_low,q$CI_high)) seg(v,-.16,v,.16,col=red,width=.65)
  txt(.86,1,'Syn.',size=7.4,just='right');txt(.86,0,'pLoF',size=7.4,just='right')
  for(v in c(1,1.2,1.4)) {
    grid.segments(unit(v,'native'),unit(-.65,'native'),unit(v,'native'),unit(-.65,'native')-unit(1.5,'points'),gp=gpar(lwd=.5,col=ink))
    grid.text(sprintf('%.1f',v),unit(v,'native'),unit(-.65,'native')-unit(3,'points'),just=c('centre','top'),gp=gpar(fontfamily=ff,fontsize=7.2,col=ink))
  }
  txt(.88,2.1,'Odds ratio (95% CI)',size=7.5)
  popViewport()
  txt(.10,.718,sprintf('Empirical P = %.4f',q$empirical_P),size=7.4,units='npc')
}
draw <- function() {
 grid.newpage()
 for(k in 1:4) {
  axis_panel(k)
  if(k==1) {
   seg(-.015,0,1.01,0,col='#898E92',width=.5);cloud(A,red)
   label(A,'10:23541774:G:T',.56,73,.54,68,1)
   label(A,'5:102790526:T:C',.51,46,.49,46,1,'WC1')

  }
  if(k==2) {
   seg(sx(0),sx(0),sx(14),sx(14),col='#A7ACB0',width=.6,lty=2)
   seg(sx(3),sx(-.15),sx(3),sx(14),col='#DEE0E2',width=.4,lty=3)
   for(st in names(cols)) {q<-B[B$status==st,];mark(sx(q$x),sx(q$y),cols[st],diam=3.7,alpha=.85)}
   label(B,'5:102790526:T:C',169,4.1,211,3.9,2,'WC1')
   label(B,'4:105599505:A:T',97,.50,144,.24,2)

  }
  if(k==3) {
   seg(-.025,2,1.04,2,col=adjustcolor(blue,.6),width=.6,lty=2);cloud(C,blue,solid_diam=3.4,ring_diam=5.2,bg_diam=1.8)
   label(C,'12_21367333',.71,3.82,.86,3.66,3,'ATP7B')
   label(C,'7_52490114',.69,2.60,.86,2.65,3,'PCDHGC3')
   label(C,'4_113113368',.31,3.12,.81,3.22,3,'BTA4:113.1\u2013113.2 Mb\n(2 pLoFs)')
  }
  if(k==4) {
   seg(-.075,1.2661,.105,1.2661,col=adjustcolor(red,.6),width=.6,lty=2)
   seg(0,-.35,0,13,col='#DEE0E2',width=.45);cloud(D,red,solid_diam=3.4,bg_diam=1.8)
   label(D,'29_41147975',.024,11.8,.021,11.75,4)
   label(D,'7_11145809',.014,9.6,.015,9.23,4)
   label(D,'5_31184185',.021,4.05,.024,3.7,4)
  }
  popViewport()
  # Legend symbols are legible at print scale and reflect mutually exclusive
  # drawing categories, while the statistical sets remain nested.
  if(k==1) {
   inset()
   key(.51,.965,'Other pLoFs',grey)
   key(.51,.905,'Depleted, FDR < 0.05',red)
   key(.51,.845,'No homozygotes, E \u2265 3',red,TRUE)
  }
  if(k==2) {
   nb<-table(getpanel('B',FALSE)$status)
   key(.40,.965,sprintf('Informative absence (%s)',nb['Informative absence']),red)
   key(.40,.905,sprintf('Homozygotes observed (%s)',nb['Homozygotes observed']),blue)
   key(.40,.845,sprintf('Underpowered absence (%s)',nb['Underpowered absence']),grey)
   txt(.435,.765,sprintf('Not evaluable: %s',nb['Not evaluable']),size=7.6,units='npc')
  }
  if(k==3) {
   key(.045,.965,'Other evaluated pLoFs',grey)
   key(.045,.905,sprintf('Regional support (%s/%s)',sum(truth(C$primary_support)),nrow(C)),blue)
   key(.045,.845,sprintf('All three scales (%s)',sum(truth(C$ring_overlay))),blue,TRUE)
  }
  if(k==4) {
   key(.045,.965,'Other pLoFs',grey)
   key(.045,.905,sprintf('CLUES: %s/%s',sum(truth(D$primary_support)),nrow(D)),red)
  }
  popViewport()
  r<-rects[[k]]
  txt(r[1]-.055,r[2]+r[4]+.020,LETTERS[k],size=12,face='bold',units='npc')
 }
}
name<-'Fig5_final'
cairo_pdf(file.path(out,paste0(name,'.pdf')),width=175/25.4,height=170/25.4,family=ff,onefile=TRUE)
draw();dev.off()
svg(file.path(out,paste0(name,'.svg')),width=175/25.4,height=170/25.4,family=ff)
draw();dev.off()
png(file.path(out,paste0(name,'.png')),width=175,height=170,units='mm',res=240,type='cairo',bg='white')
draw();dev.off()
stopifnot(identical(before,tools::md5sum(input)))

# Statistics and coordinates remain unchanged; labels are recorded in the final source.
cat('Frozen source unchanged. All supported/subset memberships preserved.\n')

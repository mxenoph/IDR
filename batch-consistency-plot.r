# 1-20-10 Qunhua Li
#
# This program first plots correspondence curve and IDR threshold plot
# (i.e. number of selected peaks vs IDR) for each pair of sample
#
# usage: 
# Rscript batch-consistency-plot-merged.r [npairs] [output.dir] [input.file.prefix 1, 2, 3 ...]
# [npairs]: integer, number of consistency analyses
#          (e.g. if 2 replicates, npairs=1, if 3 replicates, npairs=3
# [output.dir]: output directory for plot
# [input.file.prefix 1, 2, 3]: prefix for the output from batch-consistency-analysis2. They are the input files for merged analysis see below for examples (i.e. saved.file.prefix). It can be multiple files
#

args <- commandArgs(trailingOnly=T)
npair <- args[1] # number of curves to plot on the same figure
output.file.prefix <- args[2] # file name for plot, generated from script at the outer level
df.txt <- 10
ntemp <- as.numeric(npair)
saved.file.prefix <- list() # identifier of filenames that contain the em and URI results

path_to_functions = system('which functions-all-clayton-12-13.r', intern=T)
source(path_to_functions)

uri.list <- list()
uri.list.match <- list()
ez.list <- list()
legend.txt <- c()
em.output.list <- list()
uri.output.list <- list()

for(i in 1:npair){
  saved.file.prefix[i] <- args[2+i]
 
  load(paste(saved.file.prefix[i], "-uri.sav", sep=""))
  load(paste(saved.file.prefix[i], "-em.sav", sep=""))

  uri.output.list[[i]] <- uri.output
  em.output.list[[i]] <- em.output

  ez.list[[i]] <- get.ez.tt.all(em.output, uri.output.list[[i]]$data12.enrich$merge1,
                                uri.output.list[[i]]$data12.enrich$merge2) # reverse =T for error rate

  # URI for all peaks
  uri.list[[i]] <- uri.output$uri.n
  # URI for matched peaks
  uri.match <- get.uri.matched(em.output$data.pruned, df=df.txt)
  uri.list.match[[i]] <- uri.match$uri.n

  file.name <- unlist(strsplit(as.character(saved.file.prefix[i]), "/"))
  
  legend.txt[i] <- paste(i, "=", file.name[length(file.name)])

}

plot.uri.file <- paste(output.file.prefix, "-plot.ps", sep="")

############# plot and report output
# plot correspondence curve for each pair,
# plot number of selected peaks vs IDR 
# plot all into 1 file
postscript(paste(output.file.prefix, "-plot.ps", sep=""))
par(mfcol=c(2,3), mar=c(5,6,4,2)+0.1)
plot.uri.group(uri.list, NULL, file.name=NULL, c(1:npair), title.txt="all peaks")
plot.uri.group(uri.list.match, NULL, file.name=NULL, c(1:npair), title.txt="matched peaks")
plot.ez.group(ez.list, plot.dir=NULL, file.name=NULL, legend.txt=c(1:npair), y.lim=c(0, 0.6))
plot(0, 1, type="n", xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", xaxt="n", yaxt="n") # legends
legend(0, 1, legend.txt, cex=0.6)

dev.off()

# plots also the r standard plots, should really write a function to gather the data
# rather than hacking the idr functions
gg_param = modules::import('ggplots')
pdf(paste(output.file.prefix, "-ggplot.pdf", sep=""), paper='a4')

plot_diagnostics = function(points, lines, slopes, title=NULL){# {{{
    library(ggplot2)
    source("~/source/Rscripts/ggplot-functions.R")
    # Setting colour=NA outside the aes removes the black border from points
    p = ggplot(data=points, aes(x=tv, y=uri, fill=comparison))
    # x axis is all significant peaks and y is the common in each comparison so max(x) is the limit for both axes
    p = p + geom_line(data=lines, aes(x=x, y=y, colour=comparison), size=2.5) + ylim(0, max(points[,1:2]))
    p = p + ggtitle(title) + gg_param$theme_publication
    legend = g_legend(p)
    # adding points after the legend is saved so that there is no point printed in
    # legend
    p = p + geom_point(pch=21, size=2.5, colour=NA)
    p = p + geom_abline(intercept=0, slope=1, linetype="dotted")
    p = p + xlab('# significant peaks') + ylab('# peaks in common') + theme(legend.key=element_blank())
    p = p + theme(legend.position="none")

    y = ggplot(data=slopes, aes(x=x, y=y, colour=comparison)) + geom_line(size=2.5)
    y = y + geom_abline(intercept=1, slope=0, linetype="dotted")
    y = y + xlab('# significant peaks') + ylab('slope')
    y = y + gg_param$theme_publication + theme(legend.position="none") + ylim(0, 1.5) + ggtitle(title)

    return(list('p' = p, 'y' = y, 'leg'=legend))
}
# }}}

grobs= list()
summarised_data = plot.uri.group(uri.list, NULL, file.name=NULL, legend.txt, title.txt="all peaks")
x = plot_diagnostics(summarised_data$data_points, summarised_data$data_lines, summarised_data$data_slopes, 'all peaks')
grobs[[1]] = x[['p']]
grobs[[2]] = x[['y']]

summarised_data = plot.uri.group(uri.list.match, NULL, file.name=NULL, legend.txt, title.txt="matched peaks")
x = plot_diagnostics(summarised_data$data_points, summarised_data$data_lines, summarised_data$data_slopes, 'matched peaks')
grobs[[3]] = x[['p']]
grobs[[4]] = x[['y']]

summarised_data = plot.ez.group(ez.list, plot.dir=NULL, file.name=NULL, legend.txt=legend.txt, y.lim=c(0, 0.6))
grobs[[5]] = ggplot(summarised_data, aes(x=n.plot, y=IDR.plot, fill=comparison)) + geom_point(pch=21, colour=NA) + geom_line() + xlab('# significant peaks') + ylab('IDR') + gg_param$theme_publication + theme(legend.position="none")
grobs[[6]] = x[['leg']]

library(cowplot)
plot_grid(grobs[[1]], grobs[[2]],  grobs[[3]],
          grobs[[4]], grobs[[5]],  grobs[[6]],
          labels = toupper(letters[1:6]), ncol = 2, align = 'v')
dev.off()

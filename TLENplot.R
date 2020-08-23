
### Tian Huang, 08112020

library("argparse")
options(error=traceback)

parser <- ArgumentParser(description='Plot template length (TLEN) distribution for .bam files.')
parser$add_argument('input',  help="Input one .bam file or one .csv file (1st col: TLEN (sorted); 2rd col: frequency, the header should not be present)")
parser$add_argument('out_dir',  help='Output directory plus prefix')

parser$add_argument('-p','--plot', nargs="+",
                    help = 'Specify the plot type; you can choose either "ecdf" or "hist" (default: "hist")',
                    required=FALSE, default="hist")
parser$add_argument('-f','--format', nargs="+",
                    help = 'Spesify the output figure format (default: "png")',
                    required=FALSE, default="png")
parser$add_argument('-s','--size', nargs="+",
                    help = "Specify figure size (width and height); default: 100 60",
                    required=FALSE, default=c(100,60))
parser$add_argument('-u','--unit', nargs="+",
                    help = 'Specify figure size unit (default: "mm")',
                    required=FALSE, default="mm")
parser$add_argument('-x','--xlim', nargs="+",
                    help = 'Specify x axis range in the figure (default: 0 500)',
                    required=FALSE, default=c(0,500))
parser$add_argument('-t','--table', nargs="+",
                    help = 'When the input is a bam file, this parameter determines whether the script should output an extra table of the length distribution (default: TRUE)',
                    required=FALSE, default="TRUE")


args <- parser$parse_args()


library("Rsamtools")
library("ggplot2")

input<-args$input
xlim=as.numeric(args$xlim)

if(endsWith(input,".bam")){

  ### process .bam input
  bamFile<-BamFile(input)
  aln<-scanBam(bamFile)[[1]]
  isize<-abs(aln$isize)
  df<-data.frame(TLEN=isize)
  
  if(args$table=="TRUE") {
    tb=as.matrix(table(isize))
    tb=cbind(as.numeric(row.names(tb)),tb[,1])
    row.names(tb)=NULL
    
    write.table(tb, paste(args$out_dir,".csv",sep=""), row.names = F, col.names = F, sep = ",")
  }
  
  if(args$plot!="hist"){
    p <- ggplot(df,aes(x=TLEN)) +
      stat_ecdf(alpha=0.7) +
      coord_trans(xlim=c(xlim[1],xlim[2])) +
      labs(x="Template length (bp)", y="Cumulative frequency")+
      theme(axis.title.y=element_text(size=9),
            axis.title.x=element_text(size=9),
            legend.text = element_text(size=8),
            legend.title = element_text(size=9))
  } else{
    p <- ggplot(df,aes(x=TLEN)) +
      geom_histogram( #aes(y=..density..),
        color="black",fill="white", breaks=seq(xlim[1],xlim[2],round((xlim[2]-xlim[1])/50))) +
      #geom_density(size=1,color="#2D608A",bw=round((xlim[2]-xlim[1])/50)) +
      coord_trans(xlim=c(xlim[1],xlim[2])) +
      labs(x="Template length (bp)", y="Counts")+
      theme(axis.title.y=element_text(size=9),
            axis.title.x=element_text(size=9),
            legend.text = element_text(size=8),
            legend.title = element_text(size=9))
  }

} else{
  
  ### process .csv data
  library("tidyr")
  
  tb=read.csv(input,header = F, colClasses = c("integer","integer"))

  if(args$plot!="hist"){
    ### turn the counts into freq
    tb[,2]<-cumsum(tb[,2])
    tb[,2]<-tb[,2]/max(tb[,2])
    
    df=data.frame(TLEN=tb[,1],freq=tb[,2])
    
    p <- ggplot(df,aes(x=TLEN,y=freq)) +
      geom_step(alpha=0.7) +
      coord_trans(xlim=c(xlim[1],xlim[2])) +
      labs(x="Template length (bp)", y="Cumulative frequency")+
      theme(axis.title.y=element_text(size=9),
            axis.title.x=element_text(size=9),
            legend.text = element_text(size=8),
            legend.title = element_text(size=9))
  } else{
    ### uncount the table
    df=data.frame(TLEN=tb[,1],freq=tb[,2])
    df=uncount(df,freq)
    
    p <- ggplot(df,aes(x=TLEN)) +
      geom_histogram( #aes(y=..density..),
        color="black",fill="white", breaks=seq(xlim[1],xlim[2],round((xlim[2]-xlim[1])/50))) +
      #geom_density(size=1,color="#2D608A",bw=round((xlim[2]-xlim[1])/50)) +
      coord_trans(xlim=c(xlim[1],xlim[2])) +
      labs(x="Template length (bp)", y="Counts")+
      theme(axis.title.y=element_text(size=9),
            axis.title.x=element_text(size=9),
            legend.text = element_text(size=8),
            legend.title = element_text(size=9))
  }
}


ggsave(filename = paste(args$out_dir,".",args$format,sep=""),plot = p,width = as.numeric(args$size[1]),height = as.numeric(args$size[2])
       ,units = args$unit)
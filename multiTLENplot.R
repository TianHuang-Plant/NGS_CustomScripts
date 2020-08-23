
### Tian Huang, 0812020

library("argparse")
options(error=traceback)

parser <- ArgumentParser(description='Plot template length (TLEN) cumulative curve for .bam files.')
parser$add_argument('-i','--input', nargs="+", help=".bam or .csv files",required=TRUE)
parser$add_argument('-o','--out_dir', help='Output directory',required=TRUE)
parser$add_argument('-N','--figName', help='Name the ouput figure',required=TRUE)


parser$add_argument('-S','--suffix', help='Specify the number of characters in the input suffix and name the samples using trimmed file names (ex. 4 for .bam files)',required=FALSE,default=0)
parser$add_argument('-n','--name', help='Name the input samples sequentially ("-s" will override this arguments)',required=FALSE,default=c())
parser$add_argument('-p','--plot', nargs="+",
                    help = 'Specify the plot type; you can choose either "ecdf" or "hist" (default: "ecdf")',
                    required=FALSE, default="ecdf")
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
parser$add_argument('--facet', nargs="+",
                    help = 'Change the facet col number in hist output (default: 4)',
                    required=FALSE, default=4)


tlen_hist <- function(df,xlim){
  ggplot(df,aes(x=TLEN)) +
    geom_histogram(breaks=seq(xlim[1],xlim[2],round((xlim[2]-xlim[1])/30)),color="black",fill="white") +
    facet_wrap(vars(sample), ncol = facet_ncol, scales = "free") +
    coord_trans(xlim=c(xlim[1],xlim[2])) +
    scale_y_continuous(labels = function(x){100*x}) +
    labs(x="Template length (bp)", y="Counts")+
    theme(axis.title.y=element_text(size=9),
          axis.title.x=element_text(size=9),
          legend.text = element_text(size=8),
          legend.title = element_text(size=9))
}

tlen_ecdf <- function(df,xlim){
  ggplot(df,aes(x=TLEN,color=sample)) +
    stat_ecdf(size=1) +
    coord_trans(xlim=c(xlim[1],xlim[2])) +
    scale_y_continuous(labels = scales::percent) +
    labs(x="Template length (bp)", y="Cumulative frequency")+
    theme(axis.title.y=element_text(size=9),
          axis.title.x=element_text(size=9),
          legend.text = element_text(size=8),
          legend.title = element_text(size=9))
}

compress1k <- function(isize){
  tb<-as.matrix(table(isize))
  temp=data.frame(TLEN=tb[,1],freq=round(tb[,2]/1000))
  uncount(temp,freq)
}

bamToTLEN <- function(input){
  bamFile<-BamFile(input)
  aln<-scanBam(bamFile)[[1]]
  abs(aln$isize)
}

args <- parser$parse_args()


library("Rsamtools")
library("ggplot2")

print("Loading files...")

input<-args$input
xlim=as.numeric(args$xlim)
inputLength<-length(input)
facet_ncol<-args$facet

if(!endsWith(args$out_dir,"/")){
  out_dir<-paste(out_dir,"/",sep="")
} else {out_dir<-args$out_dir}

if(as.numeric(args$suffix)==0 && length(args$names)==inputLength){
  sampleName<-args$names
}else{
  sampleName<-substr(input,max(gregexpr("/",input)[[1]])+1,nchar(input)-as.numeric(args$suffix))
}

isize<-matrix(nrow=0,ncol=2)

if(sum(endsWith(input,".bam"))==inputLength){
  
  ############################## process .bam input
  for(i in 1:inputLength){
    temp<-bamToTLEN(input[i])
    print(paste(i,"file(s) loaded"))
    
    if(args$table=="TRUE") {
      tb=as.matrix(table(temp))
      tb=cbind(as.numeric(row.names(tb)),tb[,1])
      row.names(tb)=NULL
      
      write.table(tb, paste(out_dir,sampleName[i],".csv",sep=""), row.names = F, col.names = F, sep = ",")
      print(paste(i,"csv file(s) written"))
      
    }
    
    temp<-compress1k(temp)
    isize<-rbind(isize,cbind(temp,sampleName[i]))
  }
  df<-data.frame(TLEN=as.numeric(isize[,1]),
                 sample=isize[,2])
  
  if(args$plot=="ecdf"){
    p <- tlen_ecdf(df,xlim)
  } else{
    p <- tlen_hist(df,xlim)
  }
  
} else if(sum(endsWith(input,".csv"))==inputLength){
  
  ############################### process .csv data
  library("tidyr")
  
  for(i in 1:inputLength){
    ### uncount the table
    tb=read.csv(input[i],header = F, colClasses = c("integer","integer"))
    temp=data.frame(TLEN=tb[,1],freq=round(tb[,2]/1000))
    temp=uncount(temp,freq)
    isize<-rbind(isize,cbind(temp,sampleName[i]))
    
    print(paste(i,"file(s) loaded"))
  }
  df<-data.frame(TLEN=as.numeric(isize[,1]),
                 sample=isize[,2])
  if(args$plot=="ecdf"){
    p <- tlen_ecdf(df,xlim)
  } else{
    p <- tlen_hist(df,xlim)
  }
  

} else{print('Wrong files!')}


ggsave(filename = paste(out_dir,args$figName,".",args$format,sep=""),plot = p,width = as.numeric(args$size[1]),height = as.numeric(args$size[2])
       ,units = args$unit)

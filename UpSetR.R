### Tian Huang, 0823020

library("argparse")
options(error=traceback)

parser <- ArgumentParser(description='Visualization of set intersections using UpSet.')
parser$add_argument('-i','--input', nargs="+", help="Tab-deliminated tables",required=TRUE)
parser$add_argument('-c','--col', nargs="+", help="The column index",required=TRUE)
parser$add_argument('-o','--output', help='Output directory plus prefix',required=TRUE)

parser$add_argument('-n','--name', help='Name of each set',required=F)
parser$add_argument('-s','--suffix', help='Specify the number of characters in the input suffix and name the samples using trimmed file names (ex. 4 for .bam files)',required=FALSE,default=0)
parser$add_argument('--nintersects', help='Number of intersections to plot',required=F,default=NA)

args <- parser$parse_args()

library("UpSetR")

input<-args$input
col<-args$col
lt<-list()

if(as.numeric(args$suffix)==0 && length(args$names)==length(input)){
  names<-args$names
}else{
  names<-substr(input,max(gregexpr("/",input)[[1]])+1,nchar(input)-as.numeric(args$suffix))
}

for(i in length(names)){
  tb<-read.table(input[i],sep = "\t",header = F)
  lt[names[i]]<-tb[,col]
}

png(paste(args$output,".png",sep=""))
upset(fromList(lt), nintersects = args$nintersects, mainbar.y.label = "log2(Counts)",
         sets.x.label = "Counts", scale.intersections = "log2",mb.ratio=c(0.55,0.45), order.by = "freq")
dev.off()

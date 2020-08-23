### Tian Huang, 0823020

library("argparse")
options(error=traceback)

parser <- ArgumentParser(description='Visualization of set intersections using UpSet.')
parser$add_argument('-i','--input', nargs="+", help="Tab-deliminated tables",required=TRUE)
parser$add_argument('-c','--col', help="The column index",required=TRUE)
parser$add_argument('-o','--output', help='Output directory plus prefix',required=TRUE)

parser$add_argument('-n','--name',nargs="+", help='Name of each set',required=FALSE)
parser$add_argument('-S','--suffix', help='Specify the number of characters in the input suffix and name the samples using trimmed file names (ex. 4 for .bam files)',required=FALSE,default=0)
parser$add_argument('--nintersects', help='Number of intersections to plot',required=FALSE,default=10)
parser$add_argument('-s','--size', nargs="+",
                    help = "Specify figure size (width and height); default: 7 7",
                    required=FALSE, default=c(7,7))

args <- parser$parse_args()

library("UpSetR")

input<-args$input
col<-as.integer(args$col)
ni<-args$nintersects
size<-as.numeric(args$size)
lt<-list()

if(as.numeric(args$suffix)==0 && length(args$name)==length(input)){
  names<-args$name
}else{
  names<-substr(input,max(gregexpr("/",input)[[1]])+1,nchar(input)-as.numeric(args$suffix))
}

for(i in 1:length(names)){
  tb<-read.table(input[i],sep = "\t",header = T)
  lt[[names[i]]]<-tb[,col]
}

svg(file=paste(args$output,".svg",sep=""), width=size[1],height=size[2])
upset(fromList(lt), sets=names, keep.order=T, nintersects = ni, mainbar.y.label = "Counts",
         sets.x.label = "Counts", mb.ratio=c(0.55,0.45), order.by = "freq")
dev.off()

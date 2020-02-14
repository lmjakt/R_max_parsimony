## test functions in somewhere..

## use a matrix of distances to get a tree from ape
## make these from a set of intron sizes using whatever
## distance we like.

require(ape)

## we have some distance data as well as the intron sizes that were used to
## infer those distances. This if for a limited set of 20 species, but rather
## a large data number of intron sizes.

dists <- as.matrix( read.table( "mi_based_distances_cr.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE ) )


dists.nj <- nj( dists )
plot(dists.nj)

names(dists.nj)
## [1] "edge"        "edge.length" "tip.label"   "Nnode"

## edge contains the edge information..
## Nnode the number of (internal?) nodes (18)
length(unique(as.numeric(dists.nj$edge)))
## [1] 38
## that is 20 + 18, so Nnode does not include the
## leaf nodes.

## note that
length(dists.nj$edge.length)
## 37, which means we have one edge less than the number of nodes
## which indicates an unrooted tree.

## let us read in the intron size data used to make the tree
int.s <- read.table("intron_orthology_length_dr.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
dim(int.s)
## [1] 65717    20

## no we are more interested in the log transformed data than the reality, so we can do
int.sl <- log2( int.s )

range(int.sl, na.rm=TRUE)
## [1]  0.00000 19.43276

## now we want to discretise this data, into a set of values that can be represented as
## character vectors. How do we get the character vectors in R.
## use
## intToUtf8()
## utf8toInt()
## for conversions..

## Note that we can use a simple encoding starting at A
## A is 65
## @ is 64
## Then NA / inf -> 64, which can both be considered to forms of missing data
## Note that 0 values should probably be considered as missing values as the
## most likely issue may be an alternative splicing or incorrect gene annotation

range(as.integer(as.matrix(int.sl)), na.rm=TRUE )
## [1]  0 19

## we could then simply use 64 + as.integer after converting non-numeric values to
## 0

encode.dist <- function(data, base=64, conversion=as.integer){
    data <- apply(data, 2, conversion)
    data[ is.na(data) ] <- 0
    data <- data + base
    apply(data, 2, intToUtf8)
}

int.sa <- encode.dist(int.sl)
range(sapply(int.sa, function(x){ range( utf8ToInt(x) )} ))
## 64 83

## to look at them..
sapply( int.sa, substr, 1, 10 ) ## looks good.. 

## this will only return a dummy variable, but for testing purposes we can see how well
## it works.
make.sub.matrix <- function(size){
    m <- matrix(nrow=size, ncol=size)
    for(i in 1:size){
        for(j in 1:size){
            m[i,j] <- as.integer(abs(i-j))
        }
    }
    m
}

## the range of values is between 64 -> 83

sub.matrix <- make.sub.matrix( 1 + 83 - 64 )
sub.matrix[1,] <- 20L
sub.matrix[,1] <- 20L
sub.matrix[1,1] <- 0L

## missing intron lengths, or intron lengths that are less than 2 will
## here be considered as missing. That is, that there is no penalty to merge a missing
## intron with a not missing one. Note that it might be an idea for anything shorter than
## a given length to be used for this.. 

dyn.load( "src/max_parsimony.so" )

tmp <- .Call("sankoff", dists.nj$edge, c(38L, 20L), sub.matrix, c(64L, 20L), int.sa )
             

## to confirm the behaviour of this I want to draw with more labels..

nj.lines <- function(tree){
    y <- vector(mode='numeric', length=nrow(tree$edge))
    x <- matrix(0, nrow=nrow(tree$edge), ncol=2)
    v.lines <- matrix(nrow=0, ncol=3)
    leaf.b <- tree$edge[,2] < min(tree$edge[,1])
    y[leaf.b] <- 1:sum(leaf.b)
    ## start from the root...
    ## the root has no parent..
    root <- setdiff( tree$edge[,1], tree$edge[,2] )
    visit.tree <- function(root, r.x ){
        ## child.i is wrong. lets call it this.i
        this.i <- which( tree$edge[,2] == root )
        root.i <- which( tree$edge[,1] == root )
        if(!length(root.i)){
            return(y[ this.i ])
        }
        x[root.i,1] <<- r.x
        x[root.i,2] <<- x[root.i,1] + tree$edge.length[root.i]
        children <- tree$edge[root.i, 2]
        child.y <- vector(length=length(root.i))
        for(i in 1:length(root.i))
            child.y[i] <- visit.tree( children[i], x[root.i[i], 2] )
        y[ this.i ] <<- mean(child.y)
        v.lines <<- rbind(v.lines, c(x[ root.i[1], 1], min(child.y), max(child.y)))
        return( y[this.i] )
    }
    visit.tree(root, 0)
    list('x'=x, 'y'=y, 'v'=v.lines)
}

tmp.l <- nj.lines(dists.nj)
usr <- par("usr")

dev.set(3)
plot.new()
plot.window( xlim=usr[1:2], ylim=usr[3:4], xaxs='i', yaxs='i')
segments( tmp.l$x[,1], tmp.l$y, tmp.l$x[,2], tmp.l$y )

segments( tmp.l$v[,1], tmp.l$v[,2], tmp.l$v[,1], tmp.l$v[,3] )

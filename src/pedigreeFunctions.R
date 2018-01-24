############################################################################################
# Functions to analyze and visualize pedigree data, building on 'synbreed' package
#     Tom Poorten
#     Feb. 2016
#########
# Change plot.pedigree function to allow for control of colors
plot.pedigree2 = 
function (x, effect = NULL, vLabel = pedigree$ID,...) 
{
  if (any(class(x) == "gpData")) 
    pedigree <- x$pedigree
  else pedigree <- x
  if (min(pedigree$gener) != 0) 
    pedigree$gener <- pedigree$gener - min(pedigree$gener)
  if (!is.null(effect) & length(effect) != nrow(pedigree)) 
    stop("length of effect does not equal nrow(pedigree)")
  pedigree[!(pedigree$Par1 %in% pedigree$ID), ]$Par1 <- 0
  pedigree[!(pedigree$Par2 %in% pedigree$ID), ]$Par2 <- 0
  relations <- rbind(as.matrix(pedigree[pedigree$Par1 != 0, 
                                        c("Par1", "ID")], ncol = 2), as.matrix(pedigree[pedigree$Par2 != 
                                                                                          0, c("Par2", "ID")], ncol = 2))
  ped.graph <- graph.data.frame(relations, directed = TRUE, 
                                pedigree)
  gener <- pedigree$gener
  n <- nrow(pedigree)
  pos <- matrix(data = NA, nrow = n, ncol = 2)
  pos[, 2] <- max(gener) - gener
  if (is.null(effect)) 
    pos[, 1] <- order(gener, partial = order(pedigree$ID, 
                                             decreasing = TRUE)) - cumsum(c(0, table(gener)))[gener + 
                                                                                                1]
  else pos[, 1] <- effect
  myscale <- function(x) {
    if (length(x) == 1) 
      x <- 0
    else {
      x <- unlist(scale(x))
    }
    return(x)
  }
  if (is.null(effect)) 
    pos[n:1, 1] <- unlist(tapply(pos[, 1], pos[, 2], myscale))
  # cols <- rep("lightblue", n)
  if (!is.null(pedigree$sex)) 
    cols[pedigree$sex == 0] <- "palevioletred1"
  plot(ped.graph, rescale = TRUE, vertex.label = vLabel, 
       layout = pos, edge.color = 1, edge.width = 0.5, edge.arrow.size = 0.5, 
       vertex.label.family = "sans", 
       # vertex.color = cols, 
       ...)
  if (!is.null(effect)) 
    axis(side = 1, at = seq(-1, 1, length = 10), labels = round(seq(min(pos[, 
                                                                            1]), max(pos[, 1]), length = 10), 0))
  return(ped.graph)
}

############################################################################
############################################################################
## Functions to access ascendants and descendants in a pedigree object
#     getAscendants
#     getDescendants
#############
# function - get ascendants
getAscendants = function(focal=NULL,pedigree=NULL){
  parents = c(pedigree$Par1[which(pedigree$ID==focal)],pedigree$Par2[which(pedigree$ID==focal)])
  check = TRUE
  while(check){
    checkLength = length(parents)
    parents = na.omit(unique(c(parents,pedigree$Par1[match(parents, pedigree$ID)],pedigree$Par2[match(parents, pedigree$ID)])))
    # set stop criteria for recursion - stop when 'parents' is no longer growing bc founders are reached
    if(length(parents) == checkLength){
      check = FALSE
    }
  }
  parents = parents[-which(parents == "0")]
  return(parents)
}

#
#############
# function - get ascendants, version 2 - also return # of generations back of ascendants
getAscendants2 = function(focal=NULL,pedigree=NULL){
  colnames(pedigree) = c("ID","Par1","Par2")	
  parents = c(pedigree$Par1[which(pedigree$ID==focal)],pedigree$Par2[which(pedigree$ID==focal)])
  gensBack = c(1,1)
  check = TRUE
  j = 2
  while(check){
    checkLength = length(parents)
    parents = na.omit(unique(c(parents,pedigree$Par1[match(parents, pedigree$ID)],pedigree$Par2[match(parents, pedigree$ID)])))
    # print(parents)
    if(paste(parents, collapse="") == c("0")){
      parents = c("0","0")
    }
    if(length(parents) - checkLength > 0){
      gensBack = c(gensBack, rep(j, length(parents) - checkLength))  
    }
    # set stop criteria for recursion - stop when 'parents' is no longer growing bc founders are reached
    if(length(parents) == checkLength){
      check = FALSE
    }
    j = j+1
  }
  gensBack = gensBack[which(parents != "0")]
  parents = parents[which(parents != "0")]
  return(list(ascendants = parents, gensBack = gensBack))
}

#
#############

#############
# function - get descendants
getDescendants = function(focal=NULL,pedigree=NULL){
  # keep = focal
  descendants = pedigree$ID[unique(c(which(focal == pedigree$Par1), which(focal == pedigree$Par2)))]
  check = TRUE
  while(check){
    checkLength = length(descendants)
    newDescendants = NULL
    for(i in 1:length(descendants)){
      newDescendants = c(newDescendants, pedigree$ID[unique(c(which(descendants[i] == pedigree$Par1), which(descendants[i] == pedigree$Par2)))])
    }
    descendants = unique(c(descendants, newDescendants))
    (descendants)
    length(descendants)
    # set stop criteria for recursion - stop when 'descendants' is no longer growing
    if(length(descendants) == checkLength){
      check = FALSE
    }
  }
  return(descendants)
}

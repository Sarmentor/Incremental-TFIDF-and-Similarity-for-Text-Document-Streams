library("Rcpp")
library("fastmatch")
library("microbenchmark")
library("tm")
library("igraph")


'%!in%' <<- Negate('%in%')



##################################################################################
################## Incremental SIM. MATRIX Function  #############################
##################################################################################


SearchNeighbors <- function(changed.node, order=1){
  
  #search for neighbors of new word
  #or changed word's tfidf weight
  
  if(order==1){
    #1st order neighbors
    neigh <- ego(graph, 1, changed.node, mindist=1)
  }else{
    #2nd order neighbors
    neigh <- ego(graph, 2, changed.node,mindist = 2)
  }
  
  return(neigh)
}


#Function to add nodes (docs and words) to graph
my.addNode <- function(node, word){
  
  
  #Adds nodes if they do not exist in the graph
  is.node.in.graph <-fmatch(V(graph)$name, node)
  
  if(1 %!in% is.node.in.graph){
    #graph <<- graph + vertices(node1.origin)
    if(word==FALSE){
      #new node is a doc
      graph <<- graph + vertices(node, color=1, shape="rectangle", type=1)  
    }else{
      #new node is a word
      graph <<- graph + vertices(node, color=2, shape="circle", type=0)
    }
    
  }
  
  #adds connection between those nodes
  #graph <<- graph + edge(node1.origin, node2.destination)
  
  
}


ConstructTextGraph <- function(doc, node, word) { 
  
  
  
  if(word==0){
    #Add doc
    my.addNode(doc, word)
  }else{
    my.addNode(node, word)
    #Add edge
    graph <<- graph + edge(node, doc)
  }
      
  number.of.nodes <<- number.of.nodes + 1
  graph <<- graph
}

PreUpdateSim <- function(word, doc){
  
  neigh <- SearchNeighbors(word, order = 1)[[1]]$name
  
  
  if(length(neigh[-which(neigh==paste(doc))])==0){
    ##the only 2nd order neighbor is himself
    ##thus, no similarity calculation is needed
    return()
    
  }
  
  #combination of pairs of neighbors
  m.neigh <- t(combn(unique(neigh), 2))
  
  #out with loops between docs
  indexes <- which(m.neigh[,1]==m.neigh[,2])
  if(length(indexes)!=0){
    m.neigh.doc <<- unique(matrix(m.neigh[-indexes,],ncol=2,byrow=TRUE))
  }else{
    m.neigh.doc <<- m.neigh
  }
  
  #out with duplicated undirected connections between docs
  indexes.dup <- which(duplicated(t(apply(m.neigh.doc, 1, sort)))==TRUE)
  if(length(indexes.dup)!=0){
    m.neigh.doc <<- m.neigh.doc[-indexes.dup,]
  }
  
  
  
}


#Incremental Similarity Update
UpdateSim <- function(m.neigh.doc){

    if(nrow(m.neigh.doc)==0){
      return()
    }
    
    for(n in 1:nrow(m.neigh.doc)){
      #Calculate incremental similarity
      #change between pairs
      doca <- m.neigh.doc[n,1]
      docb <- m.neigh.doc[n,2]
      
      #find same words in doca and docb
      #this info is in the TF-IDF sparse list
      a <- SearchNeighbors(doca, order=1)[[1]]$name
      b <- SearchNeighbors(docb, order=1)[[1]]$name   
      same.words <- a[a %in% b]   
      
      same.words.vector.doca <- c()
      same.words.vector.docb <- c()
      words.vector.doca <- c()
      words.vector.docb <- c()
      
      sapply(same.words, FUN=function(x){
        same.words.vector.doca <<- c(same.words.vector.doca, dtm[[doca]][[x]][paste("TFIDF")]) 
        same.words.vector.docb <<- c(same.words.vector.docb, dtm[[docb]][[x]][paste("TFIDF")])
      })
      
      lapply(dtm[[doca]], FUN=function(x){
        words.vector.doca <<- c(words.vector.doca, x[paste("TFIDF")])  
      })
      lapply(dtm[[docb]], FUN=function(x){
        words.vector.docb <<- c(words.vector.docb, x[paste("TFIDF")])  
      })
      
      #sum TF-IDF Ai*Bi values
      numerator <- sum(same.words.vector.doca*same.words.vector.docb)
      
      #new denominator
      new.denominator <- sqrt(sum(words.vector.doca^2)*(sum(words.vector.docb^2)))
      
      new.sim.doca.docb <- numerator/new.denominator
      
      #if(doca==doc){
        sim[[paste(doca)]][[paste(docb)]][paste("sim")] <<- new.sim.doca.docb 
      #}else{
        sim[[paste(docb)]][[paste(doca)]][paste("sim")] <<- new.sim.doca.docb
      #}
      
      
    }
  
}



##################################################################################
################# END - Incremental SIM. MATRIX Function #########################
##################################################################################

##################################################################################
####################### Incremental TF-IDF Function ##############################
##################################################################################


UpdateIDF <- function( term, new.doc.flag){
  
  
  N.docs <- length(dtm)
  n.docs <- as.numeric(term.counter.per.doc[paste(term)])
  
  
  
  if(new.doc.flag == TRUE){
    #Add 1 to the counting of documents the term appears
    n.docs <- n.docs + 1
    
    #doc.number <<- 1
    #for every doc with term update IDF
    
    for(doc in 1:length(dtm)){
      
      ln <- names(dtm[[doc]])
      
      if(length(term.index <- which(fmatch(ln,term)==1))==0)
      {
        next
      }
      else
      {
        dtm[[doc]][[term.index]]["IDF"] <<- log(1+(N.docs/n.docs))
      }
      
    }
    
    
  }else{
    
    #return smooth idf
    return(log(1+(N.docs/n.docs)))
    
    
  }
  
}

UpdateTFIDF <- function(doc,term){
  
  #Call UpdateTF
  tf <- dtm[[paste(doc)]][[paste(term)]][paste("TF")]
  
  #Call UpdateIDF
  idf <- dtm[[paste(doc)]][[paste(term)]][paste("IDF")]
  
  #Calculate value
  
  tfidf.value <-  tf*idf
  
  dtm[[paste(doc)]][[paste(term)]][paste("TFIDF")] <<- tfidf.value
  
 
}

PreProcess <- function(text){
  
  #remove stop-words
  #remove punctuation
  #remove numbers
  #change to lower case 
  
  
 
  mycorpus <- Corpus(VectorSource(text))#, readerControl = list(reader = myReader))
  
  
  mydtm <- DocumentTermMatrix(mycorpus,control = list(weighting = function(x)
    weightTf(x),
    stemming = FALSE,
    stopwords = TRUE,
    minWordLength=3,
    stripWhitespace=TRUE,
    content_transformer(tolower),
    removePunctuation = TRUE,
    removeNumbers = TRUE))
  
  mydtm.matrix <- as.matrix(mydtm)
  
  
  
  return(mydtm.matrix)
  
  
}



#######################################################################################
####################### Incremental TF-IDF Function - END #############################
#######################################################################################

dtm <- list()
sim <- list()
term.counter.per.doc <- c()
vec.res <- c()
df.results <- data.frame(row.names =  c("Proc Time (s)") )
graph <- make_empty_graph()
number.of.nodes <- 0
m.neigh.doc <<- matrix(nrow=0,ncol=2,byrow=TRUE)

go <- function(){
  
  time.main2 <- c()
  print(Sys.time())
  
  #Example with 22 chunks of text in the stream
  #Each chunk has several Documents labeled with an ID
  for(i in 1:22){
    
    
    #Read the .csv files
    if(i ==1){
      FILE = read.csv("<start file with some docs labeled with id column>.csv", header=TRUE,sep=",")
    }else{
      FILE = read.csv(paste("<more file names with chunks of data to simulate the stream of data>",i,".csv",sep=""), header=TRUE,sep=",")
    }  
    
    #DOC id in this example is on col2 column of csv data
    docs <- unique(FILE$col2)
    
    tfidf.time.measured <<- 0
    
    #####Start time counting at 1
    tfidf.tmp <- proc.time()
    
    #for each author/doc in stream
    for(doc in docs){
      
      
     #Gets all text from a particular Document author  
     text <- paste(FILE[which(FILE$col2 == doc),"col4"], collapse = " ")
      
     #pre-process text for the author
     tf.matrix <- PreProcess(text)
    
   
    
    #Initializes list of docs with new doc
    if(length(dtm[[paste(doc)]])==0){
      dtm[[paste(doc)]] <<- list()
      sim[[paste(doc)]] <<- list()
      ConstructTextGraph(paste(doc),paste(doc),word=0)
    }
    
    tn <- names(term.counter.per.doc)
    
    for (term in  colnames(unlist(tf.matrix))){
      
      
      #register which new or updated tokens are in the stream
      if(length(term.counter.per.doc[term.index <- which(fmatch(tn,term)==1)]) != 0){
        
        term.counter.per.doc[term.index] <<- term.counter.per.doc[term.index] + 1 
        #Update dtm list
        #Call UpdateTF
        dtm[[paste(doc)]][[paste(term)]][paste("TF")] <<- 1+log(tf.matrix[,paste(term)])
        
        #Construct Graph
        ConstructTextGraph(paste(doc),node = term, word=1)
        PreUpdateSim(paste(term),paste(doc))
        #Call UpdateIDF
        UpdateIDF(paste(term),new.doc.flag = TRUE)
        #Update Incremental TF-IDF
        UpdateTFIDF(paste(doc),term)
      }else{
        new.term <- c(1)
        names(new.term) <- paste(term)
        term.counter.per.doc <<- c(term.counter.per.doc, new.term)
        dtm[[paste(doc)]][[paste(term)]][paste("TF")] <<- 1+log(tf.matrix[,paste(term)])
        
        #Construct Graph
        ConstructTextGraph(paste(doc),node = term, word=1)
        
        idf.value <- UpdateIDF(paste(term),new.doc.flag = FALSE)
        dtm[[paste(doc)]][[paste(term)]][paste("IDF")] <<- idf.value
        #Update Incremental TF-IDF
        UpdateTFIDF(paste(doc),term)
      }  
    }
    
    if(i > 1){
      UpdateSim(m.neigh.doc)
      #reset structure
      m.neigh.doc <<- matrix(nrow=0,ncol=2,byrow=TRUE)
    }
    
   
    
    
  }
  
    tfidf.time.measured <<- proc.time() - tfidf.tmp
    
    #Calls UpdateCache
    vec.res <<- c(vec.res, as.numeric(tfidf.time.measured[3]))
    names(vec.res) <<- paste("iter ",1:i,sep = "")
  
  }
  
  print(Sys.time())
  
  df.results <<- as.data.frame(vec.res, stringsAsFactors = FALSE)
  vars <- c("df.results")
  save(list=vars, file = "<file name>.RData")
  
  #reached the end of the stream
  
  
}

system.time(go())

#Output dtm list of TF-IDF
dtm

#Output sim list of similarities
sim

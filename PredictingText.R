## ----setup, include=FALSE--------------------------------------------------------------------
knitr::opts_knit$set(root.dir = "C:/Users/amanafpour/Desktop/final/en_US")
knitr::opts_chunk$set(echo = FALSE, cache = TRUE)
knitr::knit_hooks$set(inline = function(x) {
  if(!is.numeric(x)){ x }else{ prettyNum(round(x,2), big.mark=",") } 
  })
library(readr)
library(ggplot2)
library(tm)
library(quanteda)
library(hunspell)
library(tokenizers)
library(cld3)
library(data.table)
library(knitr)
library(gridExtra)

## ----data-import-funcs, echo=TRUE------------------------------------------------------------
# Function used to load text files in current directory as data.table
getImportDt <- function(linesnmax = 50000, nsamp = 2000, txtfilename = NULL, skip = 0){
  # linesnmax and nsamp can take on NULL values which mean they are infinity in this function
  # If txtfilename is not defined, all files in current directory are read
  if (is.null(txtfilename)) {
    filenames <- list.files(".") #read filenames from directory
  } else {
    filenames <- txtfilename
  }
  dtread <- NULL
  dtread <- data.table(file = character(), lines = character())
  # Only read the first linesnmax lines of each text doc
  for (filename in filenames) {
    if (is.null(linesnmax)){
      lins <- readr::read_lines(filename) 
    } else {
      lins <- readr::read_lines(filename, n_max = linesnmax, skip = )
    }
    newdt <- data.table(file = filename, lines = lins)
    dtread <- rbind(dtread, newdt)
  }
  dtread$file<-as.factor(dtread$file)
  # Take sample of nsamp size from the read data 
  if (!is.null(nsamp)) dtread <- dtread[, .SD[sample(.N, nsamp)], by = file]
  return(dtread)
}

# This function converts the data table generated from getImportDt to a quanteda corpus
dtToQcorp <- function(dtinput){
  docs <- NULL
  for (filename in levels(dtinput$file)){
    subdt <- dtinput[file == filename]
    doc_str <- paste(subdt$lines, collapse = "\n")
    tempdoc <- corpus(doc_str, docnames = filename)
    if (is.null(docs)) {docs <- tempdoc
    } else {docs <- c(docs, tempdoc)}
  }
  return(docs)
}


## ----tokenization-funcs, echo=TRUE-----------------------------------------------------------
# This function generates a document-feature matrix based on a quanteda corpus
# To avoid tokenizing words from different sentences, first all lines are tokenized into sentences
Get_qdfm <- function(qcorp, n = 1, removeStpwrds = F, addtlStpWrds = NULL){
  sens <- unlist(lapply(tokenize_lines(as.String(qcorp)), tokenize_sentences),
                 use.names = F)
  # Generate a list of explicit swear words
  stpwrds <- readLines("http://www.bannedwordlist.com/lists/swearWords.txt", warn = F)
  if (removeStpwrds == T) stpwrds <- c(stpwrds, stopwords(), addtlStpWrds)
  # The tokenize_ngrams function automatically removes punct and extra whitespace
  ngrams <- unlist(tokenize_ngrams(sens, n=n, lowercase = T, stopwords = stpwrds))
  ngrams <- ngrams[!is.na(ngrams)]
  dfm(as.tokens(list(ngrams)))
}

# This function generates an ordered term-frequency data table based on a quanteda dfm
# The function combines all docs in dfm, info about individual docs is lost
getOrderedFreqDt <- function(dfminput, spellCheck = T){
  dt <- data.table(convert(dfminput, "data.frame"))
  dt <- dt[,-c(1)] # remove first "document" column
  freqv <- colSums(dt)
  freqdt <- data.table(term = names(freqv), freq = freqv)
  if (spellCheck == T) {
    freqdt[, wrongTerms := hunspell(term)]
    freqdt[, correctSpell := identical(wrongTerms[[1]], character(0)), by= 1:nrow(freqdt)]
    freqdt <- freqdt[correctSpell==T, c("term", "freq")]
  }
  freqdt <- freqdt[!grepl("[0-9]", term)] #removes all numbers from ngrams
  setorder(freqdt, -freq)
  return(freqdt)
}

# Function for plotting terms vs occurrences
plot_occurrences <- function(freqDt, nwords = 20){
  wf <- data.frame(term = freqDt$term,
                   occurrences = freqDt$freq,
                   row.names = 1:nrow(freqDt))[1:nwords,]
  ggplot(wf, aes(term, occurrences)) +
    geom_bar(stat="identity") +
    coord_flip() +
    scale_x_discrete(limits=wf$term)
}


## ----predict-model-funcs, echo=TRUE----------------------------------------------------------
# This function splits sentence by their last n number of words and is used in following funcs
getSplitSent <- function(sen, nwrds){
  wrds <- unlist(strsplit(sen, split = " "))
  lastNwrds <- paste(tail(wrds,nwrds), collapse = " ")
  remaining <- paste(head(wrds,length(wrds)-nwrds), collapse = " ")
  c(remaining, lastNwrds)
}

# This function generates a list of frequency matrices up to an maxngram
getFreqMatrix <- function(inputDocs, maxngram = 3, coverage = 1.0) {
  freqList <- list()
  # Generate unigram separately first
  dt <- getOrderedFreqDt(Get_qdfm(inputDocs, n=1), spellCheck = T)
  if (coverage > 0.0 & coverage < 1.0) {
    dt[,cumSumFreq:=cumsum(freq)]
    # Select all extra words not required for specified coverage level
    addtlStopWords <- dt[cumSumFreq > coverage*sum(freq)]$term
    dt <- dt[sumFreq < coverage*sum(freq)]
    dt <- dt[,-c("sumFreq")]
  }
  freqList[[1]] <- dt
  
  # Convert frequency vectors into data tables and split words in ngram cols
  for (num in 2:maxngram){
    dt <- getOrderedFreqDt(Get_qdfm(inputDocs, n=num, addtlStpWrds = addtlStopWords),
                           spellCheck = T)
    dt[,remainingTerm := getSplitSent(term, 1)[1], by= 1:nrow(dt)]
    dt[,lastWrd := getSplitSent(term, 1)[2], by= 1:nrow(dt)]
    dt[,rTermFreq := sum(freq), by = .(remainingTerm)]
    freqList[[num]] <- dt
  }
  return(freqList)
}

getProbMatrix <- function(inputFreqMat) {
  pmat <- list()
  pmat[[1]] <- inputFreqMat[[1]][, p:= freq/sum(freq)]
  for (n in 2:length(inputFreqMat)){
    # Calculate probabilities
    dt <- inputFreqMat[[n]]
    dt[,p := freq / rTermFreq]
    pmat[[n]] <- dt
  }
  return(pmat)
}

# This function cleans and reformats an input sentences used to clean user input before predicting
Clean_Str <- function(inputstr, removeStpwrds = F){
  corpus <- VCorpus(VectorSource(inputstr),
                       readerControl = list(reader=readPlain, language = "en"))
  # Lowercase
  corpus <- tm_map(corpus, content_transformer(tolower))
  # Remove numbers
  corpus <- tm_map(corpus, removeNumbers)
  # Remove explicitly profane words
  profanity <- readLines("http://www.bannedwordlist.com/lists/swearWords.txt", warn = F)
  corpus <- tm_map(corpus, removeWords, profanity)
  # Remove extra whitespace BUT maintain \n line breaks
  whitespaceFUN <- content_transformer(function(x) gsub("[ ]+", " ",as.String(x)))
  corpus <- tm_map(corpus, whitespaceFUN)
  # Remove stop words if applicable
  if (removeStpwrds == T) {
    corpus <- tm_map(corpus, removeWords, words = stopwords("en"))
    }
  return(corpus[[1]]$content)
}

# This function suggest the most probable words to the user based on their input phrase
predictNxtWrd <- function(inputpmat, inputsent) {
  # determine starting n based on length of input sentence and maxngram in prob matrix
  predictTop <- character()
  numWrds <- length(strsplit(inputsent, split = " ")[[1]])
  maxngram <- length(inputpmat)
  n <- numWrds + 1
  if (n > maxngram) n <- maxngram
  # Use "Backoff" to determine suggested words
  for (i in n:2){
    lastNWrds_str <- getSplitSent(inputsent, i - 1)[2]
    # format input sentence
    lastNWrds_str <- Clean_Str(lastNWrds_str)
    # select corresponding ngram matrix
    pdt <- inputpmat[[i]]
    subpdt <- pdt[remainingTerm == lastNWrds_str]
    setorderv(subpdt, c("p"))
    predictTop <- c(predictTop,subpdt$lastWrd)
    # if 3 or less matches found then use next ngram down, otherwise break loop
    if (sum(is.na(predictTop[1:3])) == 0) break
  }
  # If still suggesting less than 3 words, fill in with unigram predictions
  if (sum(is.na(predictTop[1:3])) != 0) {
    pdt <- inputpmat[1]
    setorderv(pdt, c("p"))
    predictTop <- c(predictTop,pdt$term[1:3])
  }
  return(predictTop[1:3])
}


## ----model-eval-funcs------------------------------------------------------------------------
# This function uses the stupid back off method to calculate score for each ngram
getScoreMatrix <- function(trainProbMat, testDocs) {
  nmax <- length(trainProbMat)
  
  # Break testDocs into nmax size grams
  testFreqMatrix <- getFreqMatrix(inputDocs = testDocs, maxngram = nmax, coverage = 1.0)
  testScoreMat <- list()
  
  # For unigram, do +1 laplace smoothing for unseen words
  testScoreDt <- testFreqMatrix[[1]]
  testNRow <- nrow(testScoreDt)
  trainProbDt <- trainProbMat[[1]]
  trainSumFreq <- trainProbDt[,sum(freq)]
  testScoreDt[trainProbDt, s := (i.freq + 1)/ (trainSumFreq + testNRow), on = .(term)]
  testScoreDt[is.na(s), s := (1)/ (trainSumFreq + testNRow)]
  testScoreMat[[1]] <- testScoreDt
  
  # Loop through remaining ngrams and update scores with available training probs
  for (n in 2:nmax) {
    testScoreDt <- testFreqMatrix[[n]]
    trainProbDt <- trainProbMat[[n]]
    testScoreDt[trainProbDt, s:= p, on = .(term)]
    testScoreMat[[n]] <- testScoreDt
  }
  
  # Populate remaining NA score values by backing off to lower ngram probs
  fctr <- .4
  for (n in nmax:2) {
    testScoreDt <- testScoreMat[[n]]
    # Split last nwords from term for next n-gram down
    for (nextn in (n-1):1) {
      npower <- n - nextn
      testScoreDt[, lastnwords := getSplitSent(term, nextn)[2],
                  by = 1:nrow(testScoreDt)]
      testScoreDtNextn <- testScoreMat[[nextn]]
      # Duplicate term col as lastnwords so it can be joined
      testScoreDtNextn[, lastnwords := term]
      testScoreDt[testScoreDtNextn,
                  s := ifelse(is.na(s), fctr ^ npower * i.s, s),
                  on = "lastnwords"]
      testScoreMat[[n]] <- testScoreDt
    }
  }
  return(testScoreMat)
}

calcPerplex <- function(testScoreMat) {
  nmax <- length(testScoreMat)
  perplex <- numeric()
  for(n in 1:nmax) {
    testScoreDt <- testScoreMat[[n]]
    freqsum <- testScoreDt[, sum(freq)]
    testScoreDt[, slog := log(s)]
    slogsums <- testScoreDt[, sum(slog * freq)]
    perplex[n] <- exp(slogsums * (-1/freqsum))
  }
  return(perplex)
}

#TODO: complete function below
calcAccuracy <- function(trainProbMat, testDocs) {
  nmax <- length(testScoreMat)
  # Predict word for each ngram in scores matrix and compare to actual last word

}



memory.limit(100000)

## ----model-testing, eval=FALSE---------------------------------------------------------------
setwd("C:/Users/amanafpour/Desktop/final/en_US")
linesnmax = 50000
nsamp = 1000
set.seed(34341)
trainDt <- getImportDt(linesnmax = linesnmax, nsamp = nsamp)
# testDt <- getImportDt(linesnmax = linesnmax, nsamp = nsamp, skip = linesnmax * 2)
# 
# trainCorp <- dtToQcorp(trainDt)
# testCorp <- dtToQcorp(testDt)
# 
# trainFreqMat <- getFreqMatrix(trainCorp, maxngram = 3, coverage = 1.0)
# testFreqMat <- getFreqMatrix(testCorp, maxngram = 3, coverage = 1.0)
# 
# trainPmat <- getProbMatrix(trainFreqMat)

#testScoresMatrix <- getScoreMatrix(trainPmat, testCorp)
#perplexv <- calcPerplex(testScoresMatrix)

#format(object.size(pmat), units = "auto")
print(predictNxtWrd(trainPmat, "hello there"))

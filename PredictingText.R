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
library(beepr)

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
  stpwrds <- c(stpwrds, addtlStpWrds)
  if (removeStpwrds == T) stpwrds <- c(stpwrds, stopwords())
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
  #TODO: revise to allow appostrophes ' in words 
  freqList <- list()
  addtlStopWords <- NULL
  # Generate unigram separately first
  dt <- getOrderedFreqDt(Get_qdfm(inputDocs, n=1), spellCheck = T)
  if (coverage > 0.0 & coverage < 1.0) {
    dt[,cumSumFreq:=cumsum(freq)]
    # Select all extra words not required for specified coverage level
    addtlStopWords <- dt[cumSumFreq > coverage*sum(freq)]$term
    dt <- dt[cumSumFreq < coverage*sum(freq)]
    dt <- dt[,-c("cumSumFreq")]
  }
  freqList[[1]] <- dt
  
  # Convert frequency vectors into data tables and split words in ngram cols
  for (num in 2:maxngram){
    dt <- getOrderedFreqDt(Get_qdfm(inputDocs, n=num, addtlStpWrds = addtlStopWords),
                           spellCheck = T)
    dt[,c("remainingTerm", "lastWrd") := as.list(getSplitSent(term, 1)), by= 1:nrow(dt)]
    dt[,rTermFreq := sum(freq), by = .(remainingTerm)]
    freqList[[num]] <- dt
  }
  return(freqList)
}

getProbMatrix <- function(inputFreqMat, bareMat = F) {
  pmat <- list()
  pmat[[1]] <- inputFreqMat[[1]][, p:= freq/sum(freq)]
  if (bareMat == T) pmat[[1]] <- pmat[[1]][,c("term", "p")]
  for (n in 2:length(inputFreqMat)){
    # Calculate probabilities
    dt <- inputFreqMat[[n]]
    dt[,p := freq / rTermFreq]
    if (bareMat == T) dt <- dt[, c("remainingTerm", "lastWrd", "p")]
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
predictNxtWrd <- function(inputsent, returnScores = F) {
  # TODO: spell check the input and replace misspelled words with top suggested word before analyzing
  # format input sentence
  inputsent <- Clean_Str(inputsent)
  # determine starting n based on length of input sentence and maxngram in prob matrix
  numWrds <- length(strsplit(inputsent, split = " ")[[1]])
  maxngram <- length(inputpmat)
  n <- numWrds + 1
  if (n > maxngram) n <- maxngram
  # Use "Backoff" to determine suggested words
  predictedWordsDt <- data.table(predictedWord = character(),
                                 score = numeric(),
                                 ngram = numeric())
  fact <- 1.0
  for (i in n:1){
    if (i == 1) {
      # only use unigram if less than 3 words
      if (nrow(predictedWordsDt) < 3) { 
        pdt <- inputpmat[[1]]
        predictedWordsDt <- rbind(predictedWordsDt,
                                  data.table(predictedWord = pdt$term,
                                             score = fact * pdt$p,
                                             ngram = 1))
      } 
      break
    } 
    lastNWrds_str <- paste(tail(unlist(strsplit(inputsent, " ")), i - 1), collapse = " ")
    # select corresponding ngram matrix
    pdt <- inputpmat[[i]]
    subpdt <- pdt[remainingTerm == lastNWrds_str]
    predictedWordsDt <- rbind(predictedWordsDt,
                              data.table(predictedWord = subpdt$lastWrd,
                                         score = fact * subpdt$p,
                                         ngram = i))
    fact <- fact * 0.4 # stupid backoff factor
  }
  # order the data table and return it
  setorderv(predictedWordsDt, cols = c("score"), order = -1, na.last = T)
  # remove duplicates
  predictedWordsDt <- unique(predictedWordsDt, by = c("predictedWord"))
  if (returnScores==F) return(predictedWordsDt$predictedWord[1:3])
  return(predictedWordsDt[1:3])
}

# This function suggest the most probable words to the user based on their input phrase
# predictNxtWrdOLD <- function(inputpmat, inputsent, cleanstring = T) {
#   # TODO: spell check the input and replace misspelled words with top suggested word before analyzing
#   predictTop <- character()
#   # format input sentence
#   if(cleanstring == T) inputsent <- Clean_Str(inputsent)
#   # determine starting n based on length of input sentence and maxngram in prob matrix
#   numWrds <- length(strsplit(inputsent, split = " ")[[1]])
#   maxngram <- length(inputpmat)
#   n <- numWrds + 1
#   if (n > maxngram) n <- maxngram
#   # Use "Backoff" to determine suggested words
#   for (i in n:2){
#     lastNWrds_str <- paste(tail(unlist(strsplit(inputsent, " ")), i - 1), collapse = " ")
#     # select corresponding ngram matrix
#     pdt <- inputpmat[[i]]
#     subpdt <- pdt[remainingTerm == lastNWrds_str]
#     setorder(subpdt, -p)
#     predictTop <- c(predictTop,subpdt$lastWrd[1:3])
#     # if 3 or less matches found then use next ngram down, otherwise break loop
#     if (sum(is.na(predictTop[1:3])) == 0) break
#   }
#   # If still suggesting less than 3 words, fill in with unigram predictions
#   if (sum(is.na(predictTop[1:3])) != 0) {
#     pdt <- inputpmat[[1]]
#     setorder(pdt, -p)
#     predictTop <- c(predictTop,pdt$term[1:3])
#   }
#   return(predictTop[1:3])
# }


## ----model-eval-funcs------------------------------------------------------------------------
# This function uses the stupid back off method to calculate score for each ngram
# getScoreMatrix <- function(trainProbMat, testDocs, smoothing = TRUE) {
#   nmax <- length(trainProbMat)
#   
#   # Break testDocs into nmax size grams
#   testFreqMatrix <- getFreqMatrix(inputDocs = testDocs, maxngram = nmax, coverage = 1.0)
#   testScoreMat <- list()
#   
#   # For unigram, do +1 laplace smoothing for unseen words
#   testScoreDt <- testFreqMatrix[[1]]
#   testNRow <- nrow(testScoreDt)
#   trainProbDt <- trainProbMat[[1]]
#   trainSumFreq <- trainProbDt[,sum(freq)]
#   if (smoothing == T) {
#     testScoreDt[trainProbDt, s := (i.freq + 1)/ (trainSumFreq + testNRow), on = .(term)]
#     testScoreDt[is.na(s), s := (1)/ (trainSumFreq + testNRow)]
#     
#   } else {
#     testScoreDt[trainProbDt, s := p, on = .(term)]
#     testScoreDt[is.na(s), s := 0]
#   }
#   testScoreMat[[1]] <- testScoreDt
#   
#   # Loop through remaining ngrams and update scores with available training probs
#   for (n in 2:nmax) {
#     testScoreDt <- testFreqMatrix[[n]]
#     trainProbDt <- trainProbMat[[n]]
#     testScoreDt[trainProbDt, s:= p, on = .(term)]
#     testScoreMat[[n]] <- testScoreDt
#   }
#   
#   # Populate remaining NA score values by backing off to lower ngram probs
#   fctr <- .4
#   for (n in nmax:2) {
#     testScoreDt <- testScoreMat[[n]]
#     # Split last nwords from term for next n-gram down
#     for (nextn in (n-1):1) {
#       npower <- n - nextn
#       testScoreDt[, lastnwords := paste(tail(unlist(strsplit(term, " ")), nextn),
#                                         collapse = " "),
#                   by = 1:nrow(testScoreDt)]
#       testScoreDtNextn <- testScoreMat[[nextn]]
#       # Duplicate term col as lastnwords so it can be joined
#       testScoreDtNextn[, lastnwords := term]
#       testScoreDt[testScoreDtNextn,
#                   s := ifelse(is.na(s), fctr ^ npower * i.s, s),
#                   on = "lastnwords"]
#       testScoreDt[,lastnwords := NULL]
#       testScoreMat[[n]] <- testScoreDt
#     }
#   }
#   return(testScoreMat)
# }
# 
# calcPerplex <- function(testScoreMat) {
#   nmax <- length(testScoreMat)
#   perplex <- numeric()
#   for(n in 1:nmax) {
#     testScoreDt <- testScoreMat[[n]]
#     freqsum <- testScoreDt[, sum(freq)]
#     slogsums <- testScoreDt[, sum(log(s) * freq)]
#     perplex[n] <- exp(slogsums * (-1/freqsum))
#   }
#   return(perplex)
# }
# 
# getAccMat <- function(trainProbMat, testScoreMat) {
#   accMat <- list()
#   nmax <- length(testScoreMat)
#   # Predict word for each ngram in scores matrix and compare to actual last word
#   for (n in 2:nmax) {
#     testScoreDt <- testScoreMat[[n]]
#     # Remove variables not required
#     testScoreDt[,c("term", "rTermFreq", "s") := NULL]
#     # Try a non-loop way of updating 1st word correct
#     testScoreDt[, c("correctwrd1","correctwrd2", "correctwrd3") :=
#                   as.list(predictNxtWrd(trainProbMat, remainingTerm, cleanstring = F) == lastWrd),
#                 by = 1:nrow(testScoreDt)]
#     testScoreDt[is.na(correctwrd1), correctwrd1 := FALSE]
#     testScoreDt[is.na(correctwrd2), correctwrd2 := FALSE]
#     testScoreDt[is.na(correctwrd3), correctwrd3 := FALSE]
#     # return vector of total correct predictions
#     accVec <- c(testScoreDt[,sum(correctwrd1*freq)/sum(freq)],
#                 testScoreDt[,sum(correctwrd2*freq)/sum(freq)],
#                 testScoreDt[,sum(correctwrd3*freq)/sum(freq)])
#     names(accVec) <- c("acc1stwrd", "acc2ndwrd", "acc3rdwrd")
#     accMat[[n]] <- accVec
#   }
#   return(accMat)
# }



memory.limit(100000)

## ----model-testing, eval=FALSE---------------------------------------------------------------
setwd("C:/Users/amanafpour/Documents/GitHub/PredictingText-CapstoneProject/Saved Data")
profanity <<- readLines("http://www.bannedwordlist.com/lists/swearWords.txt", warn = F)

# setwd("C:/Users/amanafpour/Desktop/final/en_US")
# linesnmax = 500000
# trainnsamp = 100000
# testnsamp = 10000
# cover = 1.0
# maximumngrams = 5
# 
# set.seed(34341)
# trainDt <- getImportDt(linesnmax = linesnmax, nsamp = trainnsamp)
# trainCorp <- dtToQcorp(trainDt)
# 
# trainFreqMat <- getFreqMatrix(trainCorp, maxngram = maximumngrams, coverage = cover)
# trainPmat <- getProbMatrix(trainFreqMat)
# setwd("C:/Users/amanafpour/Documents/GitHub/PredictingText-CapstoneProject/Saved Data")
# saveRDS(trainPmat, file = paste("trainPmat", trainnsamp, "lines_", maximumngrams, "gram_", cover, "cover.RDS", sep=""))
# beep()

# set.seed(789)
# testDt <- getImportDt(linesnmax = linesnmax, nsamp = testnsamp, skip = linesnmax)
# setwd("C:/Users/amanafpour/Documents/GitHub/PredictingText-CapstoneProject/Saved Data")
# saveRDS(testDt, file = paste("testDt", testnsamp, "lines.RDS", sep=""))




# Read input mat from directory
inputpmat_orig <- readRDS(paste("trainPmat1e+05lines_5gram_1cover.RDS",sep = ""))

newpmat <<- inputpmat_orig
# Only leave first 3 lines for 1-gram
newpmat[[1]] <- newpmat[[1]][1:3]
# Subselect terms with freq > 1 and only selected required columns
for (n in 2:5) {
  newpmat[[n]] <- newpmat[[n]][freq>1, c("remainingTerm", "p", "lastWrd")]
}

inputpmat <<- newpmat

# Run benchmark on all ngrams
# bmlist <- list()
# for (n in 2:5){
#   inputpmat <<- inputpmat_orig[2:n]
#   bmlist[n] <- benchmark(predictNxtWrd,
#                          sent.list = list('tweets' = tweets, 
#                                           'blogs' = blogs))
# }
# 
# saveRDS(bmlist, file = paste("bmresults_1000testlines_","10000lines_5gram",".RDS", sep= ""))

# testDt <- readRDS("testDt1000lines.RDS")
# testCorp <- dtToQcorp(testDt)
# 
# trainPmatFilenames <- c("10000lines_5gram_0.5cover",
#                         "10000lines_5gram_0.8cover",
#                         "10000lines_5gram_1cover")
#  
# # Test out this section below
# accMatList <- list()
# for (filesuffix in trainPmatFilenames){
#   trainPmat <- readRDS(paste("trainPmat",filesuffix,".RDS",sep = ""))
#   testScoreMatrix <- getScoreMatrix(trainPmat, testCorp, smoothing = F)
#   timerun <- system.time(accMat <- getAccMat(trainPmat, testScoreMatrix))
#   accMatList[[filesuffix]] <- accMat
# }
# 
# saveRDS(accMatList, file = paste("accMatList_1000testlines_","10000lines_5gram",".RDS", sep= ""))

# perplexvList <- list()
# for (filesuffix in trainPmatFilenames){
#   trainPmat <- readRDS(paste("trainPmat",filesuffix,".RDS",sep = ""))
#   testScoreMat_smoothed <- getScoreMatrix(trainPmat, testCorp, smoothing = T)
#   perplexv <- calcPerplex(testScoreMat_smoothed)
#   perplexvList[[filesuffix]] <- perplexv
# }

#format(object.size(pmat), units = "auto")
#print(predictNxtWrd(trainPmat, "hello there"))

# Rprof(tmp <- tempfile(), interval = .02)
# accMat <- getAccMat(testScoreMat)
# Rprof(NULL)
# summaryRprof(tmp)
# unlink(tmp)

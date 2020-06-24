library(shiny)
library(DT)
library(readr)
library(tm)
library(quanteda)
library(hunspell)
library(tokenizers)
library(data.table)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
    
    output$predResults = DT::renderDataTable({
        DT::datatable(predictNxtWrd(input$sentInput, T), options = list(dom = "t",
                                                                        autoWidth = T,
                                                                        columnDefs = list(list(width = "100px"))))
    })
    
    inputpmat <<- readRDS(paste("finalPmat1e+05lines_3gram_1freqRemoved.RDS",sep = ""))
    profanity <<- readLines("http://www.bannedwordlist.com/lists/swearWords.txt", warn = F)
    
    # This function splits sentence by their last n number of words and is used in following funcs
    getSplitSent <- function(sen, nwrds){
        wrds <- unlist(strsplit(sen, split = " "))
        lastNwrds <- paste(tail(wrds,nwrds), collapse = " ")
        remaining <- paste(head(wrds,length(wrds)-nwrds), collapse = " ")
        c(remaining, lastNwrds)
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
})

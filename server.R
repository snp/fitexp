
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
# 
# http://www.rstudio.com/shiny/
#

library(shiny)
options(shiny.maxRequestSize=300*1024^2)
median.fixed <- function(data){
  median(data[is.finite(data)],na.rm=TRUE)
}
median.ratio <- function(data, numerator, denominator){
  # median.fixed <- median
  numerator.columns <- numerator
  denominator.columns <- denominator
  res <- data.frame(row.names=rownames(data))
  message("Numerator: ", paste(numerator.columns,collapse='; '))
  message("Denominator: ",paste(denominator.columns,collapse='; '))
  for(nc in numerator.columns){
    for (dc in denominator.columns){
      #print(paste(nc,dc))
      res <- cbind(res,data[,nc]/data[,dc])
    }
  }
  apply(res,1,median.fixed)
}
getTargets <- function(data, cell_rx, drug_rx, control_rx){
  require(dplyr)
  data <- tbl_df(data)
  print(nrow(data))
  # LFQ data from MaxQuant
  if(!("Reverse" %in% names(data)))
    data$Reverse <- ''
  if(!("Potential.contaminant" %in% names(data)))
    data$Potential.contaminant <- ''
  
  data %>%
    filter(Reverse != '+') %>%
    filter(Potential.contaminant!='+') %>% 
    as.data.frame(.)-> dd
  dd %>%
    select(Protein.IDs, Protein.names, Gene.names, Peptides) -> data.info
  dd %>%
    select(starts_with('LFQ')) %>%
    log2() -> lfq.data
  # Normalization
  lfq.data <- as.data.frame(apply(lfq.data, 2, function(x){x - median(x, na.rm=T)}))
  cells <- as.character((sub(cell_rx,'\\1',names(lfq.data))))
  drugs <- as.character((sub(drug_rx,'\\1',names(lfq.data))))
  ldata <- as.data.frame(2^lfq.data)
  rdata <- data.frame(row.names=rownames(ldata))
  print(unique(cells))
  print(unique(drugs))
  for(cell in unique(cells)){
    for(drug in unique(drugs)){
      target.cols <- names(ldata)[(cells==cell)&(drugs==drug)]
      others.cols <- names(ldata)[(cells==cell)&grepl(control_rx, names(lfq.data))]
      res <- median.ratio(
        ldata,
        target.cols,
        others.cols
      )
      if(sum(!is.na(res))>0){
        rdata[,sprintf("%s_%s.regulation", drug, cell)] <- res
      }
    }
  }
  # Specificity
  for(cell in unique(cells)){
    for(drug in unique(drugs)){
      target.cols <- names(ldata)[(cells==cell)&(drugs==drug)]
      others.cols <-names(ldata)[(cells==cell)&(!grepl(control_rx, names(lfq.data)))&(drugs != drug)]
      res <- median.ratio(
        ldata,
        target.cols,
        others.cols
      )
      if(sum(!is.na(res))>0){
        rdata[,sprintf("%s_%s.specificity", drug, cell)] <- res
      }
    }
  }
  result <- list()
  rrdata <- rdata
  for(drug in unique(drugs[!grepl(control_rx, names(lfq.data))])){
    drug.res <- data.frame(row.names=rownames(rrdata))
    tRPMAData<-log2(rrdata[,grepl(drug,names(rrdata))])
    if(is.data.frame(tRPMAData)){
      NumElements<-rowSums(!is.na(tRPMAData))
      NumReps <- max(unique(NumElements))
      RPMAownUp_pvalues<-RPMAownDown_pvalues<-NULL
      for (d in unique(NumElements)) {
        RPMAData<-tRPMAData[NumElements==d,]
        if(d>1 && length(as.matrix(RPMAData))>ncol(tRPMAData)) {
          RP.own<-0
          Rank<-NULL
          RankNAs<-0
          for (r in 1:NumReps) {
            Rank[[r]]<-rank(RPMAData[,r],na.last="keep")/(sum(!is.na(RPMAData[,r]))+1)
            names(Rank[[r]]) <-rownames(RPMAData)
            Rank[[r]][is.na(Rank[[r]])]<-1
            RP.own<-RP.own+log(Rank[[r]])
            RankNAs<-RankNAs+sum(Rank[[r]]>1)
          }
          RP.own<-exp(RP.own)
          RPownCorr<- -log(RP.own)
          RPMAownUp_pvalues<-c(RPMAownUp_pvalues,pgamma(RPownCorr,d))
          
          RP.own<-0
          Rank<-NULL
          RankNAs<-0
          for (r in 1:NumReps) {
            Rank[[r]]<-rank(-RPMAData[,r],na.last="keep")/(sum(!is.na(RPMAData[,r]))+1)
            names(Rank[[r]]) <-rownames(RPMAData)
            Rank[[r]][is.na(Rank[[r]])]<-1
            RP.own<-RP.own+log(Rank[[r]])
            RankNAs<-RankNAs+sum(Rank[[r]]>1)
          }
          RP.own<-exp(RP.own)
          RPownCorr<- -log(RP.own)
          RPMAownDown_pvalues<-c(RPMAownDown_pvalues,pgamma(RPownCorr,d))
        }
      }
      RPMAown_pvalues<-2*apply(cbind(RPMAownDown_pvalues,RPMAownUp_pvalues),1,min)
      # RPMAown_pvalues<-RPMAown_pvalues[RPMAown_pvalues<1]
      drug.res <- data.frame(
        row.names=names(RPMAown_pvalues),
        p.value=RPMAown_pvalues,
        updown=ifelse(RPMAownDown_pvalues>RPMAownUp_pvalues, 'Up','Down'))
      #    d2 <- merge(d2, drug.res, by=0, all=F, suffixes=c('',drug))
      
      drug.res <- merge(drug.res, data.info, by=0, all=F)
      drug.res$p.value[drug.res$p.value>1] <- 1
      drug.res$p.adjust <- p.adjust(drug.res$p.value)
      drug.res <- drug.res[drug.res$p.value>0,]
      
     drug.res %>%
       arrange(p.value) %>%
       select(Protein.IDs,Protein.names, Gene.names, Peptides, updown, p.value, p.adjust) -> result[[drug]]
    }
  }
  result
}

shinyServer(function(input, output) {
   
  output$contents <- renderDataTable({
    
    inFile <- input$file1
    
    if (is.null(inFile))
      return(NULL)
    data <- read.csv(inFile$datapath, head=T, sep=input$sep, comment='', na.strings=c('NA','#N/A','0.0','0'))
    rp.targets <- getTargets(data, input$cell_r, input$drug_r, input$ctrl)
    res <- data.frame()
    for(drug in names(rp.targets)){
      dframe <- rp.targets[[drug]]
      dframe$Drug <- drug
      res <- rbind(res,select(dframe, Drug, Protein.names, Gene.names, Peptides, updown, p.value, p.adjust))
    }
    as.data.frame(res)
#    print(nrow(res))
  })
})

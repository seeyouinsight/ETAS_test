setwd("C:/Users/test/Desktop/ETAS/Hawkes/how_to_use_Hawkes/")
#{r setup, include=FALSE}
if(.Platform$OS.type == "windows" ){ n.cores = 1}
knitr::opts_chunk$set(echo = TRUE, warnings = FALSE, message = FALSE)
# load the functions needed for this analysis
source('hawkes_functions.R')
list.input <- input.file.to.list('user_input_amatrice.txt')
str(list.input)
ETAS.model.fit <- Temporal.ETAS.fit(list.input)

list.output<-as.data.frame(list.output)
post.list <- get_posterior_param(input.list = list.output)
post.list$post.plot
p.samp <- post_sampling(input.list = list.output, n.samp = 1000)

dd.ama <- read.csv2(file = 'data_M3.0.csv', header = TRUE, sep = ',') %>%
  mutate(time_date = as.POSIXct(paste0(year,'-',month,'-',day,' ',hr,':',min,':',sec)),
         time.diff = as.numeric(difftime(time_date, min(time_date) - 1, units = 'days')),
         Mw = as.numeric(Mw),
         Lon = as.numeric(Lon),
         Lat = as.numeric(Lat),
         Mw.class = cut(Mw, breaks = c(M0, 5, 7))) %>%
  arrange(time_date)
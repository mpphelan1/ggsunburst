## State = Event of interest (i.e. treatment)
## Stage = When an event occurred (i.e. first treatment, second, etc...)

## Workflow :
## assuming data looks like:
#  patid start end state
#      1     1   7     A
#      1     8  10     B
#      2     1  10     C
# ..
#  Restructure to:
#  patid stage  state
#      1     1      A
#      1     2      B
#      2     1      C
#  Summarize counts at each stage/state combo w.
#  Plot sunburts

library(dplyr)
library(ggplot2)

make_sequence <- function(df, patid=patid, state=state, stage=stage){
  patid <- enquo(patid)
  state <- enquo(state)
  stg <- enquo(stage)
  
  df %>%
    arrange(!!patid, !!stg, !!state) %>%
    ## Collapse multiple !!states within a single !!stage
    group_by(!!patid, !!stg) %>%
    arrange(!!patid, !!stg, !!state) %>%
    summarize(state = paste(unique(!!state), collapse = "/")) %>%
    ungroup %>%
    group_by(!!patid) %>%
    mutate(stage = !!stg,
      stage_seq = purrr::accumulate(.x=state, 
                                    paste,sep = '->'),
      prv_seq = lag(stage_seq,default = 'Initial' )) %>% 
    group_by(stage, stage_seq, prv_seq) %>%
    summarize(ct = n()) %>%
    group_by(prv_seq) %>%
    mutate(mult1 = sum(ct)) %>%
    ungroup %>%
    {
      ringNs <- group_by(.,stage) %>% 
        summarize(ringN = sum(ct))
      left_join(.,ringNs)
    } %>%
    {
      df1 <- filter(., prv_seq == 'Initial') %>% 
        mutate(ct_prv = sum(ct))
      df2 <- group_by(., stage, prv_seq) %>%
        mutate(ct_prv = ct,
               stage = stage + 1,
               prv_seq = stage_seq) %>%
        select(prv_seq,  stage, ct_prv) 
      df2 <- filter(., prv_seq != 'Initial') %>% left_join(df2)
      bind_rows(df1,df2)
    } %>%
    group_by(stage, prv_seq) %>%
    mutate(pct = ct/ct_prv,
           pct_cum = cumsum(pct),
           pct_cum_min = lag(pct_cum, default = 0),
           pct_cum_med = (pct_cum_min + pct_cum) /2
    ) %>%
    ungroup
  
}
# 
# nPats <- 101L
# mult1 <- 0.9
# mult2 <- 0.8
# set.seed(1)
# set.seed(12345)
# ring1 <- data.frame(
#   patientID = 1:nPats,
#   state = sample(letters[1:3],nPats, replace = T),
#   stage = 1
# )
# n2 <- floor(nPats*mult1)
# ring2 <- data.frame(
#   patientID = sample(ring1$patientID, n2),
#   state = sample(letters[1:4],n2, replace = T),
#   stage = 2
# )
# n3 <- n2*mult2
# ring3 <- data.frame(
#   patientID = sample(ring2$patientID,n3),
#   state = sample(letters[1:4],n3, replace = T),
#   stage = 3
# )
substrRight <- function(x){
  sapply(x, function(xx)
    substr(xx, (nchar(xx)), nchar(xx))
  )
}
# 
# df <- bind_rows(ring1, ring2, ring3)
# 
# paste2 <- function(x, y, sep = ".") paste(x, y, sep = '+')
# max_stg <- max(df$stage)
# 
# df_sun1 <- df %>% make_sequence()
# 
# df_sun1 <- df %>%
#   arrange(patientID, stage, state) %>%
#   group_by(patientID) %>%
#   mutate(
#     stage_seq = purrr::accumulate(.x=state, 
#                                   paste,sep = '->'),
#     prv_seq = lag(stage_seq,default = 'Initial' )) %>% 
#   group_by(stage, stage_seq, prv_seq) %>%
#   summarize(ct = n()) %>%
#   group_by(prv_seq) %>%
#   mutate(mult1 = sum(ct)) %>%
#   ungroup %>%
#   {
#     ringNs <- group_by(.,stage) %>% 
#       summarize(ringN = sum(ct))
#     left_join(.,ringNs)
#   } %>%
#   {
#     df1 <- filter(., prv_seq == 'Initial') %>% 
#       mutate(ct_prv = sum(ct))
#     df2 <- group_by(., stage, prv_seq) %>%
#       mutate(ct_prv = ct,
#              stage = stage + 1,
#              prv_seq = stage_seq) %>%
#       select(prv_seq,  stage, ct_prv) 
#     df2 <- filter(., prv_seq != 'Initial') %>% left_join(df2)
#     bind_rows(df1,df2)
#   } %>%
#   group_by(stage, prv_seq) %>%
#   mutate(pct = ct/ct_prv,
#          pct_cum = cumsum(pct),
#          pct_cum_min = lag(pct_cum, default = 0),
#          pct_cum_med = (pct_cum_min + pct_cum) /2
#   ) %>%
#   ungroup
# 


## Case 1 Overlap States ##
## If there is shared window create a new start/end window when both exist 
## I.e. [0,4] & [3,6] -> [0,3,], [3,4], [4,6]
df_window <- data.frame(patid = 1,
                        start = c(0,3,100),
                        end = c(4,6,101),
                        state = c('A','B','C'))
## window_to_stage 
## A function that transforms data that has a start/end time associated with each state to a data.frame
## capturing the ordering of each window. 
## If there are overlapping windows that overlap is treated as its own stage
## 'gap' is the maximum gap between two states before an empty window is considered
## 'gap_val' is the state that is entered during that empty window

window_to_stage <- function(df, gap=NA,gap_val = 'None'){
  df %>% 
    arrange(patid, start) %>%
    group_by(patid) %>%
    mutate(ld = lead(start),
           lg = lag(end),
           dlta = start - lg) %>% {
             if(!is.na(gap)){
             empty = filter(., dlta > gap) %>%
               mutate(end = start,
                      start = lg,
                      state = gap_val)
             bind_rows(.,empty) %>% arrange(patid,start)}
             else {.}
           } %>% {
             cln  = filter(.,is.na(ld) | ld >= end,
                           is.na(lg) | lg <= start)
             splt_A = filter(.,ld < end) %>% {
               aa = mutate(., end = ld)
               ab = mutate(., start = ld)
               bind_rows(aa,ab)
             }
             splt_B = filter(., lg > start) %>% {
               aa = mutate(., end = lg)
               ab = mutate(., start = lg)
               bind_rows(aa,ab)
             }
             df_out = bind_rows(cln,splt_A) %>% 
               bind_rows(splt_B) %>% 
               select(patid, start,end, state) %>%
               arrange(patid, start, state)
             df_state = df_out %>% select(patid,start) %>% distinct() %>% group_by(patid) %>%
               mutate(stage = row_number()) 
             left_join(df_out,df_state)
           } %>%
    select(patid, state, stage)
}
# window_to_stage(df_window)

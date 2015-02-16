require(plyr)
require(dplyr)
den <- data.frame(
  Rep = rep(1:3,each=100),
  val = runif(300)
  )
d2 <-den %>%  group_by(Rep) %>% mutate(nt=ntile(Rep,10)) %>%  
  group_by(Rep,nt) %>% summarize(meanV=mean(val)) %>% 
  mutate(nt=nt*100)

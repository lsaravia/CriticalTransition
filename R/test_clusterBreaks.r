
      s<-"A B C
      1 2 1
      2 22 3
      0 0 -1
      2 12 2
      0 0 -1
      20 2 5
      1 3 1
      0 2 2"
      d<-read.delim(textConnection(s),sep=" ",header=T)

      require(dplyr)
      
      mutate(d,rn=row_number()) %>% filter(C==-1) 

      bind_rows(slice(d, 1:2) %>% mutate(grp=1),slice(d,4) %>%mutate(grp=2), slice(d,6:n()) %>% mutate(grp=3)) 

# Solution
      
      d %>% mutate(grp = cumsum(C == -1) + 1) %>% filter(C != -1)

    require(quantreg)
    data("engel")
    require(dplyr)
    engel$grp <- trunc(runif(nrow(engel), min=0, max=3))
    
    group_by(engel,grp) %>% do(summary(rq(foodexp~income,data=.,tau=c(.05, .25, .5, .75, .95)),se="boot"))
    
    rqm <- group_by(engel,grp) %>% do(mdl=rq(foodexp~income,data=.,tau=c(.05, .25, .5, .75, .95)))
    summarise(rqm, coef(summary(mdl,se="boot")))

    f <-group_by(engel,grp) %>% do(as.data.frame(lapply(summary(rq(foodexp~income,data=.,tau=c(.05, .25, .5, .75, .95)),se="boot"), coef)))
    
    f <-group_by(engel,grp) %>% 
      do(as.data.frame(do.call(rbind,
                               lapply(summary(rq(foodexp~income,data=.,tau=c(.05, .25, .5, .75, .95)), se="boot"), coef)
      ), row.names = NULL))
    
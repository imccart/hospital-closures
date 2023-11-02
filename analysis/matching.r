## Organize AHA and neighbor data for matching

aha.match <- aha.data %>% 
    left_join(aha.neighbors, by=c("ID", "year")) %>%
    mutate(closed=case_when(
                !is.na(change_type) & change_type=="Closure" ~ 1,
                !is.na(change_type) & change_type!="Closure" ~ 0,
                is.na(change_type) ~ 0))

ever.cah <- aha.match %>% 
    group_by(ID) %>% 
    summarize(ever_cah=sum(critical_access)) %>% 
    mutate(ever_cah=ifelse(ever_cah>0, 1, 0))

## Match on distance and bed size (average over time for now)
aha.id.mean <- aha.match %>% 
    group_by(ID) %>% 
    summarize(mean_bed=mean(BDTOT, na.rm=TRUE),
              mean_dist=mean(min_dist, na.rm=TRUE),
              ever_close=max(closed, na.rm=TRUE)) %>%
    left_join(ever.cah, by="ID") %>%
    filter(!is.na(mean_dist), !is.na(mean_bed))


logit.reg <- glm(ever_cah ~ mean_bed + mean_dist, data=aha.id.mean, family=binomial(link='logit'))
aha.id.mean <- aha.id.mean %>%
  mutate(ps = predict(logit.reg, type = 'response')) %>%
  filter(ps>0 & ps<1) %>%
  mutate(ipw = case_when(
    ever_cah == 1 ~ 1/ps,
    ever_cah == 0 ~ 1/(1-ps)
  ))

reg.ipw <- lm(ever_close ~ ever_cah, data=aha.id.mean, weights=ipw)
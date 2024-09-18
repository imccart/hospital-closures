## Inverse probability weighting ---------------------------------------------

logit.dat <- final.dat %>%
  group_by(ID) %>% 
  summarize(ever_cah=max(ever_cah),
            mean_bed=mean(BDTOT, na.rm=TRUE),
            mean_dist=mean(distance, na.rm=TRUE)) %>%
  filter(!is.na(mean_bed), !is.na(mean_dist))


logit.reg <- glm(ever_cah ~ mean_bed + mean_dist, data=logit.dat, 
          family=binomial(link='logit'))
id.weights <- logit.dat %>%
  mutate(ps = predict(logit.reg, type = 'response')) %>%
  filter(ps>0 & ps<1) %>%
  mutate(ipw = case_when(
    ever_cah == 1 ~ 1/ps,
    ever_cah == 0 ~ 1/(1-ps)
  ))

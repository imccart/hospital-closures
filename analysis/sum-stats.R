# Quick summary -----------------------------------------------------------

## plot closures, mergers, and both over time by treat_state
changes.plot <- est.dat %>%
  group_by(year, treat_state) %>%
  summarize(closures=sum(closed),
            mergers=sum(merged),
            all_changes=sum(closed_merged)) 

closure.plot  <- changes.plot %>%
  ggplot(aes(x=year, y=closures, color=as.factor(treat_state))) +
  geom_line() +
  geom_vline(xintercept=1999, linetype="dashed") +
  scale_color_manual(values=c("black", "red")) +
  theme_bw() +
  labs(x="Year", y="Number of closures", color="Treatment group") +
  theme(legend.position="bottom")
ggsave("results/closures.png", closure.plot, width = 6, height = 10, dpi = 300)


merger.plot  <- changes.plot %>%
  ggplot(aes(x=year, y=mergers, color=as.factor(treat_state))) +
  geom_line() +
  geom_vline(xintercept=1999, linetype="dashed") +
  scale_color_manual(values=c("black", "red")) +
  theme_bw() +
  labs(x="Year", y="Number of mergers", color="Treatment group") +
  theme(legend.position="bottom")
ggsave("results/merger.png", merger.plot, width = 6, height = 10, dpi = 300)


anychange.plot  <- changes.plot %>%
  ggplot(aes(x=year, y=all_changes, color=as.factor(treat_state))) +
  geom_line() +
  geom_vline(xintercept=1999, linetype="dashed") +
  scale_color_manual(values=c("black", "red")) +
  theme_bw() +
  labs(x="Year", y="Number of closures or mergers", color="Treatment group") +
  theme(legend.position="bottom")
ggsave("results/all-changes.png", anychange.plot, width = 6, height = 10, dpi = 300)


## count of critical access by year
view(final.dat %>% group_by(year) %>% summarize(tot_cah=sum(cah)))

## first critical access hospital year by state
view(final.dat %>% filter(cah==1) %>% group_by(MSTATE) %>% summarize(first_cah=min(year)))

## bed sizes of hospitals that closed
bed.closed <- final.dat %>% filter(change_type=="Closure") %>% 
  mutate(bed_bin=cut(BDTOT, breaks=c(0,5,10,15,20,25,30,40,50,75,100,150,200,250,300,400,500,999))) %>% 
  count(bed_bin) %>% arrange(bed_bin)

## CAH designation among hospitals that closed
cah.closed <- final.dat %>% mutate(all_cah=sum(critical_access)) %>% 
  group_by(change_type, critical_access) %>% summarize(n=n(), all_cah=first(all_cah), all_hosp=nrow(aha.combine))

## Count of CAH and Non-CAH over time
fig.hosp.type <- final.dat %>% group_by(year, critical_access) %>%
  summarize(hosp_count=n())  %>%
  ggplot(aes(x=year, y=hosp_count, group=critical_access)) + 
  geom_line() + geom_point() + theme_bw() +
  geom_text(data = final.dat %>% filter(year==2016) %>% group_by(critical_access, year) %>% summarize(hosp_count=n()), 
            aes(label = c("Non-CAH", "CAH"),
                x = year+1,
                y = hosp_count-100)) +
  scale_y_continuous(labels = comma,
                     breaks=seq(0, 8000, 500)) +
  scale_x_continuous(breaks=seq(1980, 2019, 1)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  guides(linetype="none") +
  labs(
    x="Year",
    y="Count of Hospital Types",
    caption="Limited to general, short-term, acute care community hospitals"
  )

## Share of CAH and Non-CAH closures over time
cah.change <- final.dat %>%
  filter(critical_access==1) %>%
  group_by(year, change_type) %>%
  summarize(cah_change=n()) %>%
  group_by(year) %>%
  mutate(all_cah=sum(cah_change),
         prop_cah=cah_change/all_cah) %>%
  select(year, prop_cah, change_type) %>%
  ungroup()

non.cah.change <- final.dat %>%
  filter(critical_access==0) %>%
  group_by(year, change_type) %>%
  summarize(non_cah_change=n()) %>%
  group_by(year) %>%
  mutate(all_non_cah=sum(non_cah_change),
         prop_non_cah=non_cah_change/all_non_cah) %>%
  select(year, prop_non_cah, change_type) %>%
  ungroup()

all.change <- final.dat %>%
  group_by(year, change_type) %>%
  summarize(any_change=n()) %>%
  group_by(year) %>%
  mutate(all_hosp=sum(any_change),
         prop_change=any_change/all_hosp) %>%
  select(year, change_type, prop_change) %>%
  ungroup()


fig.close <- all.change %>%
  left_join(non.cah.change,
            by=c("year","change_type")) %>%
  left_join(cah.change, 
            by=c("year", "change_type")) %>%
  filter(change_type %in% c("Closure", "Merger")) %>%
  select(year, change_type, prop_non_cah, prop_cah) %>%
  replace_na(list(prop_non_cah=0, prop_cah=0)) %>%
  pivot_longer(cols = starts_with("prop"), names_to = "type", values_to = "value") %>%
  ggplot(aes(x = year, y = value, color = type, linetype = change_type)) +
  geom_line() +
  ylim(0,0.05) + 
  scale_x_continuous(breaks=seq(1980, 2019, 1)) +
  labs(x = "Year", y = "Proportion", color = "Hospital Type", linetype = "Change Type") +
  theme_minimal() +
  scale_color_discrete(labels=c("CAH","Non-CAH")) + 
  guides(color = guide_legend(order = 2), linetype = guide_legend(order = 1)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        legend.position=c(0.8,0.65))



## Treatment timing --------------------------------------------------------

panel.view <- final.dat %>%
  mutate(cah_treat=ifelse(year>=first_year_obs & first_year_obs>0, 1, 0)) %>%
  group_by(MSTATE, year) %>%
  summarize(hospitals=n(), cah_treat=max(cah_treat), 
    closures=sum(closed), sum_cah=sum(cah, na.rm=TRUE))

panel.cah <- panelview(1~cah_treat, data=panel.view %>% mutate(MSTATE=as.character(MSTATE)), 
                      index=c("MSTATE","year"), legendOff=TRUE, 
          theme.bw=TRUE, by.timing=TRUE, xlab="Year", ylab="State",
          main="", color=c("white","gray"))
ggsave("results/panelview-cah-treat.png", panel.cah)

panel.closure <- panelview(closures ~ cah_treat,
          data = panel.view, index = c("MSTATE","year"), type = "outcome",
          main = "", by.cohort = TRUE, outcome.type="discrete",
          legend.labs = c("Control States","Treated States (before CAH)",
                          "Treated States (after CAH)"),
          theme.bw=TRUE)
ggsave("results/panelview-cah-closure.png", panel.closure)   
# Quick summary -----------------------------------------------------------

## bed sizes of hospitals that closed
bed.closed <- aha.data %>% filter(change_type=="Closure") %>% 
  mutate(bed_bin=cut(BDTOT, breaks=c(0,5,10,15,20,25,30,40,50,75,100,150,200,250,300,400,500,999))) %>% 
  count(bed_bin) %>% arrange(bed_bin)

## CAH designation among hospitals that closed
cah.closed <- aha.data %>% mutate(all_cah=sum(critical_access)) %>% 
  group_by(change_type, critical_access) %>% summarize(n=n(), all_cah=first(all_cah), all_hosp=nrow(aha.combine))

## Count of CAH and Non-CAH over time
fig.hosp.type <- aha.data %>% group_by(year, critical_access) %>%
  summarize(hosp_count=n())  %>%
  ggplot(aes(x=year, y=hosp_count, group=critical_access)) + 
  geom_line() + geom_point() + theme_bw() +
  geom_text(data = aha.data %>% filter(year==2016) %>% group_by(critical_access, year) %>% summarize(hosp_count=n()), 
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
cah.change <- aha.data %>%
  filter(critical_access==1) %>%
  group_by(year, change_type) %>%
  summarize(cah_change=n()) %>%
  group_by(year) %>%
  mutate(all_cah=sum(cah_change),
         prop_cah=cah_change/all_cah) %>%
  select(year, prop_cah, change_type) %>%
  ungroup()

non.cah.change <- aha.data %>%
  filter(critical_access==0) %>%
  group_by(year, change_type) %>%
  summarize(non_cah_change=n()) %>%
  group_by(year) %>%
  mutate(all_non_cah=sum(non_cah_change),
         prop_non_cah=non_cah_change/all_non_cah) %>%
  select(year, prop_non_cah, change_type) %>%
  ungroup()

all.change <- aha.data %>%
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


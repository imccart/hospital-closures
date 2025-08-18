# Check missing information -----------------------------------------------

non.missing.counts <- aha.data %>%
  group_by(year) %>%
  summarise(across(everything(), ~ sum(!is.na(.) & !(. %in% "")))) %>%
  pivot_longer(cols = -year, names_to = "Variable", values_to = "Count") %>%
  pivot_wider(names_from = year, values_from = Count)

check.distance <- aha.neighbors %>%
  group_by(year) %>%
  summarize(mean_dist=mean(distance, na.rm=TRUE),
            count_dist=sum(!is.na(distance)),
            count_all=n())


# Quick summary -----------------------------------------------------------

## plot closures, mergers, and both over time by treat_state
changes.plot <- est.dat %>%
  group_by(year, MSTATE, treat_state) %>%
  summarize(
    closures     = sum(closed),
    mergers      = sum(merged),
    all_changes  = sum(closed_merged),
    .groups = "drop"
  ) %>%
  group_by(year, treat_state) %>%
  summarize(
    closures     = mean(closures),
    mergers      = mean(mergers),
    all_changes  = mean(all_changes),
    n_states     = n_distinct(MSTATE),
    .groups = "drop"
  )

closure.plot  <- changes.plot %>%
  ggplot(aes(x=year, y=closures, linetype=as.factor(treat_state))) +
  geom_line(linewidth=0.8, color="black") +
  scale_linetype_manual(values=c("solid", "dashed"), labels=c("Never CAH States", "CAH States")) +
  theme_bw() +
  labs(x="Year", y="Number of closures", linetype="Group") +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.9, .95),
    legend.background = element_blank(),
    legend.key.width = unit(2, "cm") 
  )
ggsave("results/closures.png", closure.plot, width = 6.5, height = 4.25, dpi = 300)


merger.plot  <- changes.plot %>%
  ggplot(aes(x=year, y=mergers, linetype=as.factor(treat_state))) +
  geom_line(linewidth=0.8, color="black") +
  scale_linetype_manual(values=c("solid", "dashed"), labels=c("Never CAH States", "CAH States")) +
  theme_bw() +
  labs(x="Year", y="Number of mergers", linetype="Group") +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.9, .95),
    legend.background = element_blank(),
    legend.key.width = unit(2, "cm") 
  )
ggsave("results/merger.png", merger.plot, width = 6.5, height = 4.25, dpi = 300)


anychange.plot  <- changes.plot %>%
  ggplot(aes(x=year, y=all_changes, linetype=as.factor(treat_state))) +
  geom_line(linewidth=0.8, color="black") +
  scale_linetype_manual(values=c("solid", "dashed"), labels=c("Never CAH States", "CAH States")) +
  theme_bw() +
  labs(x="Year", y="Number of mergers", linetype="Group") +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.9, .95),
    legend.background = element_blank(),
    legend.key.width = unit(2, "cm") 
  )
ggsave("results/all-changes.png", anychange.plot, width = 6.5, height = 4.25, dpi = 300)


## count of critical access by year
view(est.dat %>% group_by(year) %>% summarize(tot_cah=sum(cah)))

## first critical access hospital year by state
view(est.dat %>% filter(cah==1) %>% group_by(MSTATE) %>% summarize(first_cah=min(year)))

## bed sizes of hospitals that closed
bed.closed <- est.dat %>% filter(change_type=="Closure") %>% 
  mutate(bed_bin=cut(BDTOT, breaks=c(0,5,10,15,20,25,30,40,50,75,100,150,200,250,300,400,500,999))) %>% 
  count(bed_bin) %>% arrange(bed_bin)

## CAH designation among hospitals that closed
cah.closed <- est.dat %>% mutate(all_cah=sum(cah)) %>% 
  group_by(change_type, cah) %>% summarize(n=n(), all_cah=first(all_cah), all_hosp=nrow(aha.data))

## Count of CAH and Non-CAH over time
fig.hosp.type <- est.dat %>% group_by(year, cah) %>%
  summarize(hosp_count=n())  %>%
  ggplot(aes(x=year, y=hosp_count, group=cah)) + 
  geom_line() + geom_point() + theme_bw() +
  geom_text(data = est.dat %>% filter(year==2016) %>% group_by(cah, year) %>% summarize(hosp_count=n()), 
            aes(label = c("Non-CAH", "CAH"),
                x = year+1,
                y = hosp_count-100)) +
  scale_y_continuous(labels = comma,
                     breaks=seq(0, 8000, 500)) +
  scale_x_continuous(breaks=seq(1980, 2019, 5)) +
  guides(linetype="none") +
  labs(
    x="Year",
    y="Count of Hospital Types",
    caption="Limited to general, short-term, acute care community hospitals"
  )
ggsave("results/hosp-types.png", fig.hosp.type, width = 6.5, height = 4.25, dpi = 300)

## Share of CAH and Non-CAH closures over time
cah.change <- est.dat %>%
  filter(cah==1) %>%
  group_by(year, change_type) %>%
  summarize(cah_change=n()) %>%
  group_by(year) %>%
  mutate(all_cah=sum(cah_change),
         prop_cah=cah_change/all_cah) %>%
  select(year, prop_cah, change_type) %>%
  ungroup()

non.cah.change <- est.dat %>%
  filter(cah==0) %>%
  group_by(year, change_type) %>%
  summarize(non_cah_change=n()) %>%
  group_by(year) %>%
  mutate(all_non_cah=sum(non_cah_change),
         prop_non_cah=non_cah_change/all_non_cah) %>%
  select(year, prop_non_cah, change_type) %>%
  ungroup()

all.change <- est.dat %>%
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
  ylim(0,0.02) + 
  scale_x_continuous(breaks=seq(1980, 2019, 1)) +
  labs(x = "Year", y = "Proportion", color = "Hospital Type", linetype = "Change Type") +
  theme_minimal() +
  scale_color_discrete(labels=c("CAH","Non-CAH")) + 
  guides(color = guide_legend(order = 2), linetype = guide_legend(order = 1)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        legend.position=c(0.8,0.65))


## Treatment timing --------------------------------------------------------

panel.cah <- panelview(1~treat_post, data=state.dat1 %>% mutate(MSTATE=as.character(MSTATE)), 
                      index=c("MSTATE","year"), legendOff=TRUE, 
          theme.bw=TRUE, by.timing=TRUE, xlab="Year", ylab="State",
          main="", color=c("white","gray"), axis.lab.angle=45)
ggsave("results/panelview-cah-treat.png", panel.cah, width = 6.5, height = 4.25, dpi = 300)

panel.closure <- panelview(closures ~ treat_post,
          data = state.dat1, index = c("MSTATE","year"), type = "outcome",
          main = "", by.cohort = TRUE, outcome.type="discrete",
          legend.labs = c("Never CAH States","Treated States (before CAH)",
                          "Treated States (after CAH)"),
          theme.bw=TRUE, xlab="Year", ylab="Count of Closures")          
ggsave("results/panelview-cah-closure.png", panel.closure, width = 6.5, height = 4.25, dpi = 300)
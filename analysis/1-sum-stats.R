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
  labs(x="Year", y="Closures per State per Year", linetype="Group") +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.85, .85),
    legend.background = element_blank(),
    legend.key.width = unit(2, "cm") 
  )
ggsave("results/desc-closures.png", closure.plot, width = 6.5, height = 4.25, dpi = 300, scale=1.5)


merger.plot  <- changes.plot %>%
  ggplot(aes(x=year, y=mergers, linetype=as.factor(treat_state))) +
  geom_line(linewidth=0.8, color="black") +
  scale_linetype_manual(values=c("solid", "dashed"), labels=c("Never CAH States", "CAH States")) +
  theme_bw() +
  labs(x="Year", y="Mergers per State per Year", linetype="Group") +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.85, .85),
    legend.background = element_blank(),
    legend.key.width = unit(2, "cm") 
  )
ggsave("results/desc-mergers.png", merger.plot, width = 6.5, height = 4.25, dpi = 300, scale=1.5)

anychange.plot  <- changes.plot %>%
  ggplot(aes(x=year, y=all_changes, linetype=as.factor(treat_state))) +
  geom_line(linewidth=0.8, color="black") +
  scale_linetype_manual(values=c("solid", "dashed"), labels=c("Never CAH States", "CAH States")) +
  theme_bw() +
  labs(x="Year", y="Mergers + Closures per State per Year", linetype="Group") +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.85, .85),
    legend.background = element_blank(),
    legend.key.width = unit(2, "cm") 
  )
ggsave("results/desc-all-changes.png", anychange.plot, width = 6.5, height = 4.25, dpi = 300, scale=1.5)


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
                x = year,
                y = hosp_count-110)) +
  scale_y_continuous(labels = comma,
                     breaks=seq(0, 8000, 500)) +
  scale_x_continuous(breaks=seq(1980, 2019, 5)) +
  guides(linetype="none") +
  labs(
    x="Year",
    y="Count of Hospitals"
  )
ggsave("results/desc-hosp-types.png", fig.hosp.type, width = 6.5, height = 4.25, dpi = 300, scale=1.5)


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


fig.change.type <- all.change %>%
  left_join(non.cah.change, by = c("year","change_type")) %>%
  left_join(cah.change,      by = c("year","change_type")) %>%
  filter(change_type %in% c("Closure", "Merger")) %>%
  select(year, change_type, prop_non_cah, prop_cah) %>%
  replace_na(list(prop_non_cah = 0, prop_cah = 0)) %>%
  pivot_longer(
    cols = starts_with("prop"),
    names_to = "type", values_to = "value"
  ) %>%
  mutate(
    type = dplyr::recode(type,
                         prop_non_cah = "Non-CAH",
                         prop_cah     = "CAH"),
    type = factor(type, levels = c("Non-CAH","CAH"))
  ) %>%
  ggplot(aes(x = year, y = value, linetype = change_type, group = change_type)) +
  geom_line(color = "black", size = 0.7) +
  facet_wrap(~ type, nrow = 2) +
  ylim(0, 0.0175) +
  scale_x_continuous(breaks = seq(1980, 2019, 5)) +
  scale_linetype_manual(values = c("Closure" = "solid", "Merger" = "dotted")) +
  labs(x = "Year", y = "Proportion of Hospitals Closing or Merging", linetype = "Change Type") +
  guides(linetype = guide_legend(order = 1)) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )
ggsave("results/desc-changes-by-type.png", fig.change.type, width = 6.5, height = 4.25, dpi = 300, scale=1.5)

## Treatment timing --------------------------------------------------------

panel.cah <- panelview(1~treat_post, data=state.dat %>% mutate(MSTATE=as.character(MSTATE)), 
                      index=c("MSTATE","year"), legendOff=TRUE, 
          theme.bw=TRUE, by.timing=TRUE, xlab="Year", ylab="State",
          main="", color=c("white","gray"), axis.lab.angle=45)
ggsave("results/desc-panelview-treat.png", panel.cah, width = 6.5, height = 4.25, dpi = 300, scale=2)

panel.closure <- panelview(closures ~ treat_post,
          data = state.dat, index = c("MSTATE","year"), type = "outcome",
          main = "", by.cohort = TRUE, outcome.type="discrete",
          legend.labs = c("Never CAH States","Treated States (before CAH)",
                          "Treated States (after CAH)"),
          theme.bw=TRUE, xlab="Year", ylab="Count of Closures")          
ggsave("results/desc-panelview-closures.png", panel.closure, width = 6.5, height = 4.25, dpi = 300, scale=1.5)



## Characteristics of CAH and non-CAH Hospitals ---------------------------------------------------
merge1 <- est.dat %>% filter(year>=1995 & year<=2010) %>%
  filter(ever_cah==1, merged==1) %>% select(ID, year, eff_year, ever_cah, cah)


cah.desc <- est.dat %>% filter(year>=1995 & year<=2010) %>%
  filter(ever_cah==1) %>%
  summarize(
    bed_mean   = mean(BDTOT, na.rm = TRUE),
    bed_p10    = quantile(BDTOT, 0.10, na.rm = TRUE),
    bed_p90    = quantile(BDTOT, 0.90, na.rm = TRUE),
    dist_mean  = mean(distance, na.rm = TRUE),
    dist_p10   = quantile(distance, 0.10, na.rm = TRUE),
    dist_p90   = quantile(distance, 0.90, na.rm = TRUE),
    ip_days    = mean(IPDTOT, na.rm = TRUE),
    ip_days_p10= quantile(IPDTOT, 0.10, na.rm = TRUE),
    ip_days_p90= quantile(IPDTOT, 0.90, na.rm = TRUE),
    rural      = mean(ever_rural, na.rm=TRUE),
    margin_mean=mean(margin, na.rm=TRUE),
    margin_p10 =quantile(margin, 0.10, na.rm=TRUE),
    margin_p90 =quantile(margin, 0.90, na.rm=TRUE),
    cr_mean    = mean(current_ratio, na.rm=TRUE),
    cr_p10     = quantile(current_ratio, 0.10, na.rm=TRUE),
    cr_p90     = quantile(current_ratio, 0.90, na.rm=TRUE),
    system_mean= mean(system, na.rm=TRUE),
    closures   = sum(closed, na.rm=TRUE),
    mergers    = sum(merged, na.rm=TRUE),
    hospitals  = n_distinct(ID),
    obs        = n()
  ) %>%
  mutate(rural = rural * 100,
         system_mean = system_mean * 100) %>%
  mutate(across(c(ip_days, ip_days_p10, ip_days_p90), ~round(., 2))) %>%
  mutate(across(-c(closures, mergers, hospitals, obs), ~round(., 2))) %>%
  mutate(across(c(bed_p10, bed_p90), ~round(.,0))) %>%
  rename_with(~paste0(., "_cah"))


cah.pre <- est.dat %>% filter(year>=1995 & year<=2010) %>%
  filter(ever_cah==1, eff_year >= year) %>%
  summarize(
    bed_mean   = mean(BDTOT, na.rm = TRUE),
    bed_p10    = quantile(BDTOT, 0.10, na.rm = TRUE),
    bed_p90    = quantile(BDTOT, 0.90, na.rm = TRUE),
    dist_mean  = mean(distance, na.rm = TRUE),
    dist_p10   = quantile(distance, 0.10, na.rm = TRUE),
    dist_p90   = quantile(distance, 0.90, na.rm = TRUE),
    ip_days    = mean(IPDTOT, na.rm = TRUE),
    ip_days_p10= quantile(IPDTOT, 0.10, na.rm = TRUE),
    ip_days_p90= quantile(IPDTOT, 0.90, na.rm = TRUE),
    rural      = mean(ever_rural, na.rm=TRUE),
    margin_mean=mean(margin, na.rm=TRUE),
    margin_p10 =quantile(margin, 0.10, na.rm=TRUE),
    margin_p90 =quantile(margin, 0.90, na.rm=TRUE),
    cr_mean    = mean(current_ratio, na.rm=TRUE),
    cr_p10     = quantile(current_ratio, 0.10, na.rm=TRUE),
    cr_p90     = quantile(current_ratio, 0.90, na.rm=TRUE),
    system_mean= mean(system, na.rm=TRUE),
    closures   = sum(closed, na.rm=TRUE),
    mergers    = sum(merged, na.rm=TRUE),
    hospitals  = n_distinct(ID),
    obs        = n()
  ) %>%
  mutate(rural = rural * 100,
         system_mean = system_mean * 100) %>%
  mutate(across(c(ip_days, ip_days_p10, ip_days_p90), ~round(., 2))) %>%
  mutate(across(-c(closures, mergers, hospitals, obs), ~round(., 2))) %>%
  mutate(across(c(bed_p10, bed_p90), ~round(.,0))) %>%
  rename_with(~paste0(., "_cahpre"))


cah.post <- est.dat %>% filter(year>=1995 & year<=2010) %>%
  filter(ever_cah==1, eff_year < year) %>%
  summarize(
    bed_mean   = mean(BDTOT, na.rm = TRUE),
    bed_p10    = quantile(BDTOT, 0.10, na.rm = TRUE),
    bed_p90    = quantile(BDTOT, 0.90, na.rm = TRUE),
    dist_mean  = mean(distance, na.rm = TRUE),
    dist_p10   = quantile(distance, 0.10, na.rm = TRUE),
    dist_p90   = quantile(distance, 0.90, na.rm = TRUE),
    ip_days    = mean(IPDTOT, na.rm = TRUE),
    ip_days_p10= quantile(IPDTOT, 0.10, na.rm = TRUE),
    ip_days_p90= quantile(IPDTOT, 0.90, na.rm = TRUE),
    rural      = mean(ever_rural, na.rm=TRUE),
    margin_mean=mean(margin, na.rm=TRUE),
    margin_p10 =quantile(margin, 0.10, na.rm=TRUE),
    margin_p90 =quantile(margin, 0.90, na.rm=TRUE),
    cr_mean    = mean(current_ratio, na.rm=TRUE),
    cr_p10     = quantile(current_ratio, 0.10, na.rm=TRUE),
    cr_p90     = quantile(current_ratio, 0.90, na.rm=TRUE),
    system_mean= mean(system, na.rm=TRUE),
    closures   = sum(closed, na.rm=TRUE),
    mergers    = sum(merged, na.rm=TRUE),
    hospitals  = n_distinct(ID),
    obs        = n()
  ) %>%
  mutate(rural = rural * 100,
         system_mean = system_mean * 100) %>%
  mutate(across(c(ip_days, ip_days_p10, ip_days_p90), ~round(., 2))) %>%
  mutate(across(-c(closures, mergers, hospitals, obs), ~round(., 2))) %>%
  mutate(across(c(bed_p10, bed_p90), ~round(.,0))) %>%
  rename_with(~paste0(., "_cahpost"))


non.cah.desc <- est.dat %>% filter(year>=1995 & year<=2010) %>%
  filter(ever_cah==0) %>%
  summarize(
    bed_mean   = mean(BDTOT, na.rm = TRUE),
    bed_p10    = quantile(BDTOT, 0.10, na.rm = TRUE),
    bed_p90    = quantile(BDTOT, 0.90, na.rm = TRUE),
    dist_mean  = mean(distance, na.rm = TRUE),
    dist_p10   = quantile(distance, 0.10, na.rm = TRUE),
    dist_p90   = quantile(distance, 0.90, na.rm = TRUE),
    ip_days    = mean(IPDTOT, na.rm = TRUE),
    ip_days_p10= quantile(IPDTOT, 0.10, na.rm = TRUE),
    ip_days_p90= quantile(IPDTOT, 0.90, na.rm = TRUE),
    rural      = mean(ever_rural, na.rm=TRUE),
    margin_mean=mean(margin, na.rm=TRUE),
    margin_p10 =quantile(margin, 0.10, na.rm=TRUE),
    margin_p90 =quantile(margin, 0.90, na.rm=TRUE),
    cr_mean    = mean(current_ratio, na.rm=TRUE),
    cr_p10     = quantile(current_ratio, 0.10, na.rm=TRUE),
    cr_p90     = quantile(current_ratio, 0.90, na.rm=TRUE),
    system_mean= mean(system, na.rm=TRUE),
    closures   = sum(closed, na.rm=TRUE),
    mergers    = sum(merged, na.rm=TRUE),
    hospitals  = n_distinct(ID),
    obs        = n()
  ) %>%
  mutate(rural = rural * 100,
         system_mean = system_mean * 100) %>%
  mutate(across(c(ip_days, ip_days_p10, ip_days_p90), ~round(., 2))) %>%
  mutate(across(-c(closures, mergers, hospitals, obs), ~round(., 2))) %>%
  mutate(across(c(bed_p10, bed_p90), ~round(.,0))) %>%
  rename_with(~paste0(., "_noncah"))

int <- function(lo, hi) sprintf("[%.2f, %.2f]", lo, hi)

make_col <- function(x, suf) {
  c(
    x[[paste0("rural_",       suf)]],
    x[[paste0("bed_mean_",    suf)]],
    int(x[[paste0("bed_p10_", suf)]],    x[[paste0("bed_p90_", suf)]]),
    x[[paste0("dist_mean_",   suf)]],
    int(x[[paste0("dist_p10_",suf)]],    x[[paste0("dist_p90_",suf)]]),
    x[[paste0("ip_days_",     suf)]],
    int(x[[paste0("ip_days_p10_",suf)]], x[[paste0("ip_days_p90_",suf)]]),
    x[[paste0("margin_mean_", suf)]],
    int(x[[paste0("margin_p10_",suf)]],  x[[paste0("margin_p90_",suf)]]),
    x[[paste0("cr_mean_",     suf)]],
    int(x[[paste0("cr_p10_",  suf)]],    x[[paste0("cr_p90_",  suf)]]),
    x[[paste0("system_mean_", suf)]],
    x[[paste0("closures_",    suf)]],
    x[[paste0("mergers_",     suf)]],
    x[[paste0("hospitals_",   suf)]],
    x[[paste0("obs_",         suf)]]
  )
}

desc_tab <- tibble(
  Variable = c(
    "Rural (\\%)",
    "Bed size", "",
    "Distance to nearest hospital", "",
    "Inpatient days per bed", "",
    "Operating margin", "",
    "Current ratio", "",
    "System member (\\%)",
    "Closures",
    "Mergers",
    "Hospitals",
    "Observations"
  ),
  `Pre`       = make_col(cah.pre,      "cahpre"),
  `Post`      = make_col(cah.post,     "cahpost"),
  `Ever CAH`  = make_col(cah.desc,     "cah"),
  `Never CAH` = make_col(non.cah.desc, "noncah")
)

# LaTeX table innards only (for \input in paper.tex)
tex_lines <- c(
  "\\toprule",
  "\\multicolumn{1}{c}{ } & \\multicolumn{3}{c}{CAHs} & \\multicolumn{1}{c}{ } \\\\",
  "\\cmidrule(l{3pt}r{3pt}){2-4}",
  " & Pre & Post & Ever CAH & Never CAH\\\\",
  "\\midrule"
)

for (i in seq_len(nrow(desc_tab))) {
  row_vals <- c(desc_tab$Variable[i], desc_tab$Pre[i], desc_tab$Post[i],
                desc_tab$`Ever CAH`[i], desc_tab$`Never CAH`[i])
  tex_lines <- c(tex_lines, paste(row_vals, collapse = " & ") %>% paste0(" \\\\"))
}

tex_lines <- c(tex_lines, "\\bottomrule")

writeLines(tex_lines, "results/desc_compare.tex")


# Summary Statistics -------------------------------------------------------

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


## Treatment timing --------------------------------------------------------

panel.cah <- panelview(1~treat_post, data=state.dat %>% mutate(MSTATE=as.character(MSTATE)),
                      index=c("MSTATE","year"), legendOff=TRUE,
          theme.bw=TRUE, by.timing=TRUE, xlab="Year", ylab="State",
          main="", color=c("white","gray"), axis.lab.angle=45)
ggsave("results/desc-panelview-treat.png", panel.cah, width = 6.5, height = 4.25, dpi = 300, scale=2)


## Characteristics of CAH and non-CAH Hospitals ---------------------------------------------------

cah.desc <- est.dat %>% filter(year>=1995 & year<=2010) %>%
  filter(ever_cah==1) %>%
  summarize(
    bed_mean   = mean(BDTOT, na.rm = TRUE),
    bed_p10    = quantile(BDTOT, 0.10, na.rm = TRUE),
    bed_p90    = quantile(BDTOT, 0.90, na.rm = TRUE),
    dist_mean  = mean(distance, na.rm = TRUE),
    dist_p10   = quantile(distance, 0.10, na.rm = TRUE),
    dist_p90   = quantile(distance, 0.90, na.rm = TRUE),
    ip_days    = mean(ip_per_bed, na.rm = TRUE),
    ip_days_p10= quantile(ip_per_bed, 0.10, na.rm = TRUE),
    ip_days_p90= quantile(ip_per_bed, 0.90, na.rm = TRUE),
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
    ip_days    = mean(ip_per_bed, na.rm = TRUE),
    ip_days_p10= quantile(ip_per_bed, 0.10, na.rm = TRUE),
    ip_days_p90= quantile(ip_per_bed, 0.90, na.rm = TRUE),
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
    ip_days    = mean(ip_per_bed, na.rm = TRUE),
    ip_days_p10= quantile(ip_per_bed, 0.10, na.rm = TRUE),
    ip_days_p90= quantile(ip_per_bed, 0.90, na.rm = TRUE),
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
    ip_days    = mean(ip_per_bed, na.rm = TRUE),
    ip_days_p10= quantile(ip_per_bed, 0.10, na.rm = TRUE),
    ip_days_p90= quantile(ip_per_bed, 0.90, na.rm = TRUE),
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

int <- function(lo, hi) sprintf("{}[%.2f, %.2f]", lo, hi)

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

# Complete LaTeX tabular block (LaTeX 2025 breaks \noalign/\omit inside \input within tabular)
tex_lines <- c(
  "\\begin{tabular}[t]{lcccc}",
  "\\multicolumn{1}{c}{ } & \\multicolumn{3}{c}{CAHs} & \\multicolumn{1}{c}{ } \\\\",
  "\\cmidrule(l{3pt}r{3pt}){2-4}",
  " & Pre & Post & Ever CAH & Never CAH\\\\"
)

for (i in seq_len(nrow(desc_tab))) {
  row_vals <- c(desc_tab$Variable[i], desc_tab$Pre[i], desc_tab$Post[i],
                desc_tab$`Ever CAH`[i], desc_tab$`Never CAH`[i])
  tex_lines <- c(tex_lines, paste(row_vals, collapse = " & ") %>% paste0(" \\\\"))
}

tex_lines <- c(tex_lines, "\\end{tabular}")

writeLines(tex_lines, "results/desc_compare.tex")

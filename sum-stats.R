

# Quick summary -----------------------------------------------------------

aha.combine <- read_csv('data/output/aha_data.csv')
bed.closed <- aha.combine %>% filter(change_type=="Closure") %>% 
  mutate(bed_bin=cut(BDTOT, breaks=c(0,5,10,15,20,25,30,40,50,75,100,150,200,250,300,400,500,1000))) %>% 
  count(bed_bin) %>% arrange(bed_bin)

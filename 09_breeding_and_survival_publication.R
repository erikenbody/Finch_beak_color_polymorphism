library(tidyverse)
library(lme4)
library(lmerTest)
library(kableExtra)
library(car)
library(sjPlot)
library(binom)


df.scandens.bco2.breeding <- read.csv("output/breeding/scandens_BCO2_breeding_V2.csv")
# -------------------------------------------------------------------------

df.scandens.bco2.fy <- read.csv("output/breeding/scandens_1998_1991_BCO2_breeding_V2.csv")

df.scandens.bco2.fy$First.year.min <- as.factor(df.scandens.bco2.fy$First.year.min)

df.scandens.bco2.fy$BCO2 = relevel(as.factor(df.scandens.bco2.fy$BCO2), ref = "G/A")


m.scandens.bco2.first.year <- glm(first.year.surv ~ BCO2 * First.year.min, 
                                  data = df.scandens.bco2.fy, family = binomial(link = "logit"))

ci.method <- "agresti-coull"
df.scandens.bco2.fy.cohort.S <- df.scandens.bco2.fy %>% 
  mutate(SURVIVAL = ifelse(first.year.surv == 0, "DIED", 
                           ifelse(first.year.surv == 1, "SURVIVED", NA))) %>% 
  group_by(BCO2, SURVIVAL, First.year.min) %>% 
  summarise(n = n()) %>% pivot_wider(names_from = SURVIVAL, values_from = n) %>% 
  mutate(n = DIED + SURVIVED)

df.x <- binom.confint(x=df.scandens.bco2.fy.cohort.S$SURVIVED, n=df.scandens.bco2.fy.cohort.S$n, methods=ci.method) %>% as.data.frame() %>% select(-n)
df.scandens.bco2.fy.cohort.S <- cbind(as.data.frame(df.scandens.bco2.fy.cohort.S), df.x)
df.scandens.bco2.fy.cohort.S$BCO2  <- factor(df.scandens.bco2.fy.cohort.S$BCO2, levels = c("A/A","G/A","G/G"))

ggplot(df.scandens.bco2.fy.cohort.S, aes(x = First.year.min, y = mean, color = BCO2, group = BCO2)) + 
  geom_point(aes(color = BCO2), size = 4,  position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = lower, ymax = upper, color = BCO2), size = 1, position = position_dodge(width = 0.5)) +
  geom_text(aes(y = 1.1, label = n), color = "black", position = position_dodge(width = 0.5)) +
  scale_color_manual(values = colors2) + theme_classic() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14)) + 
  ylab("First-year survival") + xlab(NULL) + 
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1.1, by = .25)) +
  labs(color = NULL)

ggsave("output/breeding/simplified_output/first_year_survival_scandens.pdf", 
       width = 6, height = 5, units = "in", useDingbats=FALSE)

tab_model(m.scandens.bco2.first.year)

# -------------------------------------------------------------------------

df.scandens.bco2.hatch.rate <- read.csv("output/breeding/scandens_hatch_rate_breeding_BCO2_V2.csv")
df.scandens.bco2.hatch.rate.n <- df.scandens.bco2.hatch.rate %>% 
  filter(!is.na(NESTOTAL.mean.hatch.rate)) %>% 
  group_by(BCO2) %>% 
  summarise(mean = mean(NESTOTAL.mean.hatch.rate, na.rm = T),
            std = sqrt(var(NESTOTAL.mean.hatch.rate, na.rm = T)),
            lower = mean(NESTOTAL.mean.hatch.rate, na.rm = T) - qnorm(.975)*std/sqrt(n()),
            upper = mean(NESTOTAL.mean.hatch.rate, na.rm = T) + qnorm(.975)*std/sqrt(n()),
            se = std/sqrt(n()),
            n = n(), n.nests = sum(NESTOTAL.nests)) %>% 
  mutate(Sex = "Female", label = paste0(n," (",n.nests,")"))

#pretty plot
df.scandens.bco2.hatch.rate.n.BOTH <- df.scandens.bco2.hatch.rate.n
df.scandens.bco2.hatch.rate.n.BOTH$BCO2 = factor(df.scandens.bco2.hatch.rate.n.BOTH$BCO2, levels = c("A/A","G/A","G/G")) 
df.scandens.bco2.hatch.rate.n.BOTH <- df.scandens.bco2.hatch.rate.n.BOTH %>% mutate(code = ifelse(Sex == "Male", "Male", as.character(BCO2)),
                                                                                    label = ifelse(Sex == "Male", NA, as.character(label)),
                                                                                    upper = ifelse(upper > 1, 1, upper)) #cutt off error bars at 1 for visualization purposes

p.df.scandens.bco2.hatch.rate.n.BOTH <- ggplot(df.scandens.bco2.hatch.rate.n.BOTH, aes(BCO2, mean, color = code)) +
  geom_point(size = 4, position = position_dodge(width = 0.7)) +
  geom_linerange(aes(ymin = lower, ymax = upper), size = 1, position = position_dodge(width = 0.7)) +
  geom_text(aes(y = 1.1, label = label, group = code), color = "black", position = position_dodge(width = 0.7)) +
  scale_color_manual(values = c( 'orange2', "grey40",'#e9a3c9', "grey90")) + theme_classic() +
  ylab("Hatching success") + xlab(NULL) + 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14)) +
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1.1, by = .25)) +
  labs(color = NULL)


p.df.scandens.bco2.hatch.rate.n.BOTH

ggsave("output/breeding/simplified_output/publication_hatch_plot.pdf", 
       width = 6, height = 5, units = "in", dpi = "retina", useDingbats=FALSE)

df.scandens.bco2.hatch.rate$BCO2 = relevel(as.factor(df.scandens.bco2.hatch.rate$BCO2), ref = "G/A") 

m.scandens.bco2.hatch.rate<-lmerTest::lmer(NESTOTAL.mean.hatch.rate ~ BCO2 + (1 | First.year.min), 
                                           data = df.scandens.bco2.hatch.rate)

tab_model(m.scandens.bco2.hatch.rate)


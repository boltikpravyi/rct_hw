# Load packages
# library(readxl)
# library(tidyverse)
# # library(caret)
# # library(psych)
# # library(stargazer)
# library(censReg)
# # library(AER)
# # library(rpart)
# # library(rpart.plot)
# library(MatchIt) # propensity score matching
# library(sandwich) # heterokedasticity robust standard errors
# library(lmtest) # coefficient testing

# Load data
rawdata <- read_excel(path = 'data.xlsx')

# Transform data
df <- rawdata %>% 
  rename(ID = 'ID пациента',
         Tr = 'T') %>% 
  mutate(Yr = round(Y, 2)) %>% # for censreg estimates
  # mutate(Xr = round(X, 2)) %>% # for PS estimates
  mutate(Xc = X - mean(X)) %>% # NB! X*Tr interaction term always improves ATE precion estimate (from Causal ML book)
  # arrange(Y) %>% 
  # mutate(Treatment = as_factor(Treatment)) #%>% 
  mutate(Label = as_factor(Tr)) #%>%
  # mutate(Label = case_when(Treatment == 0 ~ 'P', Treatment == 1 ~ 'T')) %>% 
  # mutate(Label = as_factor(Label))

# 1. Exploratory data analysis -----------------------------------------------------------

# Observed outcome across the patients
df %>% 
  ggplot(mapping = aes(x = ID, y = Y, color = Label)) +
  geom_point() +
  theme_bw() +
  theme(legend.position = 'top') +
  labs(x = 'Patient', y = 'Observed health outcome') +
  scale_color_discrete(name = '', labels = c('Placebo', 'Treatment'))

# Distribution of observed outcome
df %>% 
  ggplot(mapping = aes(x = Y, fill = Label)) +
    geom_histogram(bins = 250) +
    theme_bw() +
    theme(legend.position = 'top') + 
    labs(x = 'Observed health outcome', y = 'No of evidence') +
    scale_fill_discrete(name = '', labels = c('Placebo', 'Treatment'))

# Descriptive statistics
descriptive_stats <- df %>% 
  select(-ID) %>% 
  group_by(Tr) %>% 
  summarise(n = n(),
            distinct_y = n_distinct(Y),
            mean_y = mean(Y), sd_y = sd(Y),
            min_y = min(Y), max_y = max(Y),
            median_y = median(Y),
            distinct_x = n_distinct(X),
            mean_x = mean(X), sd_x = sd(X),
            min_x = min(X), max_x = max(X),
            median_x = median(X))
descriptive_stats %>% 
  select(Tr, n, mean_y, sd_y, min_y, max_y, median_y)
descriptive_stats %>% 
  select(Tr, n, mean_x, sd_x, min_x, max_x, median_x)

# Data on pre-treatment patients' characteristics
df %>% 
  ggplot(mapping = aes(x = X, fill = Label)) +
  geom_histogram(aes(y = ..density..), position = 'identity', alpha = 0.33, bins = 50) +
  theme_bw() +
  theme(legend.position = 'top') + 
  labs(x = 'Pre-treatment K steps per day', y = 'Density') +
  scale_fill_discrete(name = '', labels = c('Placebo', 'Treatment'))

# No statistical difference between treated and control group characteristics
lm(data = df, formula = X ~ Tr) %>% summary()

# How observed outcome is related to pre-treatment patients' characteristics
df %>% 
  ggplot(mapping = aes(x = X, y = Y, color = Label)) +
  geom_point() +
  stat_smooth(method = 'lm', formula = y ~ x, geom = 'smooth') +
  theme_bw() +
  theme(legend.position = 'top') + 
  labs(x = 'Pre-treatment K steps per day', y = 'Observed health outcome') +
  scale_color_discrete(name = '', labels = c('Placebo', 'Treatment'))

# Association between observed outcome and pre-treatment patients' characteristics
lm(Y ~ X, data = df) %>% summary()

# 2. Estimation of ATE -------------------------------------------------------------------

# OLS-based estimator based on linear regression (= Diff-in-means estimator)
estOLS <- lm(Y ~ Tr, data = df)
# estOLS <- lm(Y ~ Tr + Xc*Tr, data = df)
summary(estOLS)
# confint(estOLS) # NB! Use standard errors robusted to misspec, hc and ac (from Causal ML book)
# Heteroskedasticity of error term
resOLS <- estOLS$residuals
plot(y = resOLS^2, x = df$Tr)
plot(y = resOLS^2, x = df$Xc)
lm(resOLS^2 ~ df$Tr + df$Xc) %>% summary() # significant presence of hc
coeftest(estOLS, vcov = vcovHC(estOLS, type = 'HC1')) # hc robust inference (Eicker-Huber-White)
summary(estOLS) # classical (non-robust) inference

# Censored regression estimator (removes downward bias of OLS estimator)
estCens <- censReg(Yr ~ Tr, data = df, left = min(df$Yr), right = max(df$Yr))
# estCens <- censReg(Yr ~ Tr + Xc*Tr, data = df, left = min(df$Yr), right = max(df$Yr))
summary(estCens)
margEff(estCens) # Marginal effect is different from coeficient estimate
confint(estCens) # NB! Use standard errors robusted to misspec, hc and ac (from Causal ML book)
# NB! Add delta method for correct s.e.
# coeftest(estCens, vcov = vcovHC(estCens, type = 'HC1')) # hc robust inference (Eicker-Huber-White)

# Matching treatment and control
estMatch <- matchit(Tr ~ Xc, data = df, method = 'nearest') # 1:1 greedy NN on the PS
summary(estMatch)
dfMatch <- (df %>% 
  filter(ID %in% as.numeric(estMatch$match.matrix[,1]))) %>% 
  bind_rows(
    df %>% 
      filter(Tr == 1)
  )

# EDA for matched observations ----

# Observed outcome across the patients
dfMatch %>% 
  ggplot(mapping = aes(x = ID, y = Y, color = Label)) +
  geom_point() +
  theme_bw() +
  theme(legend.position = 'top') +
  labs(x = 'Patient', y = 'Observed health outcome') +
  scale_color_discrete(name = '', labels = c('Placebo', 'Treatment'))

# Distribution of observed outcome
dfMatch %>% 
  ggplot(mapping = aes(x = Y, fill = Label)) +
  geom_histogram(bins = 250) +
  theme_bw() +
  theme(legend.position = 'top') + 
  labs(x = 'Observed health outcome', y = 'No of evidence') +
  scale_fill_discrete(name = '', labels = c('Placebo', 'Treatment'))

# Descriptive statistics
descriptive_stats <- dfMatch %>% 
  select(-ID) %>% 
  group_by(Tr) %>% 
  summarise(n = n(),
            distinct_y = n_distinct(Y),
            mean_y = mean(Y), sd_y = sd(Y),
            min_y = min(Y), max_y = max(Y),
            median_y = median(Y),
            distinct_x = n_distinct(X),
            mean_x = mean(X), sd_x = sd(X),
            min_x = min(X), max_x = max(X),
            median_x = median(X))
descriptive_stats %>% 
  select(Tr, n, mean_y, sd_y, min_y, max_y, median_y)
descriptive_stats %>% 
  select(Tr, n, mean_x, sd_x, min_x, max_x, median_x)

# Data on pre-treatment patients' characteristics
dfMatch %>% 
  ggplot(mapping = aes(x = X, fill = Label)) +
  geom_histogram(aes(y = ..density..), position = 'identity', alpha = 0.33, bins = 50) +
  theme_bw() +
  theme(legend.position = 'top') + 
  labs(x = 'Pre-treatment K steps per day', y = 'Density') +
  scale_fill_discrete(name = '', labels = c('Placebo', 'Treatment'))

# No statistical difference between treated and control group characteristics
lm(X ~ Tr, data = dfMatch) %>% summary()

# How observed outcome is related to pre-treatment patients' characteristics
dfMatch %>% 
  ggplot(mapping = aes(x = X, y = Y, color = Label)) +
  geom_point() +
  stat_smooth(method = 'lm', formula = y ~ x, geom = 'smooth') +
  theme_bw() +
  theme(legend.position = 'top') + 
  labs(x = 'Pre-treatment K steps per day', y = 'Observed health outcome') +
  scale_color_discrete(name = '', labels = c('Placebo', 'Treatment'))

# Association between observed outcome and pre-treatment patients' characteristics
lm(Y ~ X, data = dfMatch) %>% summary()

# OLS estimator w/ matched observations
estOLSMatch <- lm(Y ~ Tr, data = dfMatch)
# estOLSMatch <- lm(Y ~ Tr + Xc*Tr, data = dfMatch)
summary(estOLSMatch)
confint(estOLSMatch) # NB! Use standard errors robusted to misspec, hc and ac (from Causal ML book)
# Heteroskedasticity of error term
resOLSMatch <- estOLSMatch$residuals
plot(y = resOLSMatch^2, x = dfMatch$Tr)
plot(y = resOLSMatch^2, x = dfMatch$Xc)
lm(resOLSMatch^2 ~ dfMatch$Tr + dfMatch$Xc) %>% summary() # significant presence of hc

# Censored regression estimator w/ matched observations
estCensMatch <- censReg(Yr ~ Tr, data = dfMatch, left = min(dfMatch$Yr), right = max(dfMatch$Yr))
# estCensMatch <- censReg(Yr ~ Tr + Xc*Tr, data = dfMatch, left = min(dfMatch$Yr), right = max(dfMatch$Yr))
summary(estCensMatch)
margEff(estCensMatch) # Marginal effect is different from coeficient estimate
confint(estCensMatch) # NB! Use standard errors robusted to misspec, hc and ac (from Causal ML book)
# NB! Add delta method for correct s.e.

# Estimates of ATE and 95% confidence interval
ATEest <- tibble(OLS = coef(estOLS)[2], CensReg = margEff(estCens)[1],
                 OLSMatch = coef(estOLSMatch)[2], CensRegMatch = margEff(estCensMatch)[1])
ATEconfint <- tibble(OLS = confint(estOLS)[2,], CensReg = confint(estCens)[2,],
                     OLSMatch = confint(estOLSMatch)[2,], CensRegMatch = confint(estCensMatch)[2,])
ATEprec <- ATEconfint %>% 
  mutate(OLS = diff(OLS), CensReg = diff(CensReg),
         OLSMatch = diff(OLSMatch), CensRegMatch = diff(CensRegMatch)) %>% .[1,]

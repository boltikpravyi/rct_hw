# 0. Preliminary steps -------------------------------------------------------------------

# Load data
rawdata <- read_excel(path = 'data.xlsx')

# Transform data
df <- rawdata %>% 
  rename(ID = 'ID пациента',
         D = 'T') %>% 
  mutate(Yr = round(Y, 2)) %>% # for censored regression estimates
  mutate(Label = as_factor(D)) %>%
  mutate(X2 = X^2, X3 = X^3, X4 = X^4, X5 = X^5, X6 = X^6) %>% 
  mutate(Xc = X - mean(X), X2c = X2 - mean(X2), X3c = X3 - mean(X3), 
         X4c = X4 - mean(X4), X5c = X5 - mean(X5), X6c = X6 - mean(X6))

# Matching treatment and control
estMatch <- matchit(D ~ Xc, data = df, method = 'nearest') # 1:1 greedy NN on the PS
dfMatch <- (df %>% dplyr::filter(ID %in% as.numeric(estMatch$match.matrix[,1]))) %>%
  bind_rows(df %>% dplyr::filter(D == 1)) %>% 
  mutate(X2 = X^2, X3 = X^3, X4 = X^4, X5 = X^5, X6 = X^6) %>% 
  mutate(Xc = X - mean(X), X2c = X2 - mean(X2), X3c = X3 - mean(X3), 
         X4c = X4 - mean(X4), X5c = X5 - mean(X5), X6c = X6 - mean(X6))

# 1. Exploratory data analysis -----------------------------------------------------------

# Number of treated (1) and untreated (0)
df %>% 
  group_by(D) %>% 
  summarise(Y = n())

# Distribution of observed outcomes
df %>% 
  ggplot(mapping = aes(x = Y, fill = Label)) +
    geom_histogram(bins = 250) +
    theme_bw() +
    theme(legend.position = 'top', plot.title = element_text(hjust = 0.5)) + 
    labs(x = 'Observed health outcome (Y)', y = 'No of evidence', title = 'No matching') +
    scale_fill_discrete(name = '', labels = c('Control (D=0)', 'Treatment (D=1)'))
dfMatch %>% 
  ggplot(mapping = aes(x = Y, fill = Label)) +
  geom_histogram(bins = 250) +
  theme_bw() +
  theme(legend.position = 'top', plot.title = element_text(hjust = 0.5)) + 
  labs(x = 'Observed health outcome (Y)', y = 'No of evidence', title = 'Exact matching') +
  scale_fill_discrete(name = '', labels = c('Control (D=0)', 'Treatment (D=1)'))

# Descriptive statistics
desc <- df %>% 
  select(-ID) %>% 
  group_by(D) %>% 
  summarise(n = n(),
            distinct_y = n_distinct(Y),
            mean_y = mean(Y), sd_y = sd(Y),
            min_y = min(Y), max_y = max(Y),
            median_y = median(Y),
            distinct_x = n_distinct(X),
            mean_x = mean(X), sd_x = sd(X),
            min_x = min(X), max_x = max(X),
            median_x = median(X))
desc %>% 
  select(D, n, mean_y, sd_y, min_y, max_y, median_y)
desc %>% 
  select(D, n, mean_x, sd_x, min_x, max_x, median_x)
descMatch <- dfMatch %>% 
  select(-ID) %>% 
  group_by(D) %>% 
  summarise(n = n(),
            distinct_y = n_distinct(Y),
            mean_y = mean(Y), sd_y = sd(Y),
            min_y = min(Y), max_y = max(Y),
            median_y = median(Y),
            distinct_x = n_distinct(X),
            mean_x = mean(X), sd_x = sd(X),
            min_x = min(X), max_x = max(X),
            median_x = median(X))
descMatch %>% 
  select(D, n, mean_y, sd_y, min_y, max_y, median_y)
descMatch %>% 
  select(D, n, mean_x, sd_x, min_x, max_x, median_x)

# Pre-treatment patients' characteristics
df %>% 
  ggplot(mapping = aes(x = X, fill = Label)) +
  geom_histogram(aes(y = after_stat(density)), position = 'identity', alpha = 0.33, bins = 50) +
  theme_bw() +
  theme(legend.position = 'top', plot.title = element_text(hjust = 0.5)) + 
  labs(x = 'Pre-treatment thousands of steps per day (X)', y = 'Density', title = 'No matching') +
  scale_fill_discrete(name = '', labels = c('Control (D=0)', 'Treatment (D=1)'))
dfMatch %>% 
  ggplot(mapping = aes(x = X, fill = Label)) +
  geom_histogram(aes(y = after_stat(density)), position = 'identity', alpha = 0.33, bins = 50) +
  theme_bw() +
  theme(legend.position = 'top', plot.title = element_text(hjust = 0.5)) + 
  labs(x = 'Pre-treatment thousands of steps per day (X)', y = 'Density', title = 'Exact matching') +
  scale_fill_discrete(name = '', labels = c('Control (D=0)', 'Treatment (D=1)'))

# No statistical difference between treated and control group characteristics
lm(D ~ X, data = df) %>% 
  coeftest(vcov = vcovHC(., type = 'HC3')) # hc-consistent s.e. for robust statistical inference
lm(D ~ X, data = dfMatch) %>%
  coeftest(vcov = vcovHC(., type = 'HC3'))
models <- list('Panel A: No Matching' = list('OLS' = lm(D ~ X, data = df)),
               'Panel B: Exact Matching' = list('OLS' = lm(D ~ X, data = dfMatch)))
modelsummary(models, shape = 'cbind', vcov = 'HC3', coef_rename = c('X' = 'No. of steps (X)'),
             statistic = c('conf.int', 's.e. = {std.error}', 't = {statistic}', 'p = {p.value}'),
             conf_level = .99, gof_omit = c('IC|RMSE|Log|R2$'))

# How observed outcome is related to pre-treatment patients' characteristics
df %>% 
  ggplot(mapping = aes(x = X, y = Y, color = Label)) +
  geom_point() +
  stat_smooth(method = 'lm', formula = y ~ x, geom = 'smooth') +
  theme_bw() +
  theme(legend.position = 'top', plot.title = element_text(hjust = 0.5)) + 
  labs(x = 'Pre-treatment thousands of steps per day (X)', y = 'Observed health outcome (Y)', title = 'No matching') +
  scale_color_discrete(name = '', labels = c('Control (D=0)', 'Treatment (D=1)'))
dfMatch %>% 
  ggplot(mapping = aes(x = X, y = Y, color = Label)) +
  geom_point() +
  stat_smooth(method = 'lm', formula = y ~ x, geom = 'smooth') +
  theme_bw() +
  theme(legend.position = 'top', plot.title = element_text(hjust = 0.5)) + 
  labs(x = 'Pre-treatment thousands of steps per day (X)', y = 'Observed health outcome (Y)', title = 'Exact matching') +
  scale_color_discrete(name = '', labels = c('Control (D=0)', 'Treatment (D=1)'))

# Make causal discovery procedure to understand connections between variables
suff_stat <- list(C = cor(df[,2:4]), n = nrow(df[,2:4]))
estPC <- pc(suff_stat, indepTest = gaussCItest, labels = colnames(df[,2:4]), alpha = 0.01, skel.method = 'stable.fast')
plot(estPC, main = 'Directed acyclic graph (based on PC algorithm)')

# 2. Estimation of ATE -------------------------------------------------------------------

# OLS estimators (no matching)
estOLS <- lm(Y ~ D, data = df) # Diff-in-means estimator
summary(estOLS)
coeftest(estOLS, vcov = vcovHC(estOLS, type = 'HC3')) # hc robust inference (Eicker-Huber-White)
mefOLS <- avg_slopes(estOLS, variables = 'D', vcov = 'HC3')
mefOLS
estOLSControl <- lm(Y ~ D + Xc*D + X2c*D + X3c*D, data = df) # Interaction term always improves the precision of ATE
coeftest(estOLSControl, vcov = vcovHC(estOLSControl, type = 'HC3')) # hc robust inference (Eicker-Huber-White)
mefOLSControl <- avg_slopes(estOLSControl, variables = 'D', vcov = 'HC3')
mefOLSControl

# OLS estimators (exact matching)
estOLSMatch <- lm(Y ~ D, data = dfMatch)
coeftest(estOLSMatch, vcov = vcovHC(estOLSMatch, type = 'HC3'))
mefOLSMatch <- avg_slopes(estOLSMatch, variables = 'D', vcov = 'HC3')
mefOLSMatch
estOLSMatchControl <- lm(Y ~ D + Xc*D + X2c*D + X3c*D, data = dfMatch)
coeftest(estOLSMatchControl, vcov = vcovHC(estOLSMatchControl, type = 'HC3'))
mefOLSMatchControl <- avg_slopes(estOLSMatchControl, variables = 'D', vcov = 'HC3')
mefOLSMatchControl

# Censored ML estimators (no matching)
estCens <- tobit(Yr ~ D, data = df, left = min(df$Yr), right = max(df$Yr), robust = T)
summary(estCens)
mefCens <- avg_slopes(estCens, variables = 'D')
mefCens
estCensControl <- tobit(Yr ~ D + Xc*D + X2c*D + X3c*D, data = df, left = min(df$Yr), right = max(df$Yr), robust = T)
summary(estCensControl)
mefCensControl <- avg_slopes(estCensControl, variables = 'D')
mefCensControl

# Censored ML estimators (exact matching)
estCensMatch <- tobit(Yr ~ D, data = dfMatch, left = min(dfMatch$Yr), right = max(dfMatch$Yr), robust = T)
summary(estCensMatch)
mefCensMatch <- avg_slopes(estCensMatch, variables = 'D')
mefCensMatch
estCensMatchControl <- tobit(Yr ~ D + Xc*D + X2c*D + X3c*D, data = dfMatch, left = min(dfMatch$Yr), right = max(dfMatch$Yr), robust = T)
summary(estCensMatchControl)
mefCensMatchControl <- avg_slopes(estCensMatchControl, variables = 'D')
mefCensMatchControl

# Censored ML Estimators w/ heteroscedastic errors
estCensHC <- crch(Yr ~ D | D + Xc + X2c + X3c + X4c, data = df, left = min(df$Yr), right = max(df$Yr))
summary(estCensHC)
mefCensHC <- avg_slopes(estCensHC, variables = 'D')
mefCensHC
estCensHCMatch <- crch(Yr ~ D | D + Xc + X2c + X3c + X4c, data = dfMatch, left = min(dfMatch$Yr), right = max(dfMatch$Yr))
summary(estCensHCMatch)
mefCensHCMatch <- avg_slopes(estCensHCMatch, variables = 'D')
mefCensHCMatch

# 3. Results -----------------------------------------------------------------------------

mefs <- list(
  'Panel A: No matching' = list(
    'OLS' = mefOLS,
    'OLS w/ X' = mefOLSControl,
    'CR' = mefCens,
    'CR w/ X' = mefCensControl,
    'CR w/ het. err.' = mefCensHC),
  'Panel B: Exact matching' = list(
    'OLS' = mefOLSMatch,
    'OLS w/ X' = mefOLSMatchControl,
    'CR' = mefCensMatch,
    'CR w/ X' = mefCensMatchControl,
    'CR w/ het. err.' = mefCensHCMatch))
modelsummary(mefs, shape = 'rbind', stars = T, statistic = c('s.e. = {std.error}', 'conf.int'),
             coef_rename = c('D' = 'PATE, point estimate'), conf_level = .99,
             gof_omit = c('IC|Log|R2|F|Std|r2'))

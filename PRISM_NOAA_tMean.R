library(dplyr)
library(data.table)
library(lfe)
library(broom)
library(magrittr)
library(expm)
library(zoo)

# Load PRISM overall data
 suicide_month <- read.csv('~/Desktop/county_by_month/PRISM_Lung_cancer_1960_2019.csv')

# Load PRISM Firearm data
# suicide_month <- read.csv('~/Desktop/county_by_month/county_to_month_PRISM_firearm_1970_2019.csv')

# Load PRISM Nonfirearm data
# suicide_month <- read.csv('~/Desktop/county_by_month/county_to_month_PRISM_nonfirearm_1970_2019.csv')

# Load NOAA overall data
# suicide_month <- read.csv('~/Desktop/county_by_month/county_to_month_NOAA_overall_1970_2020.csv')

# Load NOAA Firearm data
# suicide_month <- read.csv('~/Desktop/county_by_month/county_to_month_NOAA_firearm_1970_2020.csv')

# Load NOAA Nonfirearm data
# suicide_month <- read.csv('~/Desktop/county_by_month/county_to_month_NOAA_nonfirearm_1970_2020.csv')


#subset_month <- suicide_month[suicide_month$year >= 1970 & data$year <= 1981, ]
# suicide_month <- suicide_month %>%
# filter(year >= 2000 & year <= 2019)


 


# Create fixed effects variables
suicide_month <- suicide_month %>%
  mutate(cm_group = interaction(fips, month),
         sm_group = interaction(statefips, month),
         sy_group = interaction(statefips, year))

# Create ym variable
suicide_month$ym <- as.yearmon(paste(suicide_month$year, suicide_month$month, sep="-"), "%Y-%m")

# Average county population, weight cannot vary within panel for felm
suicide_month <- suicide_month %>%
  group_by(fips) %>%
  mutate(meanpop = mean(pop))


# Choose your temperature bins
temp_bins <- c("tMean_under30", "tMean_3040", "tMean_4050", "tMean_5060", "tMean_7080", "tMean_over80")

################# THI #################################
# Create lags in temperature bins
for (i in 1:6) {
  for (var in temp_bins) {
    suicide_month <- suicide_month %>%
      arrange(fips, year, month) %>%
      group_by(fips) %>%
      mutate(!!paste(var, "_lag", i, sep="") := lag(!!sym(var), i))
  }
}

# Calculate sum of suicides and store mean
suicide_month$sumsuic <- sum(suicide_month$deaths)
num <- mean(suicide_month$sumsuic)

# Set up matrix to store results
A <- matrix(NA, nrow = 1, ncol = 6)
B <- matrix(NA, nrow = 1, ncol = 6)
C <- matrix(NA, nrow = 1, ncol = 6)
D <- matrix(NA, nrow = 1, ncol = 6)

# Calculate mean of dependent variable
varmean <- weighted.mean(suicide_month$lung_cancer_mortality_rate, w = suicide_month$meanpop)

##############################Model1########################################
# Estimate model with felm (fixed effects lm)
form <- as.formula(paste('lung_cancer_mortality_rate ~', paste(temp_bins, collapse = '+'), '+ lowp + highp | fips + year + month | 0 | statefips'))
mod <- felm(form, data = suicide_month, weights = suicide_month$meanpop)

# Store the estimation results
spec_test_full_1 <- tidy(mod)

# Fill matrices A and B with the estimates and standard errors
i <- 0
for (var in temp_bins) {
  estimate <- spec_test_full_1$estimate[spec_test_full_1$term == var]
  se <- spec_test_full_1$std.error[spec_test_full_1$term == var]
  A[1, 1 + i] <- estimate
  B[1, 1 + i] <- se
  i <- i + 1
}

# Calculate C and D matrices
C <- A / varmean
D <- B / varmean

# Name the columns
colnames(A) <- c("<30", "30-40", "40-50", "50-60", "70-80", ">80")
colnames(B) <- c("<30", "30-40", "40-50", "50-60", "70-80", ">80")
colnames(C) <- c("<30", "30-40", "40-50", "50-60", "70-80", ">80")
colnames(D) <- c("<30", "30-40", "40-50", "50-60", "70-80", ">80")

# Store matrices and scalar num in the list with the model results
spec_test_full_1$A <- A
spec_test_full_1$B <- B
spec_test_full_1$C <- C
spec_test_full_1$D <- D
spec_test_full_1$num <- num


# Drop lowp and highp from the dataframe
spec_test_full_1_filtered <- spec_test_full_1[!spec_test_full_1$term %in% c("lowp", "highp"),]

# Convert term to a factor and specify the order
spec_test_full_1_filtered$term <- factor(spec_test_full_1_filtered$term, levels = temp_bins)

new_row <- data.frame(term = "tMean_6070",  estimate = 0, std.error = 0, statistic = 0, p.value = 0, A = 0, B = 0, C = 0, D = 0, num = 0)

df1 <- spec_test_full_1_filtered[1:4, ]    # rows before
df2 <- spec_test_full_1_filtered[5:nrow(spec_test_full_1_filtered), ]  # rows after

df <- rbind(df1, new_row, df2)

bins_updated <- c("tMean_under30", "tMean_3040", "tMean_4050", "tMean_5060", "tMean_6070", "tMean_7080", "tMean_over80")
# Convert term to a factor and specify the order
df$term <- factor(df$term, levels = bins_updated)



library(ggplot2)
ggplot(df, aes(x = term, y = estimate, group = 1)) +
  geom_line(color = "steelblue") +
  geom_point(color = "steelblue") +
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error), width = 0.2) +
  theme_minimal() +
  labs(x = "Variable", y = "Estimated Coefficient") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


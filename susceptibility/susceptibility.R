library(tidyverse)
library(zoo)
library(reshape2)

# Read data
partially = read.csv("data/share-people-vaccinated-covid.csv")
fully = read.csv("data/share-people-fully-vaccinated-covid.csv")
infections_by_country = read.csv("data/infections.csv")
infections_uk = read.csv("data/infections_uk.csv")
seroprevalence = read.csv("data/seroprevalence.csv")

infections = rbind(infections_by_country, infections_uk) # Combine infections data

# Rename columns
partially = rename(partially, date = Day, area = Entity, partially_vaccinated = people_vaccinated_per_hundred)
fully = rename(fully, date = Day, area = Entity, fully_vaccinated = people_fully_vaccinated_per_hundred)
infections = rename(infections, area = areaName, cumulative_cases = cumCasesBySpecimenDate)

# Retain only relevant observations
uk = c("England", "Wales", "Scotland", "Northern Ireland", "United Kingdom")
partially = partially[partially$area %in% uk, names(partially) %in% c("area", "date", "partially_vaccinated")]
fully = fully[fully$area %in% uk, names(fully) %in% c("area", "date", "fully_vaccinated")]
infections = infections[, names(infections) %in% c("area", "date", "cumulative_cases")]

# Join datasets
data = merge(partially, fully, by=c("area", "date"), all=TRUE)
data = merge(data, infections, by=c("area", "date"), all=TRUE)

data$date = as.Date(data$date, format = "%Y-%m-%d") # Obtain dates from strings
seroprevalence$date = as.Date(seroprevalence$date, format="%d %B %Y")
data = merge(data, seroprevalence, by=c("area", "date"), all=TRUE)

# Reformat data
data$partially_vaccinated = data$partially_vaccinated / 100 # Convert to percentages
data$fully_vaccinated = data$fully_vaccinated / 100
data[order(data$area, data$date), ]

# Interpolate missing values
data = data %>% 
    group_by(area) %>% 
        mutate(partially_vaccinated = na.approx(partially_vaccinated, date, maxgap=7, na.rm=FALSE),
               fully_vaccinated = na.approx(fully_vaccinated, date, maxgap=7, na.rm=FALSE),
               seroprevalence = na.approx(seroprevalence, date, maxgap=7, na.rm=FALSE)) %>% 
    ungroup()
data[, names(data) != "seroprevalence"][is.na(data[, names(data) != "seroprevalence"])] = 0

# Compute percentage of people infected
percentage_of_population = function(country, n) {
    denominator = switch(country, "United Kingdom" = 67081234, "England" = 56550138, "Wales" = 3169586, "Scotland" = 5466000, "Northern Ireland" = 1895510)
    return(n / denominator)
}

data$infected = apply(data, 1, function(x) percentage_of_population(x["area"], strtoi(x["cumulative_cases"])))

# Estimate percentage of people immune
# p_imm: vaccine efficacy 21 days after first dose
# f_imm: vaccine efficacy 7 days after second dose
# i_imm: immunity conferred by past infection
percentage_immune = function(partially, fully, infected, p_imm=0.5, f_imm=0.95, i_imm=0.85, pi_imm=0.95, fi_imm=0.95) {
    # Compute percentages of people in each vaccinated*infected category
    partially = partially - fully # people who are partially but not fully vaccinated
    pi = partially * infected
    fi = fully * infected
    p = partially - pi 
    f = fully - fi 
    i = infected - pi - fi 

    return(p * p_imm + f * f_imm + i * i_imm + pi * pi_imm + fi * fi_imm)
}

data = data %>%
    group_by(area) %>%
    mutate(lagged_p = lag(partially_vaccinated, n=21, default=0), lagged_f = lag(fully_vaccinated, n=7, default=0))

data$est_immune = apply(data, 1, function(x) percentage_immune(as.double(x["lagged_p"]), as.double(x["lagged_f"]), as.double(x["infected"])))

# Clean up dataset
# data = data[, !(names(data) %in% c("lagged_p", "lagged_f", "cumulative_cases"))]
data = data %>% 
    mutate_if(is.double, function(x) round(x, digits=4))

# Plot data
plot = ggplot(data, aes(x=date, y=est_immune)) + geom_line(aes(color=area))
ggsave(file="plot_imm.svg", plot=plot, width=5, height=4)

plot2 = ggplot(data, aes(x=date, y=seroprevalence)) + geom_line(aes(color=area))
ggsave(file="plot_ser.svg", plot=plot2, width=5, height=4)

# Compare with model without infections and with underascertainment bias
uk_data = data[data$area == "United Kingdom", ]
uk_data$est_immune_wo_inf = apply(uk_data, 1, function(x) percentage_immune(as.double(x["lagged_p"]), as.double(x["lagged_f"]), as.double(x["infected"]), i_imm=0, pi_imm=0.5))

source("../under_ascertainment_bias.R")
owid_data = read.csv("data/owid-covid-data.csv")
owid_data$date = as.Date(owid_data$date, format="%Y-%m-%d")
bias = compute_underascertainment_bias(as.Date("01-30-2020", format="%m-%d-%Y"), "United Kingdom", owid_data)
uk_data = uk_data %>% 
    merge(bias, by="date", all=TRUE) %>% 
    rename(bias=value)
uk_data$bias[is.na(uk_data$bias)] = 1
uk_data$lagged_c = lag(uk_data$cumulative_cases, n=1, default=0)
uk_data$true_daily_cases = as.integer((uk_data$cumulative_cases - uk_data$lagged_c) / uk_data$bias)
uk_data$true_cases = cumsum(uk_data$true_daily_cases)
uk_data$true_infected = apply(uk_data, 1, function(x) percentage_of_population(x["area"], strtoi(x["true_cases"])))
uk_data$est_immune_w_bias = apply(uk_data, 1, function(x) percentage_immune(as.double(x["lagged_p"]), as.double(x["lagged_f"]), as.double(x["true_infected"])))

uk_data = uk_data[, names(uk_data) %in% c("date", "est_immune", "seroprevalence", "est_immune_wo_inf", "est_immune_w_bias")]
uk_data = melt(data = uk_data, id.vars = "date")
plot3 = ggplot(data = uk_data, aes(x = date, y = value, color = variable)) + geom_line()
ggsave(file="plot_cmp.svg", plot=plot3, width=5, height=4)
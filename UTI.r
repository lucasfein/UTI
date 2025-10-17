library(readxl)
library(dplyr)
library(tidyr)
library(forcats)

library(ggplot2)

data <- read_excel("UTI Model Bacterial and phage counts.xlsx", sheet = "MM02", range = "B2:E13", col_names = FALSE)

data <- data %>% rename(name = `...1`, R1 = `...2`, R2 = `...3`, R3 = `...4`) %>%
    pivot_longer(cols=R1:R3, names_to="technical", names_transform=function(name) as.integer(sub("R([0-9]+)", "\\1", name)), values_to="count") %>%
    separate_wider_regex(name, c(biological="V[0-9]+", " ", name=".*")) %>%
    mutate(name = fct_relevel(fct_recode(name, `UPEC 8923 3h` = "Cells + UPEC 3h", `UPEC 8923 + MM02 3h` = "Cells + UPEC + Phage 3h", `UPEC 8923 24h` = "Cells + UPEC 24h", `UPEC 8923 + MM02 24/19h` = "Cells + UPEC + Phage 24/19h"), "UPEC 8923 3h", "UPEC 8923 + MM02 3h", "UPEC 8923 24h", "UPEC 8923 + MM02 24/19h")) %>%
    relocate(biological, .after=name) %>%
    mutate(biological = as.integer(sub("V([0-9]+)", "\\1", biological)))

data <- data %>% group_by(name, biological) %>% summarize(count=mean(count), .groups="drop")

ggplot(data, aes(x=name, y=count)) + stat_summary(aes(y=count), fun=mean, geom="bar") + geom_jitter(aes(shape=factor(biological)), width=0.4, height=0) + expand_limits(y=0) + labs(y="CFU/ml", shape="biological\nreplicate") + theme_bw() + theme(axis.title.x=element_blank())

ggplot(data, aes(x=name, y=count)) + geom_boxplot() + geom_jitter(width=0.4, height=0) + expand_limits(y=0) + labs(y="CFU/ml") + theme_bw() + theme(axis.title.x=element_blank())

data <- read_excel("UTI Model Bacterial and phage counts.xlsx", sheet = "MM02", range = "B21:E36", col_names = FALSE)

data <- data %>% rename(name = `...1`, R1 = `...2`, R2 = `...3`, R3 = `...4`) %>%
    pivot_longer(cols=R1:R3, names_to="technical", names_transform=function(name) as.integer(sub("R([0-9]+)", "\\1", name)), values_to="count") %>%
    separate_wider_regex(name, c(biological="V[0-9]+", " ", name=".*")) %>%
    mutate(name = fct_relevel(fct_recode(name, `UPEC 8923 + MM02 3h` = "Cells + UPEC + Phage 3h", `MM02 only 3h` = "Cells only + Phage 3h", `UPEC 8923 + MM02 24/19h` = "Cells + UPEC + Phage 24/19h", `NC` = "Phage only 3 h (no cells, NC)"), "UPEC 8923 + MM02 3h", "MM02 only 3h", "UPEC 8923 + MM02 24/19h", "NC")) %>%
    relocate(biological, .after=name) %>%
    mutate(biological = as.integer(sub("V([0-9]+)", "\\1", biological)))

data <- data %>% group_by(name, biological) %>% summarize(count=mean(count), .groups="drop")

ggplot(data, aes(x=name, y=count)) + stat_summary(aes(y=count), fun=mean, geom="bar") + geom_jitter(aes(shape=factor(biological)), width=0.4, height=0) + expand_limits(y=0) + labs(y="PFU/ml", shape="biological\nreplicate") + theme_bw() + theme(axis.title.x=element_blank())

ggplot(data, aes(x=name, y=count)) + geom_boxplot() + geom_jitter(width=0.4, height=0) + expand_limits(y=0) + labs(y="PFU/ml") + theme_bw() + theme(axis.title.x=element_blank())

data <- read_excel("UTI Model Bacterial and phage counts.xlsx", sheet = "P00", range = "B2:E13", col_names = FALSE)

data <- data %>% rename(name = `...1`, R1 = `...2`, R2 = `...3`, R3 = `...4`) %>%
    pivot_longer(cols=R1:R3, names_to="technical", names_transform=function(name) as.integer(sub("R([0-9]+)", "\\1", name)), values_to="count") %>%
    separate_wider_regex(name, c(biological="V[0-9]+", " ", name=".*")) %>%
    mutate(name = fct_relevel(fct_recode(name, `UPEC 7958 3h` = "Cells + UPEC 3h", `UPEC 7958 + P00 3h` = "Cells + UPEC + Phage 3h", `UPEC 7958 24h` = "Cells + UPEC 24 h", `UPEC 7958 + P00 24/19h` = "Cells + UPEC + Phage 24/19 h"), "UPEC 7958 3h", "UPEC 7958 + P00 3h", "UPEC 7958 24h", "UPEC 7958 + P00 24/19h")) %>%
    relocate(biological, .after=name) %>%
    mutate(biological = as.integer(sub("V([0-9]+)", "\\1", biological)))

data <- data %>% group_by(name, biological) %>% summarize(count=mean(count), .groups="drop")

ggplot(data, aes(x=name, y=count)) + stat_summary(aes(y=count), fun=mean, geom="bar") + geom_jitter(aes(shape=factor(biological)), width=0.4, height=0) + expand_limits(y=0) + labs(y="CFU/ml", shape="biological\nreplicate") + theme_bw() + theme(axis.title.x=element_blank())

ggplot(data, aes(x=name, y=count)) + geom_boxplot() + geom_jitter(width=0.4, height=0) + expand_limits(y=0) + labs(y="CFU/ml") + theme_bw() + theme(axis.title.x=element_blank())

data <- read_excel("UTI Model Bacterial and phage counts.xlsx", sheet = "P00", range = "B23:E36", col_names = FALSE)

data <- data %>% rename(name = `...1`, R1 = `...2`, R2 = `...3`, R3 = `...4`) %>%
    pivot_longer(cols=R1:R3, names_to="technical", names_transform=function(name) as.integer(sub("R([0-9]+)", "\\1", name)), values_to="count") %>%
    separate_wider_regex(name, c(biological="V[0-9]+", " ", name=".*")) %>%
    mutate(name = fct_relevel(fct_recode(name, `UPEC 7958 + P00 3h` = "Cells + UPEC + Phage 3h", `P00 only 3h` = "Cells only + Phage 3h", `UPEC 7958 + P00 24/19h` = "Cells + UPEC + Phage 24/19 h", `NC` = "Phage only 3 h (no cells, NC)", `NC`="Phage only"), "UPEC 7958 + P00 3h", "P00 only 3h", "UPEC 7958 + P00 24/19h", "NC")) %>%
    relocate(biological, .after=name) %>%
    mutate(biological = as.integer(sub("V([0-9]+)", "\\1", biological)))

data <- data %>% group_by(name, biological) %>% summarize(count=mean(count), .groups="drop")

ggplot(data, aes(x=name, y=count)) + stat_summary(aes(y=count), fun=mean, geom="bar") + geom_jitter(aes(shape=factor(biological)), width=0.4, height=0) + expand_limits(y=0) + labs(y="PFU/ml", shape="biological\nreplicate") + theme_bw() + theme(axis.title.x=element_blank())

ggplot(data, aes(x=name, y=count)) + geom_boxplot() + geom_jitter(width=0.4, height=0) + expand_limits(y=0) + labs(y="PFU/ml") + theme_bw() + theme(axis.title.x=element_blank())

data <- read_excel("LDH.xlsx", sheet = "LDH 24h", range = "B2:E5", col_names = FALSE)

data <- data %>% rename(name = `...1`, `LDH V1` = `...2`, `LDH V2` = `...3`, `LDH V3` = `...4`) %>%
    pivot_longer(cols=`LDH V1`:`LDH V3`, names_to="replicate", names_transform=function(name) as.integer(sub("LDH V([0-9]+)", "\\1", name)), values_to="count") %>%
    mutate(name = fct_relevel(fct_recode(name, `UPEC 8923 24h` = "UPEC 8923 24h",  `UPEC 8923 + MM02 24h` = "UPEC 8923 + MM02 24h", `MM02 24h` = "MM02 24h", `Cells only 24h` = "Cells only 24h"), "UPEC 8923 24h", "UPEC 8923 + MM02 24h", "MM02 24h", "Cells only 24h"))

ggplot(data, aes(x=name, y=count)) + geom_boxplot() + geom_jitter(width=0.4, height=0) + expand_limits(y=0) + labs(y="LDH") + theme_bw() + theme(axis.title.x=element_blank())

data <- read_excel("LDH.xlsx", sheet = "LDH 24h", range = "B7:E10", col_names = FALSE)

data <- data %>% rename(name = `...1`, `LDH V1` = `...2`, `LDH V2` = `...3`, `LDH V3` = `...4`) %>%
    pivot_longer(cols=`LDH V1`:`LDH V3`, names_to="replicate", names_transform=function(name) as.integer(sub("LDH V([0-9]+)", "\\1", name)), values_to="count") %>%
    mutate(name = fct_relevel(fct_recode(name, `UPEC 7958 24h` = "UPEC 7958", `UPEC 7958 + G10400 24h` = "UPEC 7958 + G10400", `G10400 24h` = "G10400 24h", `Cells only 24h` = "Cells only 24h"), "UPEC 7958 24h", "UPEC 7958 + G10400 24h", "G10400 24h", "Cells only 24h"))

ggplot(data, aes(x=name, y=count)) + geom_boxplot() + geom_jitter(width=0.4, height=0) + expand_limits(y=0) + labs(y="LDH") + theme_bw() + theme(axis.title.x=element_blank())

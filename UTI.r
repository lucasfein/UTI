library(readxl)
library(dplyr)
library(tidyr)
library(forcats)

library(ggplot2)

library(dunn.test)

library(ggsignif)
library(scales)
library(cowplot)

pdf(width=14)

data1 <- read_excel("data/UTI Model Bacterial and phage counts.xlsx", sheet = "MM02", range = "B2:E13", col_names = FALSE)

data1 <- data1 %>% rename(name = `...1`, R1 = `...2`, R2 = `...3`, R3 = `...4`) %>%
    pivot_longer(cols=R1:R3, names_to="technical", names_transform=function(name) as.integer(sub("R([0-9]+)", "\\1", name)), values_to="value") %>%
    separate_wider_regex(name, c(biological="V[0-9]+", " ", name=".*")) %>%
    mutate(name = fct_relevel(fct_recode(name, `UPEC 8923 3h` = "Cells + UPEC 3h", `UPEC 8923 + MM02 3h` = "Cells + UPEC + Phage 3h", `UPEC 8923 24h` = "Cells + UPEC 24h", `UPEC 8923 + MM02 24/19h` = "Cells + UPEC + Phage 24/19h"), "UPEC 8923 3h", "UPEC 8923 + MM02 3h", "UPEC 8923 24h", "UPEC 8923 + MM02 24/19h")) %>%
    relocate(biological, .after=name) %>%
    mutate(biological = as.integer(sub("V([0-9]+)", "\\1", biological)))

data1 <- data1 %>% group_by(name, biological) %>% summarize(value=mean(value), .groups="drop")

data2 <- read_excel("data/UTI Model Bacterial and phage counts.xlsx", sheet = "MM02", range = "B21:E36", col_names = FALSE)

data2 <- data2 %>% rename(name = `...1`, R1 = `...2`, R2 = `...3`, R3 = `...4`) %>%
    pivot_longer(cols=R1:R3, names_to="technical", names_transform=function(name) as.integer(sub("R([0-9]+)", "\\1", name)), values_to="value") %>%
    separate_wider_regex(name, c(biological="V[0-9]+", " ", name=".*")) %>%
    mutate(name = fct_relevel(fct_recode(name, `UPEC 8923 + MM02 3h` = "Cells + UPEC + Phage 3h", `MM02 only 3h` = "Cells only + Phage 3h", `UPEC 8923 + MM02 24/19h` = "Cells + UPEC + Phage 24/19h", `NC` = "Phage only 3 h (no cells, NC)"), "UPEC 8923 + MM02 3h", "MM02 only 3h", "UPEC 8923 + MM02 24/19h", "NC")) %>%
    relocate(biological, .after=name) %>%
    mutate(biological = as.integer(sub("V([0-9]+)", "\\1", biological)))

data2 <- data2 %>% group_by(name, biological) %>% summarize(value=mean(value), .groups="drop")

data3 <- read_excel("data/UTI Model Bacterial and phage counts.xlsx", sheet = "P00", range = "B2:E13", col_names = FALSE)

data3 <- data3 %>% rename(name = `...1`, R1 = `...2`, R2 = `...3`, R3 = `...4`) %>%
    pivot_longer(cols=R1:R3, names_to="technical", names_transform=function(name) as.integer(sub("R([0-9]+)", "\\1", name)), values_to="value") %>%
    separate_wider_regex(name, c(biological="V[0-9]+", " ", name=".*")) %>%
    mutate(name = fct_relevel(fct_recode(name, `UPEC 7958 3h` = "Cells + UPEC 3h", `UPEC 7958 + G10400 3h` = "Cells + UPEC + Phage 3h", `UPEC 7958 24h` = "Cells + UPEC 24 h", `UPEC 7958 + G10400 24/19h` = "Cells + UPEC + Phage 24/19 h"), "UPEC 7958 3h", "UPEC 7958 + G10400 3h", "UPEC 7958 24h", "UPEC 7958 + G10400 24/19h")) %>%
    relocate(biological, .after=name) %>%
    mutate(biological = as.integer(sub("V([0-9]+)", "\\1", biological)))

data3 <- data3 %>% group_by(name, biological) %>% summarize(value=mean(value), .groups="drop")

data4 <- read_excel("data/UTI Model Bacterial and phage counts.xlsx", sheet = "P00", range = "B23:E36", col_names = FALSE)

data4 <- data4 %>% rename(name = `...1`, R1 = `...2`, R2 = `...3`, R3 = `...4`) %>%
    pivot_longer(cols=R1:R3, names_to="technical", names_transform=function(name) as.integer(sub("R([0-9]+)", "\\1", name)), values_to="value") %>%
    separate_wider_regex(name, c(biological="V[0-9]+", " ", name=".*")) %>%
    mutate(name = fct_relevel(fct_recode(name, `UPEC 7958 + G10400 3h` = "Cells + UPEC + Phage 3h", `G10400 only 3h` = "Cells only + Phage 3h", `UPEC 7958 + G10400 24/19h` = "Cells + UPEC + Phage 24/19 h", `NC` = "Phage only 3 h (no cells, NC)", `NC`="Phage only"), "UPEC 7958 + G10400 3h", "G10400 only 3h", "UPEC 7958 + G10400 24/19h", "NC")) %>%
    relocate(biological, .after=name) %>%
    mutate(biological = as.integer(sub("V([0-9]+)", "\\1", biological)))

data4 <- data4 %>% group_by(name, biological) %>% summarize(value=mean(value), .groups="drop")

data5 <- read_excel("data/LDH.xlsx", sheet = "LDH 24h", range = "B2:E5", col_names = FALSE)

data5 <- data5 %>% rename(name = `...1`, `LDH V1` = `...2`, `LDH V2` = `...3`, `LDH V3` = `...4`) %>%
    pivot_longer(cols=`LDH V1`:`LDH V3`, names_to="replicate", names_transform=function(name) as.integer(sub("LDH V([0-9]+)", "\\1", name)), values_to="value") %>%
    mutate(name = fct_relevel(fct_recode(name, `UPEC 8923 24h` = "UPEC 8923 24h",  `UPEC 8923 + MM02 24h` = "UPEC 8923 + MM02 24h", `MM02 24h` = "MM02 24h", `Cells only 24h` = "Cells only 24h"), "UPEC 8923 24h", "UPEC 8923 + MM02 24h", "MM02 24h", "Cells only 24h"))

data6 <- read_excel("data/LDH.xlsx", sheet = "LDH 24h", range = "B7:E10", col_names = FALSE)

data6 <- data6 %>% rename(name = `...1`, `LDH V1` = `...2`, `LDH V2` = `...3`, `LDH V3` = `...4`) %>%
    pivot_longer(cols=`LDH V1`:`LDH V3`, names_to="replicate", names_transform=function(name) as.integer(sub("LDH V([0-9]+)", "\\1", name)), values_to="value") %>%
    mutate(name = fct_relevel(fct_recode(name, `UPEC 7958 24h` = "UPEC 7958", `UPEC 7958 + G10400 24h` = "UPEC 7958 + G10400", `G10400 24h` = "G10400 24h", `Cells only 24h` = "Cells only 24h"), "UPEC 7958 24h", "UPEC 7958 + G10400 24h", "G10400 24h", "Cells only 24h"))

if (p.adjust(c(
        kruskal.test(value ~ name, data = data1)$p.value,
        kruskal.test(value ~ name, data = data3)$p.value),
        method="holm")[1] < 0.05) {
    dunn.test(data1 %>% pull(value), data1 %>% pull(name), method="holm", table=FALSE, list=TRUE)
}
if (p.adjust(c(
        kruskal.test(value ~ name, data = data2)$p.value,
        kruskal.test(value ~ name, data = data4)$p.value),
        method="holm")[1] < 0.05) {
    dunn.test(data2 %>% pull(value), data2 %>% pull(name), method="holm", table=FALSE, list=TRUE)
}
if (p.adjust(c(
        kruskal.test(value ~ name, data = data1)$p.value,
        kruskal.test(value ~ name, data = data3)$p.value),
        method="holm")[2] < 0.05) {
    dunn.test(data3 %>% pull(value), data3 %>% pull(name), method="holm", table=FALSE, list=TRUE)
}
if (p.adjust(c(
        kruskal.test(value ~ name, data = data2)$p.value,
        kruskal.test(value ~ name, data = data4)$p.value),
        method="holm")[2] < 0.05) {
    dunn.test(data4 %>% pull(value), data4 %>% pull(name), method="holm", table=FALSE, list=TRUE)
}
if (p.adjust(c(
        kruskal.test(value ~ name, data = data5)$p.value,
        kruskal.test(value ~ name, data = data6)$p.value),
        method="holm")[1] < 0.05) {
    dunn.test(data5 %>% pull(value), data5 %>% pull(name), method="holm", table=FALSE, list=TRUE)
}
if (p.adjust(c(
        kruskal.test(value ~ name, data = data5)$p.value,
        kruskal.test(value ~ name, data = data6)$p.value),
        method="holm")[2] < 0.05) {
    dunn.test(data6 %>% pull(value), data6 %>% pull(name), method="holm", table=FALSE, list=TRUE)
}

plot_grid(
    ggplot(data1, aes(x=name, y=value)) + geom_boxplot() + geom_jitter(width=0.35, height=0) + expand_limits(y=0) + labs(y="CFU/mL") + theme_bw() + theme(axis.title.x=element_blank()),
    ggplot(data3, aes(x=name, y=value)) + geom_boxplot() + geom_jitter(width=0.35, height=0) + expand_limits(y=0) + labs(y="CFU/mL") + theme_bw() + theme(axis.title.x=element_blank()),
    labels=c("A", "B"), align="v")

plot_grid(
    ggplot(data2, aes(x=name, y=value)) + geom_boxplot() + geom_jitter(width=0.35, height=0) + expand_limits(y=0) + scale_y_continuous(labels=scales::label_number_auto()) + labs(y="PFU/mL") + theme_bw() + theme(axis.title.x=element_blank()) + suppressWarnings(geom_signif(data = data.frame(start="UPEC 8923 + MM02 3h", end="NC", label="p < 0.01", y=175000), aes(xmin = start, xmax = end, annotations = label, y_position=y), manual = TRUE, inherit.aes = FALSE)),
    ggplot(data4, aes(x=name, y=value)) + geom_boxplot() + geom_jitter(width=0.35, height=0) + expand_limits(y=0) + scale_y_continuous(labels=scales::label_number_auto()) + labs(y="PFU/mL") + theme_bw() + theme(axis.title.x=element_blank()) + suppressWarnings(geom_signif(data = data.frame(start="UPEC 7958 + G10400 3h", end="NC", label="p < 0.01", y=3.5e05), aes(xmin = start, xmax = end, annotations = label, y_position=y), manual = TRUE, inherit.aes = FALSE)),
    labels=c("A", "B"), align="v")

plot_grid(
    ggplot(data5, aes(x=name, y=value)) + geom_boxplot() + geom_jitter(width=0.35, height=0) + expand_limits(y=0) + labs(y="LDH (U/L)") + theme_bw() + theme(axis.title.x=element_blank()) + suppressWarnings(geom_signif(data = data.frame(start="MM02 24h", end="Cells only 24h", label="p = 0.01", y=150), aes(xmin = start, xmax = end, annotations = label, y_position=y), manual = TRUE, inherit.aes = FALSE)),
    ggplot(data6, aes(x=name, y=value)) + geom_boxplot() + geom_jitter(width=0.35, height=0) + expand_limits(y=0) + labs(y="LDH (U/L)") + theme_bw() + theme(axis.title.x=element_blank()) + suppressWarnings(geom_signif(data = data.frame(start="UPEC 7958 24h", end="UPEC 7958 + G10400 24h", label="p = 0.01", y=300), aes(xmin = start, xmax = end, annotations = label, y_position=y), manual = TRUE, inherit.aes = FALSE)),
    labels=c("A", "B"), align="v")

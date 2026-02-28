library(readxl)

library(dplyr)
library(tidyr)

library(forcats)

library(ggplot2)

library(scales)

library(ggsignif)
library(patchwork)

library(dunn.test)

library(ggokabeito)

library(coin)

data1 <- read_excel(
    "data/LDH.xlsx",
    sheet = "LDH 24h",
    range = "B2:E5",
    col_names = FALSE,
    .name_repair = "unique_quiet"
) %>%
    rename(
        name = `...1`,
        `LDH V1` = `...2`,
        `LDH V2` = `...3`,
        `LDH V3` = `...4`
    ) %>%
    pivot_longer(
        cols = matches("LDH V[0-9]"),
        names_to = "replicate",
        values_to = "value"
    ) %>%
    mutate(
        name = fct_relevel(
            name,
            "UPEC 8923 24h",
            "UPEC 8923 + MM02 24h",
            "MM02 24h",
            "Cells only 24h"
        )
    )

data1 %>% count(name)

data1 %>% group_by(name) %>% summarise(SD = sd(value))

dunn.test(
    data1 %>% pull(value),
    data1 %>% pull(name),
    method = "holm",
    list = TRUE,
    altp = TRUE
)

data2 <- read_excel(
    "data/LDH.xlsx",
    sheet = "LDH 24h",
    range = "B7:E10",
    col_names = FALSE,
    .name_repair = "unique_quiet"
) %>%
    rename(
        name = `...1`,
        `LDH V1` = `...2`,
        `LDH V2` = `...3`,
        `LDH V3` = `...4`
    ) %>%
    pivot_longer(
        cols = matches("LDH V[0-9]"),
        names_to = "replicate",
        values_to = "value"
    ) %>%
    mutate(
        name = fct_relevel(
            name,
            "UPEC 7958",
            "UPEC 7958 + G10400",
            "G10400 24h",
            "Cells only 24h"
        )
    )

data2 %>% count(name)

data2 %>% group_by(name) %>% summarise(SD = sd(value))

dunn.test(
    data2 %>% pull(value),
    data2 %>% pull(name),
    method = "holm",
    list = TRUE,
    altp = TRUE
)

(ggplot(data1, aes(name, value, fill = name)) +
    stat_summary(fun = "mean", geom = "bar", show.legend = FALSE) +
    stat_summary(
        fun.data = "mean_cl_normal",
        geom = "errorbar",
        linewidth = 0.25,
        show.legend = FALSE
    ) +
    geom_jitter(position = position_jitter(height = 0, seed = 5)) +
    ylab("LDH (U/L)") +
    scale_x_discrete(
        labels = c(
            "UPEC 8923 24h" = "UPEC 8923 24h",
            "UPEC 8923 + MM02 24h" = expression(
                "UPEC 8923 +" ~ Phi ~ "MM02 24h"
            ),
            "MM02 24h" = expression(Phi ~ "MM02 24h"),
            "Cells only 24h" = "Cells only 24h"
        )
    ) +
    scale_y_continuous(
        labels = label_number(),
        expand = expansion(mult = c(0, 0.1)),
        breaks = seq(0, 1000, 200)
    ) +
    theme_bw(base_size = 10) +
    theme(axis.title.x = element_blank(), legend.position = "none") +
    geom_signif(
        comparisons = list(c("UPEC 8923 + MM02 24h", "UPEC 8923 24h")),
        annotations = "*",
        size = 0.25,
        textsize = 10 * 0.8 / .pt,
        y_position = 1000
    ) +
    scale_fill_okabe_ito()) /
    (ggplot(data2, aes(name, value, fill = name)) +
        stat_summary(fun = "mean", geom = "bar", show.legend = FALSE) +
        stat_summary(
            fun.data = "mean_cl_normal",
            geom = "errorbar",
            linewidth = 0.25,
            show.legend = FALSE
        ) +
        geom_jitter(position = position_jitter(height = 0, seed = 5)) +
        ylab("LDH (U/L)") +
        scale_x_discrete(
            labels = c(
                "UPEC 7958" = "UPEC 7958 24h",
                "UPEC 7958 + G10400" = expression(
                    "UPEC 7958 +" ~ Phi ~ "G10400 24h"
                ),
                "G10400 24h" = expression(Phi ~ "G10400 24h"),
                "Cells only 24h" = "Cells only 24h"
            )
        ) +
        scale_y_continuous(
            labels = label_number(),
            expand = expansion(mult = c(0, 0.1)),
            breaks = seq(0, 500, 100)
        ) +
        theme_bw(base_size = 10) +
        theme(axis.title.x = element_blank(), legend.position = "none") +
        geom_signif(
            comparisons = list(c("UPEC 7958", "UPEC 7958 + G10400")),
            annotations = "*",
            size = 0.25,
            textsize = 10 * 0.8 / .pt,
            y_position = 500
        ) +
        scale_fill_okabe_ito()) +
    plot_annotation(tag_levels = "A")

data3 <- bind_rows(
    read_excel(
        "data/UTI Model Bacterial and phage counts.16.1.26xlsx.xlsx",
        sheet = "MM02",
        range = "B33:N40",
        col_names = FALSE,
        .name_repair = "unique_quiet"
    ) %>%
        select(1, 13) %>%
        rename(name = `...1`, value = `...13`),
    read_excel(
        "data/UTI Model Bacterial and phage counts.16.1.26xlsx.xlsx",
        sheet = "MM02",
        range = "B50:N55",
        col_names = FALSE,
        .name_repair = "unique_quiet"
    ) %>%
        select(1, 13) %>%
        rename(name = `...1`, value = `...13`)
) %>%
    mutate(
        name = fct_relevel(
            name,
            "MM02 + UPEC 7958 3h",
            "MM02 cells only 3h",
            "MM02 + UPEC 7958 24h",
            "MM02 cells only 24h"
        )
    )

data3 %>%
    filter(name == "MM02 + UPEC 7958 3h" | name == "MM02 cells only 3h") %>%
    count(name)

data3 %>%
    filter(name == "MM02 + UPEC 7958 3h" | name == "MM02 cells only 3h") %>%
    group_by(name) %>%
    summarise(SD = sd(value))

result1 <- t.test(
    I(log10(value)) ~ fct_rev(name),
    data = data3 %>%
        filter(name == "MM02 + UPEC 7958 3h" | name == "MM02 cells only 3h")
)

round(unname(10^-diff(result1$estimate)), 2)

round(as.numeric(10^result1$conf.int), 2)

data3 %>%
    filter(name == "MM02 + UPEC 7958 24h" | name == "MM02 cells only 24h") %>%
    count(name)

data3 %>%
    filter(name == "MM02 + UPEC 7958 24h" | name == "MM02 cells only 24h") %>%
    group_by(name) %>%
    summarise(SD = sd(value))


result2 <- t.test(
    I(log10(value)) ~ fct_rev(name),
    data = data3 %>%
        filter(name == "MM02 + UPEC 7958 24h" | name == "MM02 cells only 24h")
)

round(unname(10^-diff(result2$estimate)), 2)

round(as.numeric(10^result2$conf.int), 2)

p.adjust(
    c(
        result1$p.value,
        result2$p.value
    ),
    method = "holm"
)

data4 <- bind_rows(
    read_excel(
        "data/UTI Model Bacterial and phage counts.16.1.26xlsx.xlsx",
        sheet = "P00",
        range = "B23:N30",
        col_names = FALSE,
        .name_repair = "unique_quiet"
    ) %>%
        select(1, 13) %>%
        rename(name = `...1`, value = `...13`),
    read_excel(
        "data/UTI Model Bacterial and phage counts.16.1.26xlsx.xlsx",
        sheet = "P00",
        range = "B39:N44",
        col_names = FALSE,
        .name_repair = "unique_quiet"
    ) %>%
        select(1, 13) %>%
        rename(name = `...1`, value = `...13`)
) %>%
    mutate(
        name = fct_relevel(
            name,
            "G10400 + UPEC 8923 3h",
            "G10400 cells only 3h",
            "G10400 + UPEC 8923 24h",
            "G10400 cells only 24h"
        )
    )

data4 %>%
    filter(
        name == "G10400 + UPEC 8923 3h" | name == "G10400 cells only 3h"
    ) %>%
    count(name)

data4 %>%
    filter(
        name == "G10400 + UPEC 8923 3h" | name == "G10400 cells only 3h"
    ) %>%
    group_by(name) %>%
    summarise(SD = sd(value))

result1 <- t.test(
    I(log10(value)) ~ fct_rev(name),
    data = data4 %>%
        filter(
            name == "G10400 + UPEC 8923 3h" | name == "G10400 cells only 3h"
        )
)

round(unname(10^-diff(result1$estimate)), 2)

round(as.numeric(10^result1$conf.int), 2)

data4 %>%
    filter(
        name == "G10400 + UPEC 8923 24h" | name == "G10400 cells only 24h"
    ) %>%
    count(name)

data4 %>%
    filter(
        name == "G10400 + UPEC 8923 24h" | name == "G10400 cells only 24h"
    ) %>%
    group_by(name) %>%
    summarise(SD = sd(value))

result2 <- t.test(
    I(log10(value)) ~ fct_rev(name),
    data = data4 %>%
        filter(
            name == "G10400 + UPEC 8923 24h" | name == "G10400 cells only 24h"
        )
)

round(unname(10^-diff(result2$estimate)), 2)

round(as.numeric(10^result2$conf.int), 2)

p.adjust(
    c(
        result1$p.value,
        result2$p.value
    ),
    method = "holm"
)

data3 <- data3 %>%
    mutate(
        group = case_when(
            grepl("3h", name) ~ "1",
            grepl("24h", name) ~ "2"
        )
    )

data4 <- data4 %>%
    mutate(
        group = case_when(
            grepl("3h", name) ~ "1",
            grepl("24h", name) ~ "2"
        )
    )

(ggplot(data3, aes(name, value, fill = grepl("cells only", name))) +
    stat_summary(fun = "mean", geom = "bar", show.legend = FALSE) +
    stat_summary(
        fun.data = "mean_cl_normal",
        geom = "errorbar",
        linewidth = 0.25,
        show.legend = FALSE
    ) +
    geom_jitter(position = position_jitter(height = 0, seed = 1)) +
    facet_grid(
        cols = vars(group),
        scales = "free_x",
        labeller = as_labeller(c(`1` = "3 h", `2` = "24 h"))
    ) +
    expand_limits(y = 0) +
    ylab("PFU/ml") +
    scale_x_discrete(
        labels = c(
            "MM02 + UPEC 7958 3h" = expression(
                "UPEC 7958 +" ~ Phi ~ "MM02"
            ),
            "MM02 cells only 3h" = expression(Phi ~ "MM02"),
            "MM02 + UPEC 7958 24h" = expression(
                "UPEC 7958 +" ~ Phi ~ "MM02"
            ),
            "MM02 cells only 24h" = expression(Phi ~ "MM02")
        )
    ) +
    scale_y_continuous(
        transform = transform_pseudo_log(base = 10),
        breaks = c(0, 10^3, 10^6),
        labels = label_log(base = 10),
        expand = expansion(mult = c(0, 0.1))
    ) +
    theme_bw(base_size = 10) +
    theme(axis.title.x = element_blank(), legend.position = "none") +
    suppressWarnings(geom_signif(
        data = data.frame(
            group = as.character(1:2),
            start = c("MM02 + UPEC 7958 3h", "MM02 + UPEC 7958 24h"),
            end = c("MM02 cells only 3h", "MM02 cells only 24h"),
            label = c("**", "×0.20"),
            y = c(7, 7)
        ),
        aes(xmin = start, xmax = end, annotations = label, y_position = y),
        manual = TRUE,
        inherit.aes = FALSE,
        size = 0.25,
        textsize = 10 * 0.8 / .pt
    )) +
    scale_fill_okabe_ito()) /
    (ggplot(data4, aes(name, value, fill = grepl("cells only", name))) +
        stat_summary(fun = "mean", geom = "bar", show.legend = FALSE) +
        stat_summary(
            fun.data = "mean_cl_normal",
            geom = "errorbar",
            linewidth = 0.25,
            show.legend = FALSE
        ) +
        geom_jitter(position = position_jitter(height = 0, seed = 3)) +
        facet_grid(
            cols = vars(group),
            scales = "free_x",
            labeller = as_labeller(c(`1` = "3 h", `2` = "24 h"))
        ) +
        expand_limits(y = 0) +
        ylab("PFU/ml") +
        scale_x_discrete(
            labels = c(
                "G10400 + UPEC 8923 3h" = expression(
                    "UPEC 8923 +" ~ Phi ~ "G10400"
                ),
                "G10400 cells only 3h" = expression(Phi ~ "G10400"),
                "G10400 + UPEC 8923 24h" = expression(
                    "UPEC 8923 +" ~ Phi ~ "G10400"
                ),
                "G10400 cells only 24h" = expression(Phi ~ "G10400")
            )
        ) +
        scale_y_continuous(
            transform = transform_pseudo_log(base = 10),
            breaks = c(0, 10^3, 10^6),
            labels = label_log(base = 10),
            expand = expansion(mult = c(0, 0.1))
        ) +
        theme_bw(base_size = 10) +
        theme(axis.title.x = element_blank(), legend.position = "none") +
        suppressWarnings(geom_signif(
            data = data.frame(
                group = as.character(1:2),
                start = c("G10400 + UPEC 8923 3h", "G10400 + UPEC 8923 24h"),
                end = c("G10400 cells only 3h", "G10400 cells only 24h"),
                label = c("***", "***"),
                y = c(7, 7)
            ),
            aes(xmin = start, xmax = end, annotations = label, y_position = y),
            manual = TRUE,
            inherit.aes = FALSE,
            size = 0.25,
            textsize = 10 * 0.8 / .pt
        )) +
        scale_fill_okabe_ito()) +
    plot_annotation(tag_levels = "A")

data5 <- bind_rows(
    read_excel(
        "data/UTI Model Bacterial and phage counts.16.1.26xlsx.xlsx",
        sheet = "MM02",
        range = "B15:G22",
        col_names = FALSE,
        .name_repair = "unique_quiet"
    ) %>%
        select(1, 6) %>%
        rename(name = `...1`, value = `...6`),
    read_excel(
        "data/UTI Model Bacterial and phage counts.16.1.26xlsx.xlsx",
        sheet = "MM02",
        range = "B24:G29",
        col_names = FALSE,
        .name_repair = "unique_quiet"
    ) %>%
        select(1, 6) %>%
        rename(name = `...1`, value = `...6`)
) %>%
    mutate(
        name = fct_relevel(
            sub("V[0-9] +(.+)", "\\1", name),
            "UPEC 1h",
            "MM02 after 1,5h",
            "UPEC 0,5h",
            "MM02 after 1h"
        )
    )

data5 %>%
    filter(name == "UPEC 1h" | name == "MM02 after 1,5h") %>%
    count(name)

data5 %>%
    filter(name == "UPEC 1h" | name == "MM02 after 1,5h") %>%
    group_by(name) %>%
    summarise(SD = sd(value))

result1 <- t.test(
    I(log10(value)) ~ fct_rev(name),
    data = data5 %>%
        filter(name == "UPEC 1h" | name == "MM02 after 1,5h")
)

round(unname(10^-diff(result1$estimate)), 2)

round(as.numeric(10^result1$conf.int), 2)

data5 %>%
    filter(name == "UPEC 0,5h" | name == "MM02 after 1h") %>%
    count(name)

data5 %>%
    filter(name == "UPEC 0,5h" | name == "MM02 after 1h") %>%
    group_by(name) %>%
    summarise(SD = sd(value))

result2 <- t.test(
    I(log10(value)) ~ fct_rev(name),
    data = data5 %>%
        filter(name == "UPEC 0,5h" | name == "MM02 after 1h")
)

round(unname(10^-diff(result2$estimate)), 2)

round(as.numeric(10^result2$conf.int), 2)

p.adjust(
    c(
        result1$p.value,
        result2$p.value
    ),
    method = "holm"
)

data5 <- data5 %>%
    mutate(
        group = case_when(
            grepl("UPEC 1h|MM02 after 1,5h", name) ~ "1",
            grepl("UPEC 0,5h|MM02 after 1h", name) ~ "2"
        )
    )

ggplot(data5, aes(name, value, fill = grepl("MM02", name))) +
    stat_summary(fun = "mean", geom = "bar", show.legend = FALSE) +
    stat_summary(
        fun.data = "mean_cl_normal",
        geom = "errorbar",
        linewidth = 0.25,
        show.legend = FALSE
    ) +
    geom_jitter(position = position_jitter(height = 0, seed = 2)) +
    facet_grid(
        cols = vars(group),
        scales = "free_x",
        labeller = as_labeller(
            c(
                `1` = "'1/24 h (' * varphi ~ '1.5/24 h)'",
                `2` = "'0.5/24 h (' * varphi ~ '1/24 h)'"
            ),
            default = label_parsed
        )
    ) +
    ylab("CFU/ml") +
    scale_x_discrete(
        labels = c(
            "UPEC 1h" = "UPEC 8923",
            "MM02 after 1,5h" = expression(
                "UPEC 8923 +" ~ Phi ~ "MM02 1.5h PI"
            ),
            "UPEC 0,5h" = "UPEC 8923",
            "MM02 after 1h" = expression(
                "UPEC 8923 +" ~ Phi ~ "MM02 1h PI"
            )
        )
    ) +
    scale_y_continuous(
        transform = transform_pseudo_log(base = 10),
        breaks = c(0, 10^1, 10^2, 10^3),
        labels = label_log(base = 10),
        expand = expansion(mult = c(0, 0.05))
    ) +
    theme_bw(base_size = 10) +
    theme(axis.title.x = element_blank(), legend.position = "none") +
    suppressWarnings(geom_signif(
        data = data.frame(
            group = as.character(1:2),
            start = c("UPEC 1h", "UPEC 0,5h"),
            end = c("MM02 after 1,5h", "MM02 after 1h"),
            label = c("×0.74", "×0.64"),
            y = c(3, 3)
        ),
        aes(xmin = start, xmax = end, annotations = label, y_position = y),
        manual = TRUE,
        inherit.aes = FALSE,
        size = 0.25,
        textsize = 10 * 0.8 / .pt
    )) +
    scale_fill_okabe_ito()

data6 <- read_excel(
    "data/UTI Model Bacterial and phage counts.16.1.26xlsx.xlsx",
    sheet = "MM02",
    range = "B2:G13",
    col_names = FALSE,
    .name_repair = "unique_quiet"
) %>%
    select(1, 6) %>%
    rename(name = `...1`, value = `...6`) %>%
    mutate(
        name = fct_relevel(
            sub("V[0-9] +(.+)", "\\1", name),
            "Cells + UPEC 3h",
            "Cells + UPEC + Phage 3h",
            "Cells + UPEC 24h",
            "Cells + UPEC + Phage 24/19h"
        )
    )

data6 %>%
    filter(name == "Cells + UPEC 3h" | name == "Cells + UPEC + Phage 3h") %>%
    count(name)

data6 %>%
    filter(name == "Cells + UPEC 3h" | name == "Cells + UPEC + Phage 3h") %>%
    group_by(name) %>%
    summarise(SD = sd(value))

data6 %>%
    filter(
        name == "Cells + UPEC 24h" | name == "Cells + UPEC + Phage 24/19h"
    ) %>%
    count(name)

data6 %>%
    filter(
        name == "Cells + UPEC 24h" | name == "Cells + UPEC + Phage 24/19h"
    ) %>%
    group_by(name) %>%
    summarise(SD = sd(value))

result2 <- t.test(
    I(log10(value)) ~ fct_rev(name),
    data = data6 %>%
        filter(
            name == "Cells + UPEC 24h" | name == "Cells + UPEC + Phage 24/19h"
        )
)

round(unname(10^-diff(result2$estimate)), 2)

round(as.numeric(10^result2$conf.int), 2)

result2$p.value

data7 <- read_excel(
    "data/UTI Model Bacterial and phage counts.16.1.26xlsx.xlsx",
    sheet = "P00",
    range = "B2:G13",
    col_names = FALSE,
    .name_repair = "unique_quiet"
) %>%
    select(1, 6) %>%
    rename(name = `...1`, value = `...6`) %>%
    mutate(
        name = fct_relevel(
            sub("V[0-9] +(.+)", "\\1", name),
            "Cells + UPEC 3h",
            "Cells + UPEC + Phage 3h",
            "Cells + UPEC 24 h",
            "Cells + UPEC + Phage 24/19 h"
        )
    )

data7 %>%
    filter(name == "Cells + UPEC 3h" | name == "Cells + UPEC + Phage 3h") %>%
    count(name)

data7 %>%
    filter(name == "Cells + UPEC 3h" | name == "Cells + UPEC + Phage 3h") %>%
    group_by(name) %>%
    summarise(SD = sd(value))

data7 %>%
    filter(
        name == "Cells + UPEC 24 h" | name == "Cells + UPEC + Phage 24/19 h"
    ) %>%
    count(name)

data7 %>%
    filter(
        name == "Cells + UPEC 24 h" | name == "Cells + UPEC + Phage 24/19 h"
    ) %>%
    group_by(name) %>%
    summarise(SD = sd(value))

result2 <- t.test(
    I(log10(value)) ~ fct_rev(name),
    data = data7 %>%
        filter(
            name == "Cells + UPEC 24 h" | name == "Cells + UPEC + Phage 24/19 h"
        )
)

round(unname(10^-diff(result2$estimate)), 2)

round(as.numeric(10^result2$conf.int), 2)

result2$p.value

data6 <- data6 %>%
    mutate(
        group = case_when(
            grepl("3h", name) ~ "1",
            grepl("24h|24/19h", name) ~ "2"
        )
    )

data7 <- data7 %>%
    mutate(
        group = case_when(
            grepl("3h", name) ~ "1",
            grepl("24 h|24/19 h", name) ~ "2"
        )
    )

(ggplot(data6, aes(name, value, fill = grepl("Phage", name))) +
    stat_summary(fun = "mean", geom = "bar", show.legend = FALSE) +
    stat_summary(
        fun.data = "mean_cl_normal",
        geom = "errorbar",
        linewidth = 0.25,
        show.legend = FALSE
    ) +
    geom_jitter(position = position_jitter(height = 0, seed = 5)) +
    facet_grid(
        cols = vars(group),
        scales = "free_x",
        labeller = as_labeller(
            c(
                `1` = "'3/5 h (' * varphi ~ '3/5 h)'",
                `2` = "'3/24 h (' * varphi ~ '5/24 h)'"
            ),
            default = label_parsed
        )
    ) +
    ylab("CFU/ml") +
    scale_x_discrete(
        labels = c(
            "Cells + UPEC 3h" = "UPEC 8923",
            "Cells + UPEC + Phage 3h" = expression(
                "UPEC 8923 +" ~ Phi ~ "MM02 SIM"
            ),
            "Cells + UPEC 24h" = "UPEC 8923",
            "Cells + UPEC + Phage 24/19h" = expression(
                "UPEC 8923 +" ~ Phi ~ "MM02 5h PI"
            )
        )
    ) +
    scale_y_continuous(
        transform = transform_pseudo_log(base = 10),
        breaks = c(0, 10^3, 10^6),
        labels = label_log(base = 10),
        expand = expansion(mult = c(0, 0.1))
    ) +
    expand_limits(y = 10^6) +
    theme_bw(base_size = 10) +
    theme(axis.title.x = element_blank(), legend.position = "none") +
    scale_fill_okabe_ito() +
    suppressWarnings(geom_signif(
        data = data.frame(
            group = as.character(1:2),
            start = c("Cells + UPEC 3h", "Cells + UPEC 24h"),
            end = c("Cells + UPEC + Phage 3h", "Cells + UPEC + Phage 24/19h"),
            label = c("×0", "×1.38"),
            y = c(5.25, 5.25)
        ),
        aes(xmin = start, xmax = end, annotations = label, y_position = y),
        manual = TRUE,
        inherit.aes = FALSE,
        size = 0.25,
        textsize = 10 * 0.8 / .pt
    ))) /
    (ggplot(data7, aes(name, value, fill = grepl("Phage", name))) +
        stat_summary(fun = "mean", geom = "bar", show.legend = FALSE) +
        stat_summary(
            fun.data = "mean_cl_normal",
            geom = "errorbar",
            linewidth = 0.25,
            show.legend = FALSE
        ) +
        geom_jitter(position = position_jitter(height = 0, seed = 5)) +
        facet_grid(
            cols = vars(group),
            scales = "free_x",
            labeller = as_labeller(
                c(
                    `1` = "'3/5 h (' * varphi ~ '3/5 h)'",
                    `2` = "'3/24 h (' * varphi ~ '5/24 h)'"
                ),
                default = label_parsed
            )
        ) +
        ylab("CFU/ml") +
        scale_x_discrete(
            labels = c(
                "Cells + UPEC 3h" = "UPEC 7958",
                "Cells + UPEC + Phage 3h" = expression(
                    "UPEC 7958 +" ~ Phi ~ "G10400 SIM"
                ),
                "Cells + UPEC 24 h" = "UPEC 7958",
                "Cells + UPEC + Phage 24/19 h" = expression(
                    "UPEC 7958  +" ~ Phi ~ "G10400 5h PI"
                )
            )
        ) +
        scale_y_continuous(
            transform = transform_pseudo_log(base = 10),
            breaks = c(0, 10^3, 10^6),
            labels = label_log(base = 10),
            expand = expansion(mult = c(0, 0.1))
        ) +
        expand_limits(y = 10^6) +
        theme_bw(base_size = 10) +
        theme(axis.title.x = element_blank(), legend.position = "none") +
        scale_fill_okabe_ito() +
        suppressWarnings(geom_signif(
            data = data.frame(
                group = as.character(1:2),
                start = c("Cells + UPEC 3h", "Cells + UPEC 24 h"),
                end = c(
                    "Cells + UPEC + Phage 3h",
                    "Cells + UPEC + Phage 24/19 h"
                ),
                label = c("×0", "×0.92"),
                y = c(5.25, 5.25)
            ),
            aes(xmin = start, xmax = end, annotations = label, y_position = y),
            manual = TRUE,
            inherit.aes = FALSE,
            size = 0.25,
            textsize = 10 * 0.8 / .pt
        ))) +
    plot_annotation(tag_levels = "A")

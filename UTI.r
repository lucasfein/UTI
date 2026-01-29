library(readxl)

library(dplyr)
library(tidyr)

library(forcats)

library(ggplot2)

library(scales)

library(ggsignif)
library(patchwork)

library(dunn.test)

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

dunn.test(
    data1 %>% pull(value),
    data1 %>% pull(name),
    method = "holm",
    list = TRUE,
    altp = TRUE
)
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
        fun.min = "min",
        fun.max = "max",
        geom = "errorbar",
        linewidth = 0.25,
        show.legend = FALSE
    ) +
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
        expand = expansion(mult = c(0, 0.1))
    ) +
    theme_bw(base_size = 9) +
    theme(axis.title.x = element_blank(), legend.position = "none") +
    geom_signif(
        comparisons = list(c("UPEC 8923 + MM02 24h", "UPEC 8923 24h")),
        annotations = "p = 0.01",
        size = 0.25,
        textsize = 9 * 0.8 / .pt
    )) /
    (ggplot(data2, aes(name, value, fill = name)) +
        stat_summary(fun = "mean", geom = "bar", show.legend = FALSE) +
        stat_summary(
            fun.min = "min",
            fun.max = "max",
            geom = "errorbar",
            linewidth = 0.25,
            show.legend = FALSE
        ) +
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
            expand = expansion(mult = c(0, 0.1))
        ) +
        theme_bw(base_size = 9) +
        theme(axis.title.x = element_blank(), legend.position = "none") +
        geom_signif(
            comparisons = list(c("UPEC 7958", "UPEC 7958 + G10400")),
            annotations = "p = 0.01",
            size = 0.25,
            textsize = 9 * 0.8 / .pt
        )) +
    plot_annotation(tag_levels = "A")

data3 <- bind_rows(
    read_excel(
        "data/UTI Model Bacterial and phage counts.16.1.26xlsx.xlsx",
        sheet = "MM02",
        range = "B33:N40",
        col_names = FALSE,
        .name_repair = "unique_quiet"
    ) %>%
        select(-(2:12)) %>%
        rename(name = `...1`, value = `...13`),
    read_excel(
        "data/UTI Model Bacterial and phage counts.16.1.26xlsx.xlsx",
        sheet = "MM02",
        range = "B50:N55",
        col_names = FALSE,
        .name_repair = "unique_quiet"
    ) %>%
        select(-(2:12)) %>%
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

data4 <- bind_rows(
    read_excel(
        "data/UTI Model Bacterial and phage counts.16.1.26xlsx.xlsx",
        sheet = "P00",
        range = "B23:N30",
        col_names = FALSE,
        .name_repair = "unique_quiet"
    ) %>%
        select(-(2:12)) %>%
        rename(name = `...1`, value = `...13`),
    read_excel(
        "data/UTI Model Bacterial and phage counts.16.1.26xlsx.xlsx",
        sheet = "P00",
        range = "B39:N44",
        col_names = FALSE,
        .name_repair = "unique_quiet"
    ) %>%
        select(-(2:12)) %>%
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

(ggplot(data3, aes(name, value, fill = grepl("cells only", name))) +
    stat_summary(fun = "mean", geom = "bar", show.legend = FALSE) +
    stat_summary(
        fun.min = "min",
        fun.max = "max",
        geom = "errorbar",
        linewidth = 0.25,
        show.legend = FALSE
    ) +
    expand_limits(y = 0) +
    ylab("PFU/ml") +
    scale_x_discrete(
        labels = c(
            "MM02 + UPEC 7958 3h" = expression(
                "UPEC 7958 +" ~ Phi ~ "MM02 3h"
            ),
            "MM02 cells only 3h" = expression(Phi ~ "MM02 3h"),
            "MM02 + UPEC 7958 24h" = expression(
                "UPEC 7958" ~ Phi ~ "MM02 24h"
            ),
            "MM02 cells only 24h" = expression(Phi ~ "MM02 24h")
        )
    ) +
    scale_y_continuous(
        transform = transform_pseudo_log(base = 10),
        breaks = c(0, 10^3, 10^6),
        labels = label_log(base = 10),
        expand = expansion(mult = c(0, 0.1))
    ) +
    theme_bw(base_size = 9) +
    theme(axis.title.x = element_blank(), legend.position = "none")) /
    (ggplot(data4, aes(name, value, fill = grepl("cells only", name))) +
        stat_summary(fun = "mean", geom = "bar", show.legend = FALSE) +
        stat_summary(
            fun.min = "min",
            fun.max = "max",
            geom = "errorbar",
            linewidth = 0.25,
            show.legend = FALSE
        ) +
        expand_limits(y = 0) +
        ylab("PFU/ml") +
        scale_x_discrete(
            labels = c(
                "G10400 + UPEC 8923 3h" = expression(
                    "UPEC 8923 +" ~ Phi ~ "G10400 3h"
                ),
                "G10400 cells only 3h" = expression(Phi ~ "G10400 3h"),
                "G10400 + UPEC 8923 24h" = expression(
                    "UPEC 8923" ~ Phi ~ "G10400 24h"
                ),
                "G10400 cells only 24h" = expression(Phi ~ "G10400 24h")
            )
        ) +
        scale_y_continuous(
            transform = transform_pseudo_log(base = 10),
            breaks = c(0, 10^3, 10^6),
            labels = label_log(base = 10),
            expand = expansion(mult = c(0, 0.1))
        ) +
        theme_bw(base_size = 9) +
        theme(axis.title.x = element_blank(), legend.position = "none")) +
    plot_annotation(tag_levels = "A")

data5 <- bind_rows(
    read_excel(
        "data/UTI Model Bacterial and phage counts.16.1.26xlsx.xlsx",
        sheet = "MM02",
        range = "B15:G22",
        col_names = FALSE,
        .name_repair = "unique_quiet"
    ) %>%
        select(-(2:5)) %>%
        rename(name = `...1`, value = `...6`),
    read_excel(
        "data/UTI Model Bacterial and phage counts.16.1.26xlsx.xlsx",
        sheet = "MM02",
        range = "B24:G29",
        col_names = FALSE,
        .name_repair = "unique_quiet"
    ) %>%
        select(-(2:5)) %>%
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

wilcox.test(
    data5 %>% filter(name == "UPEC 1h") %>% pull(value),
    data5 %>% filter(name == "MM02 after 1,5h") %>% pull(value)
)
wilcox.test(
    data5 %>% filter(name == "UPEC 0,5h") %>% pull(value),
    data5 %>% filter(name == "MM02 after 1h") %>% pull(value)
)

ggplot(data5, aes(name, value, fill = grepl("MM02", name))) +
    stat_summary(fun = "mean", geom = "bar", show.legend = FALSE) +
    stat_summary(
        fun.min = "min",
        fun.max = "max",
        geom = "errorbar",
        linewidth = 0.25,
        show.legend = FALSE
    ) +
    ylab("CFU/ml") +
    scale_x_discrete(
        labels = c(
            "UPEC 1h" = "UPEC 8923 1/24h",
            "MM02 after 1,5h" = expression(
                "UPEC 8923 1/24h +" ~ Phi ~ "MM02 1.5/24h"
            ),
            "UPEC 0,5h" = "UPEC 8923 0.5/24h",
            "MM02 after 1h" = expression(
                "UPEC 8923 0.5/24h +" ~ Phi ~ "MM02 1/24h"
            )
        )
    ) +
    scale_y_continuous(
        labels = label_number(),
        expand = expansion(mult = c(0, 0.1))
    ) +
    theme_bw(base_size = 9) +
    theme(axis.title.x = element_blank(), legend.position = "none")

data6 <- read_excel(
    "data/UTI Model Bacterial and phage counts.16.1.26xlsx.xlsx",
    sheet = "MM02",
    range = "B2:G13",
    col_names = FALSE,
    .name_repair = "unique_quiet"
) %>%
    select(-(2:5)) %>%
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

data7 <- read_excel(
    "data/UTI Model Bacterial and phage counts.16.1.26xlsx.xlsx",
    sheet = "P00",
    range = "B2:G13",
    col_names = FALSE,
    .name_repair = "unique_quiet"
) %>%
    select(-(2:5)) %>%
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

(ggplot(data6, aes(name, value, fill = grepl("Phage", name))) +
    stat_summary(fun = "mean", geom = "bar", show.legend = FALSE) +
    stat_summary(
        fun.min = "min",
        fun.max = "max",
        geom = "errorbar",
        linewidth = 0.25,
        show.legend = FALSE
    ) +
    ylab("CFU/ml") +
    scale_x_discrete(
        labels = c(
            "Cells + UPEC 3h" = "UPEC 8923 3/5h",
            "Cells + UPEC + Phage 3h" = expression(
                "UPEC 8923 +" ~ Phi ~ "MM02 3/5h"
            ),
            "Cells + UPEC 24h" = "UPEC 8923 3/24h",
            "Cells + UPEC + Phage 24/19h" = expression(
                "UPEC 8923 3/24 h +" ~ Phi ~ "MM02 5/24h"
            )
        )
    ) +
    scale_y_continuous(
        labels = label_number(),
        expand = expansion(mult = c(0, 0.1))
    ) +
    theme_bw(base_size = 9) +
    theme(axis.title.x = element_blank(), legend.position = "none")) /
    (ggplot(data7, aes(name, value, fill = grepl("Phage", name))) +
        stat_summary(fun = "mean", geom = "bar", show.legend = FALSE) +
        stat_summary(
            fun.min = "min",
            fun.max = "max",
            geom = "errorbar",
            linewidth = 0.25,
            show.legend = FALSE
        ) +
        ylab("CFU/ml") +
        scale_x_discrete(
            labels = c(
                "Cells + UPEC 3h" = "UPEC 7958 3/5h",
                "Cells + UPEC + Phage 3h" = expression(
                    "UPEC 7958 +" ~ Phi ~ "G10400 3/5h"
                ),
                "Cells + UPEC 24 h" = "UPEC 7958 3/24h",
                "Cells + UPEC + Phage 24/19 h" = expression(
                    "UPEC 7958 3/24 h +" ~ Phi ~ "G10400 5/24h"
                )
            )
        ) +
        scale_y_continuous(
            labels = label_number(),
            expand = expansion(mult = c(0, 0.1))
        ) +
        theme_bw(base_size = 9) +
        theme(axis.title.x = element_blank(), legend.position = "none")) +
    plot_annotation(tag_levels = "A")

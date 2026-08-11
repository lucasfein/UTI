library(readxl)

library(dplyr)
library(stringr)
library(tidyr)

library(forcats)

library(ggplot2)

library(scales)

library(ggsignif)
library(patchwork)

library(dunn.test)

library(ggokabeito)

library(ragg)
library(ggtext)

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

ggsave(
    "Figure-6.png",
    (ggplot(data1, aes(name, value, fill = name)) +
        stat_summary(fun = "mean", geom = "bar", show.legend = FALSE) +
        stat_summary(
            fun.min = "min",
            fun.max = "max",
            geom = "errorbar",
            linewidth = 0.25,
            show.legend = FALSE
        ) +
        geom_jitter(position = position_jitter(height = 0, seed = 5)) +
        ylab("LDH (U/L)") +
        scale_x_discrete(
            labels = c(
                "UPEC 8923 24h" = "UPEC 8923 24h",
                "UPEC 8923 + MM02 24h" = "UPEC 8923 + &Phi; MM02 24h",
                "MM02 24h" = "&Phi; MM02 24h",
                "Cells only 24h" = "Cells only 24h"
            )
        ) +
        scale_y_continuous(
            labels = label_number(),
            expand = expansion(mult = c(0, 0.1)),
            breaks = seq(0, 1000, 250)
        ) +
        expand_limits(y = 1000) +
        theme_bw(base_size = 10) +
        theme(
            text = element_text(family = "Arial"),
            axis.title.x = element_blank(),
            legend.position = "none",
            axis.text.x = element_markdown()
        ) +
        geom_signif(
            comparisons = list(c("UPEC 8923 + MM02 24h", "UPEC 8923 24h")),
            annotations = "*",
            size = 0.25,
            textsize = 10 * 0.8 / .pt,
            family = "Arial",
            y_position = 1000,
            tip_length = 0.025
        ) +
        scale_fill_okabe_ito()) /
        (ggplot(data2, aes(name, value, fill = name)) +
            stat_summary(fun = "mean", geom = "bar", show.legend = FALSE) +
            stat_summary(
                fun.min = "min",
                fun.max = "max",
                geom = "errorbar",
                linewidth = 0.25,
                show.legend = FALSE
            ) +
            geom_jitter(position = position_jitter(height = 0, seed = 5)) +
            ylab("LDH (U/L)") +
            scale_x_discrete(
                labels = c(
                    "UPEC 7958" = "UPEC 7958 24h",
                    "UPEC 7958 + G10400" = "UPEC 7958 + &Phi; G10400 24h",
                    "G10400 24h" = "&Phi; G10400 24h",
                    "Cells only 24h" = "Cells only 24h"
                )
            ) +
            scale_y_continuous(
                labels = label_number(),
                expand = expansion(mult = c(0, 0.1)),
                breaks = seq(0, 1000, 250)
            ) +
            expand_limits(y = 1000) +
            theme_bw(base_size = 10) +
            theme(
                text = element_text(family = "Arial"),
                axis.title.x = element_blank(),
                legend.position = "none",
                axis.text.x = element_markdown()
            ) +
            geom_signif(
                comparisons = list(c("UPEC 7958", "UPEC 7958 + G10400")),
                annotations = "*",
                size = 0.25,
                textsize = 10 * 0.8 / .pt,
                family = "Arial",
                y_position = 1000,
                tip_length = 0.025
            ) +
            scale_fill_okabe_ito()) +
        plot_annotation(tag_levels = "A"),
    width = 7,
    height = 7,
    units = "in",
    dpi = 1200,
    device = ragg::agg_png()
)

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
    log10(value) ~ fct_rev(name),
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
    log10(value) ~ fct_rev(name),
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
    log10(value) ~ fct_rev(name),
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
    log10(value) ~ fct_rev(name),
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

ggsave(
    "Figure-1.png",
    (ggplot(
        data3 %>% filter(group == 1),
        aes(name, value, fill = grepl("cells only", name))
    ) +
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
            labeller = as_labeller(c(`1` = "3 h"))
        ) +
        expand_limits(y = 0) +
        ylab("PFU/ml") +
        scale_x_discrete(
            labels = c(
                "MM02 + UPEC 7958 3h" = "UPEC 7958 + &Phi; MM02",
                "MM02 cells only 3h" = "&Phi; MM02",
                "MM02 + UPEC 7958 24h" = "UPEC 7958 + &Phi; MM02",
                "MM02 cells only 24h" = "&Phi; MM02"
            )
        ) +
        scale_y_continuous(
            transform = transform_pseudo_log(base = 10),
            breaks = c(0, 10^3, 10^6),
            labels = label_log(base = 10),
            expand = expansion(mult = c(0, 0.1))
        ) +
        expand_limits(y = 10^6.5) +
        theme_bw(base_size = 10) +
        theme(
            text = element_text(family = "Arial"),
            axis.title.x = element_blank(),
            legend.position = "none",
            axis.text.x = element_markdown()
        ) +
        suppressWarnings(geom_signif(
            data = data.frame(
                group = as.character(1),
                start = "MM02 + UPEC 7958 3h",
                end = "MM02 cells only 3h",
                label = "**",
                y = 6.75,
                tip_length = 0.025
            ),
            aes(xmin = start, xmax = end, annotations = label, y_position = y),
            manual = TRUE,
            inherit.aes = FALSE,
            size = 0.25,
            textsize = 10 * 0.8 / .pt,
            family = "Arial"
        )) +
        scale_fill_okabe_ito() |
        ggplot(
            data4 %>% filter(group == 1),
            aes(name, value, fill = grepl("cells only", name))
        ) +
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
                labeller = as_labeller(c(`1` = "3 h"))
            ) +
            ylab("PFU/ml") +
            scale_x_discrete(
                labels = c(
                    "G10400 + UPEC 8923 3h" = "UPEC 8923 + &Phi; G10400",
                    "G10400 cells only 3h" = "&Phi; G10400",
                    "G10400 + UPEC 8923 24h" = "UPEC 8923 + &Phi; G10400",
                    "G10400 cells only 24h" = "&Phi; G10400"
                )
            ) +
            scale_y_continuous(
                transform = transform_pseudo_log(base = 10),
                breaks = c(0, 10^3, 10^6),
                labels = label_log(base = 10),
                expand = expansion(mult = c(0, 0.1))
            ) +
            expand_limits(y = 10^6.5) +
            theme_bw(base_size = 10) +
            theme(
                text = element_text(family = "Arial"),
                axis.title.x = element_blank(),
                legend.position = "none",
                axis.text.x = element_markdown()
            ) +
            suppressWarnings(geom_signif(
                data = data.frame(
                    group = as.character(1),
                    start = "G10400 + UPEC 8923 3h",
                    end = "G10400 cells only 3h",
                    label = "***",
                    y = 6.75,
                    tip_length = 0.025
                ),
                aes(
                    xmin = start,
                    xmax = end,
                    annotations = label,
                    y_position = y
                ),
                manual = TRUE,
                inherit.aes = FALSE,
                size = 0.25,
                textsize = 10 * 0.8 / .pt,
                family = "Arial"
            )) +
            scale_fill_okabe_ito() +
            plot_layout(tag_level = 'new')) /
        (ggplot(
            data3 %>% filter(group == 2),
            aes(name, value, fill = grepl("cells only", name))
        ) +
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
                labeller = as_labeller(c(`2` = "24 h"))
            ) +
            expand_limits(y = 0) +
            ylab("PFU/ml") +
            scale_x_discrete(
                labels = c(
                    "MM02 + UPEC 7958 3h" = "UPEC 7958 + &Phi; MM02",
                    "MM02 cells only 3h" = "&Phi; MM02",
                    "MM02 + UPEC 7958 24h" = "UPEC 7958 + &Phi; MM02",
                    "MM02 cells only 24h" = "&Phi; MM02"
                )
            ) +
            scale_y_continuous(
                transform = transform_pseudo_log(base = 10),
                breaks = c(0, 10^3, 10^6),
                labels = label_log(base = 10),
                expand = expansion(mult = c(0, 0.1))
            ) +
            expand_limits(y = 10^6.5) +
            theme_bw(base_size = 10) +
            theme(
                text = element_text(family = "Arial"),
                axis.title.x = element_blank(),
                legend.position = "none",
                axis.text.x = element_markdown()
            ) +
            suppressWarnings(geom_signif(
                data = data.frame(
                    group = as.character(2),
                    start = "MM02 + UPEC 7958 24h",
                    end = "MM02 cells only 24h",
                    label = "ns",
                    y = 6.75,
                    tip_length = 0.025
                ),
                aes(
                    xmin = start,
                    xmax = end,
                    annotations = label,
                    y_position = y
                ),
                manual = TRUE,
                inherit.aes = FALSE,
                size = 0.25,
                textsize = 10 * 0.8 / .pt,
                family = "Arial"
            )) +
            scale_fill_okabe_ito() |
            ggplot(
                data4 %>% filter(group == 2),
                aes(name, value, fill = grepl("cells only", name))
            ) +
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
                    labeller = as_labeller(c(`2` = "24 h"))
                ) +
                ylab("PFU/ml") +
                scale_x_discrete(
                    labels = c(
                        "G10400 + UPEC 8923 3h" = "UPEC 8923 + &Phi; G10400",
                        "G10400 cells only 3h" = "&Phi; G10400",
                        "G10400 + UPEC 8923 24h" = "UPEC 8923 + &Phi; G10400",
                        "G10400 cells only 24h" = "&Phi; G10400"
                    )
                ) +
                scale_y_continuous(
                    transform = transform_pseudo_log(base = 10),
                    breaks = c(0, 10^3, 10^6),
                    labels = label_log(base = 10),
                    expand = expansion(mult = c(0, 0.1))
                ) +
                expand_limits(y = 10^6.5) +
                theme_bw(base_size = 10) +
                theme(
                    text = element_text(family = "Arial"),
                    axis.title.x = element_blank(),
                    legend.position = "none",
                    axis.text.x = element_markdown()
                ) +
                scale_fill_okabe_ito() +
                suppressWarnings(geom_signif(
                    data = data.frame(
                        group = as.character(2),
                        start = "G10400 + UPEC 8923 24h",
                        end = "G10400 cells only 24h",
                        label = "***",
                        y = 6.75,
                        tip_length = 0.025
                    ),
                    aes(
                        xmin = start,
                        xmax = end,
                        annotations = label,
                        y_position = y
                    ),
                    manual = TRUE,
                    inherit.aes = FALSE,
                    size = 0.25,
                    textsize = 10 * 0.8 / .pt,
                    family = "Arial"
                )) +
                plot_layout(tag_level = 'new')) +
        plot_annotation(tag_levels = "A"),
    width = 7,
    height = 7,
    units = "in",
    dpi = 1200,
    device = ragg::agg_png()
)

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
    log10(value) ~ fct_rev(name),
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
    log10(value) ~ fct_rev(name),
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

ggsave(
    "Figure-4.png",
    (ggplot(
        data5 %>% filter(group == 1),
        aes(name, value, fill = grepl("MM02", name))
    ) +
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
                    `1` = "1/24 h (&phi; 1.5/24 h)"
                )
            )
        ) +
        ylab("CFU/ml") +
        scale_x_discrete(
            labels = c(
                "UPEC 1h" = "UPEC 8923",
                "MM02 after 1,5h" = "UPEC 8923 + &Phi; MM02 1.5h PI",
                "UPEC 0,5h" = "UPEC 8923",
                "MM02 after 1h" = "UPEC 8923 + &Phi; MM02 1h PI"
            )
        ) +
        scale_y_continuous(
            transform = transform_pseudo_log(base = 10),
            breaks = c(0, 10^3, 10^6),
            labels = label_log(base = 10),
            expand = expansion(mult = c(0, 0.1))
        ) +
        expand_limits(y = 10^6.5) +
        theme_bw(base_size = 10) +
        theme(
            text = element_text(family = "Arial"),
            axis.title.x = element_blank(),
            legend.position = "none",
            strip.text = element_markdown(),
            axis.text.x = element_markdown()
        ) +
        suppressWarnings(geom_signif(
            data = data.frame(
                group = as.character(1),
                start = "UPEC 1h",
                end = "MM02 after 1,5h",
                label = "×0.74",
                y = 6.75
            ),
            aes(xmin = start, xmax = end, annotations = label, y_position = y),
            manual = TRUE,
            inherit.aes = FALSE,
            size = 0.25,
            textsize = 10 * 0.8 / .pt,
            family = "Arial"
        )) +
        scale_fill_okabe_ito() |
        ggplot(
            data5 %>% filter(group == 2),
            aes(name, value, fill = grepl("MM02", name))
        ) +
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
                        `2` = "0.5/24 h (&phi; 1/24 h)"
                    )
                )
            ) +
            ylab("CFU/ml") +
            scale_x_discrete(
                labels = c(
                    "UPEC 1h" = "UPEC 8923",
                    "MM02 after 1,5h" = "UPEC 8923 + &Phi; MM02 1.5h PI",
                    "UPEC 0,5h" = "UPEC 8923",
                    "MM02 after 1h" = "UPEC 8923 + &Phi; MM02 1h PI"
                )
            ) +
            scale_y_continuous(
                transform = transform_pseudo_log(base = 10),
                breaks = c(0, 10^3, 10^6),
                labels = label_log(base = 10),
                expand = expansion(mult = c(0, 0.1))
            ) +
            expand_limits(y = 10^6.5) +
            theme_bw(base_size = 10) +
            theme(
                text = element_text(family = "Arial"),
                axis.title.x = element_blank(),
                legend.position = "none",
                strip.text = element_markdown(),
                axis.text.x = element_markdown()
            ) +
            suppressWarnings(geom_signif(
                data = data.frame(
                    group = as.character(2),
                    start = "UPEC 0,5h",
                    end = "MM02 after 1h",
                    label = "×0.64",
                    y = 6.75
                ),
                aes(
                    xmin = start,
                    xmax = end,
                    annotations = label,
                    y_position = y
                ),
                manual = TRUE,
                inherit.aes = FALSE,
                size = 0.25,
                textsize = 10 * 0.8 / .pt,
                family = "Arial"
            )) +
            scale_fill_okabe_ito()) +
        plot_annotation(tag_levels = "A"),
    width = 7,
    height = 3.5,
    units = "in",
    dpi = 1200,
    device = ragg::agg_png()
)

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
    log10(value) ~ fct_rev(name),
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
    log10(value) ~ fct_rev(name),
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

ggsave(
    "Figure-3.png",
    (ggplot(
        data6 %>% filter(group == 1),
        aes(name, value, fill = grepl("Phage", name))
    ) +
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
                    `1` = "3/5 h (&phi; 3/5 h)"
                )
            )
        ) +
        ylab("CFU/ml") +
        scale_x_discrete(
            labels = c(
                "Cells + UPEC 3h" = "UPEC 8923",
                "Cells + UPEC + Phage 3h" = "UPEC 8923 + &Phi; MM02 SIM",
                "Cells + UPEC 24h" = "UPEC 8923",
                "Cells + UPEC + Phage 24/19h" = "UPEC 8923 + &Phi; MM02 5h PI"
            )
        ) +
        scale_y_continuous(
            transform = transform_pseudo_log(base = 10),
            breaks = c(0, 10^3, 10^6),
            labels = label_log(base = 10),
            expand = expansion(mult = c(0, 0.1))
        ) +
        expand_limits(y = 10^6.5) +
        theme_bw(base_size = 10) +
        theme(
            text = element_text(family = "Arial"),
            axis.title.x = element_blank(),
            legend.position = "none",
            strip.text = element_markdown(),
            axis.text.x = element_markdown()
        ) +
        scale_fill_okabe_ito() +
        suppressWarnings(geom_signif(
            data = data.frame(
                group = as.character(1),
                start = "Cells + UPEC 3h",
                end = "Cells + UPEC + Phage 3h",
                label = "×0",
                y = 6.75,
                tip_length = 0.025
            ),
            aes(
                xmin = start,
                xmax = end,
                annotations = label,
                y_position = y
            ),
            manual = TRUE,
            inherit.aes = FALSE,
            size = 0.25,
            textsize = 10 * 0.8 / .pt,
            family = "Arial"
        )) |
        ggplot(
            data7 %>% filter(group == 1),
            aes(name, value, fill = grepl("Phage", name))
        ) +
            stat_summary(
                fun = "mean",
                geom = "bar",
                show.legend = FALSE
            ) +
            stat_summary(
                fun.data = "mean_cl_normal",
                geom = "errorbar",
                linewidth = 0.25,
                show.legend = FALSE
            ) +
            geom_jitter(
                position = position_jitter(height = 0, seed = 5)
            ) +
            facet_grid(
                cols = vars(group),
                scales = "free_x",
                labeller = as_labeller(
                    c(
                        `1` = "3/5 h (&phi; 3/5 h)"
                    )
                )
            ) +
            ylab("CFU/ml") +
            scale_x_discrete(
                labels = c(
                    "Cells + UPEC 3h" = "UPEC 7958",
                    "Cells + UPEC + Phage 3h" = "UPEC 7958 + &Phi; G10400 SIM",
                    "Cells + UPEC 24 h" = "UPEC 7958",
                    "Cells + UPEC + Phage 24/19 h" = "UPEC 7958 + &Phi; G10400 5h PI"
                )
            ) +
            scale_y_continuous(
                transform = transform_pseudo_log(base = 10),
                breaks = c(0, 10^3, 10^6),
                labels = label_log(base = 10),
                expand = expansion(mult = c(0, 0.1))
            ) +
            expand_limits(y = 10^6.5) +
            theme_bw(base_size = 10) +
            theme(
                text = element_text(family = "Arial"),
                axis.title.x = element_blank(),
                legend.position = "none",
                strip.text = element_markdown(),
                axis.text.x = element_markdown()
            ) +
            scale_fill_okabe_ito() +
            suppressWarnings(geom_signif(
                data = data.frame(
                    group = as.character(1),
                    start = "Cells + UPEC 3h",
                    end = "Cells + UPEC + Phage 3h",
                    label = "×0",
                    y = 6.75,
                    tip_length = 0.025
                ),
                aes(
                    xmin = start,
                    xmax = end,
                    annotations = label,
                    y_position = y
                ),
                manual = TRUE,
                inherit.aes = FALSE,
                size = 0.25,
                textsize = 10 * 0.8 / .pt,
                family = "Arial"
            )) +
            plot_layout(tag_level = 'new')) /
        (ggplot(
            data6 %>% filter(group == 2),
            aes(name, value, fill = grepl("Phage", name))
        ) +
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
                        `2` = "3/24 h (&phi; 5/24 h)"
                    )
                )
            ) +
            ylab("CFU/ml") +
            scale_x_discrete(
                labels = c(
                    "Cells + UPEC 3h" = "UPEC 8923",
                    "Cells + UPEC + Phage 3h" = "UPEC 8923 + &Phi; MM02 SIM",
                    "Cells + UPEC 24h" = "UPEC 8923",
                    "Cells + UPEC + Phage 24/19h" = "UPEC 8923 + &Phi; MM02 5h PI"
                )
            ) +
            scale_y_continuous(
                transform = transform_pseudo_log(base = 10),
                breaks = c(0, 10^3, 10^6),
                labels = label_log(base = 10),
                expand = expansion(mult = c(0, 0.1))
            ) +
            expand_limits(y = 10^6.5) +
            theme_bw(base_size = 10) +
            theme(
                text = element_text(family = "Arial"),
                axis.title.x = element_blank(),
                legend.position = "none",
                strip.text = element_markdown(),
                axis.text.x = element_markdown()
            ) +
            scale_fill_okabe_ito() +
            suppressWarnings(geom_signif(
                data = data.frame(
                    group = as.character(2),
                    start = "Cells + UPEC 24h",
                    end = "Cells + UPEC + Phage 24/19h",
                    label = "×1.38",
                    y = 6.75,
                    tip_length = 0.025
                ),
                aes(
                    xmin = start,
                    xmax = end,
                    annotations = label,
                    y_position = y
                ),
                manual = TRUE,
                inherit.aes = FALSE,
                size = 0.25,
                textsize = 10 * 0.8 / .pt,
                family = "Arial"
            )) |
            ggplot(
                data7 %>% filter(group == 2),
                aes(name, value, fill = grepl("Phage", name))
            ) +
                stat_summary(
                    fun = "mean",
                    geom = "bar",
                    show.legend = FALSE
                ) +
                stat_summary(
                    fun.data = "mean_cl_normal",
                    geom = "errorbar",
                    linewidth = 0.25,
                    show.legend = FALSE
                ) +
                geom_jitter(
                    position = position_jitter(height = 0, seed = 5)
                ) +
                facet_grid(
                    cols = vars(group),
                    scales = "free_x",
                    labeller = as_labeller(
                        c(
                            `2` = "3/24 h (&phi; 5/24 h)"
                        )
                    )
                ) +
                ylab("CFU/ml") +
                scale_x_discrete(
                    labels = c(
                        "Cells + UPEC 3h" = "UPEC 7958",
                        "Cells + UPEC + Phage 3h" = "UPEC 7958 + &Phi; G10400 SIM",
                        "Cells + UPEC 24 h" = "UPEC 7958",
                        "Cells + UPEC + Phage 24/19 h" = "UPEC 7958 + &Phi; G10400 5h PI"
                    )
                ) +
                scale_y_continuous(
                    transform = transform_pseudo_log(base = 10),
                    breaks = c(0, 10^3, 10^6),
                    labels = label_log(base = 10),
                    expand = expansion(mult = c(0, 0.1))
                ) +
                expand_limits(y = 10^6.5) +
                theme_bw(base_size = 10) +
                theme(
                    text = element_text(family = "Arial"),
                    axis.title.x = element_blank(),
                    legend.position = "none",
                    strip.text = element_markdown(),
                    axis.text.x = element_markdown()
                ) +
                scale_fill_okabe_ito() +
                suppressWarnings(geom_signif(
                    data = data.frame(
                        group = as.character(2),
                        start = "Cells + UPEC 24 h",
                        end = "Cells + UPEC + Phage 24/19 h",
                        label = "×0.92",
                        y = 6.75,
                        tip_length = 0.025
                    ),
                    aes(
                        xmin = start,
                        xmax = end,
                        annotations = label,
                        y_position = y
                    ),
                    manual = TRUE,
                    inherit.aes = FALSE,
                    size = 0.25,
                    textsize = 10 * 0.8 / .pt,
                    family = "Arial"
                )) +
                plot_layout(tag_level = 'new')) +
        plot_annotation(tag_levels = "A"),
    width = 7,
    height = 7,
    units = "in",
    dpi = 1200,
    device = ragg::agg_png()
)

data8 <- bind_rows(
    read_excel(
        "data/Lysosome conditions.xlsx",
        sheet = "UPEC 58 + P00",
        range = "B3:F5",
        col_names = FALSE
    ) %>%
        select(1, 5) %>%
        mutate(condition = "1"),
    read_excel(
        "data/Lysosome conditions.xlsx",
        sheet = "UPEC 58 + P00",
        range = "B9:F11",
        col_names = FALSE
    ) %>%
        select(1, 5) %>%
        mutate(condition = "1"),
    read_excel(
        "data/Lysosome conditions.xlsx",
        sheet = "UPEC 58 + P00",
        range = "B13:F15",
        col_names = FALSE
    ) %>%
        select(1, 5) %>%
        mutate(condition = "1"),
    read_excel(
        "data/Lysosome conditions.xlsx",
        sheet = "UPEC 58 + P00",
        range = "B18:F20",
        col_names = FALSE
    ) %>%
        select(1, 5) %>%
        mutate(condition = "2"),
    read_excel(
        "data/Lysosome conditions.xlsx",
        sheet = "UPEC 58 + P00",
        range = "B24:F26",
        col_names = FALSE
    ) %>%
        select(1, 5) %>%
        mutate(condition = "2"),
    read_excel(
        "data/Lysosome conditions.xlsx",
        sheet = "UPEC 58 + P00",
        range = "B28:F30",
        col_names = FALSE
    ) %>%
        select(1, 5) %>%
        mutate(condition = "2")
) %>%
    rename(name = `...1`, value = `...2`) %>%
    mutate(
        name = fct_relevel(
            case_when(
                str_detect(
                    str_squish(name),
                    "PC [1-3] UPEC 7958 pH(7| 4,5)"
                ) ~ "PC UPEC 7958",
                str_detect(
                    str_squish(name),
                    "UPEC 58 \\+ P00 MOI 1 BR[1-3]"
                ) ~ "UPEC 7958 + Phage G10400 MOI 1",
                str_detect(
                    str_squish(name),
                    "UPEC 58 \\+ P00 MOI 0[\\.,]01 BR[1-3]"
                ) ~ "UPEC 7958 + Phage G10400 MOI 0.01"
            ),
            "PC UPEC 7958",
            "UPEC 7958 + Phage G10400 MOI 1",
            "UPEC 7958 + Phage G10400 MOI 0.01"
        ),
        condition = fct_relevel(condition, "1", "2")
    )

data9 <- bind_rows(
    read_excel(
        "data/Lysosome conditions.xlsx",
        sheet = "UPEC 23 + MM02",
        range = "B3:F5",
        col_names = FALSE
    ) %>%
        select(1, 5) %>%
        mutate(condition = "1"),
    read_excel(
        "data/Lysosome conditions.xlsx",
        sheet = "UPEC 23 + MM02",
        range = "B9:F11",
        col_names = FALSE
    ) %>%
        select(1, 5) %>%
        mutate(condition = "1"),
    read_excel(
        "data/Lysosome conditions.xlsx",
        sheet = "UPEC 23 + MM02",
        range = "B13:F15",
        col_names = FALSE
    ) %>%
        select(1, 5) %>%
        mutate(condition = "1"),
    read_excel(
        "data/Lysosome conditions.xlsx",
        sheet = "UPEC 23 + MM02",
        range = "B18:F20",
        col_names = FALSE
    ) %>%
        select(1, 5) %>%
        mutate(condition = "2"),
    read_excel(
        "data/Lysosome conditions.xlsx",
        sheet = "UPEC 23 + MM02",
        range = "B24:F26",
        col_names = FALSE
    ) %>%
        select(1, 5) %>%
        mutate(condition = "2"),
    read_excel(
        "data/Lysosome conditions.xlsx",
        sheet = "UPEC 23 + MM02",
        range = "B28:F30",
        col_names = FALSE
    ) %>%
        select(1, 5) %>%
        mutate(condition = "2")
) %>%
    rename(name = `...1`, value = `...2`) %>%
    mutate(
        name = fct_relevel(
            case_when(
                str_detect(
                    str_squish(name),
                    "PC [1-3]"
                ) ~ "PC UPEC 8923",
                str_detect(
                    str_squish(name),
                    "UPEC 23 \\+ MM02 MOI 1 BR[1-3]"
                ) ~ "UPEC 8923 + Phage MM02 MOI 1",
                str_detect(
                    str_squish(name),
                    "UPEC 23 \\+ MM02 MOI 0[\\.,]01 BR[1-3]"
                ) ~ "UPEC 8923 + Phage MM02 MOI 0.01"
            ),
            "PC UPEC 8923",
            "UPEC 8923 + Phage MM02 MOI 1",
            "UPEC 8923 + Phage MM02 MOI 0.01"
        ),
        condition = fct_relevel(condition, "1", "2")
    )

p.adjust(
    c(
        t.test(
            I(log10(value + 1)) ~ name,
            data = data8 %>%
                filter(
                    name %in%
                        c("PC UPEC 7958", "UPEC 7958 + Phage G10400 MOI 1") &
                        condition == "1"
                )
        )$p.value,
        t.test(
            I(log10(value + 1)) ~ name,
            data = data8 %>%
                filter(
                    name %in%
                        c("PC UPEC 7958", "UPEC 7958 + Phage G10400 MOI 0.01") &
                        condition == "1"
                )
        )$p.value,
        t.test(
            I(log10(value + 1)) ~ name,
            data = data8 %>%
                filter(
                    name %in%
                        c("PC UPEC 7958", "UPEC 7958 + Phage G10400 MOI 1") &
                        condition == "2"
                )
        )$p.value,
        t.test(
            I(log10(value + 1)) ~ name,
            data = data8 %>%
                filter(
                    name %in%
                        c("PC UPEC 7958", "UPEC 7958 + Phage G10400 MOI 0.01") &
                        condition == "2"
                )
        )$p.value
    ),
    method = "holm"
)

t.test(
    I(log10(value + 1)) ~ name,
    data = data8 %>%
        filter(
            name %in%
                c("PC UPEC 7958", "UPEC 7958 + Phage G10400 MOI 0.01") &
                condition == "1"
        )
)

p.adjust(
    c(
        t.test(
            I(log10(value + 1)) ~ name,
            data = data9 %>%
                filter(
                    name %in%
                        c("PC UPEC 8923", "UPEC 8923 + Phage MM02 MOI 1") &
                        condition == "1"
                )
        )$p.value,
        t.test(
            I(log10(value + 1)) ~ name,
            data = data9 %>%
                filter(
                    name %in%
                        c("PC UPEC 8923", "UPEC 8923 + Phage MM02 MOI 0.01") &
                        condition == "1"
                )
        )$p.value,
        t.test(
            I(log10(value + 1)) ~ name,
            data = data9 %>%
                filter(
                    name %in%
                        c("PC UPEC 8923", "UPEC 8923 + Phage MM02 MOI 1") &
                        condition == "2"
                )
        )$p.value,
        t.test(
            I(log10(value + 1)) ~ name,
            data = data9 %>%
                filter(
                    name %in%
                        c("PC UPEC 8923", "UPEC 8923 + Phage MM02 MOI 0.01") &
                        condition == "2"
                )
        )$p.value
    ),
    method = "holm"
)

t.test(
    I(log10(value + 1)) ~ name,
    data = data9 %>%
        filter(
            name %in%
                c("PC UPEC 8923", "UPEC 8923 + Phage MM02 MOI 1") &
                condition == "1"
        )
)

t.test(
    I(log10(value + 1)) ~ name,
    data = data9 %>%
        filter(
            name %in%
                c("PC UPEC 8923", "UPEC 8923 + Phage MM02 MOI 0.01") &
                condition == "1"
        )
)

ggsave(
    "Figure-5.png",
    (ggplot(data8, aes(name, value, fill = name)) +
        stat_summary(fun = "mean", geom = "bar", show.legend = FALSE) +
        stat_summary(
            fun.data = "mean_cl_normal",
            geom = "errorbar",
            linewidth = 0.25,
            show.legend = FALSE
        ) +
        geom_jitter(
            position = position_jitter(height = 0, seed = 5),
            show.legend = FALSE
        ) +
        facet_grid(
            cols = vars(condition),
            labeller = as_labeller(c(`1` = "pH 7", `2` = "pH 4.5"))
        ) +
        scale_x_discrete(
            name = NULL,
            labels = c(
                "UPEC 7958 + Phage G10400 MOI 1" = "UPEC 7958 +<br>&Phi; G10400 MOI 1",
                "UPEC 7958 + Phage G10400 MOI 0.01" = "UPEC 7958 +<br>&Phi; G10400 MOI 0.01"
            )
        ) +
        scale_y_continuous(
            name = "CFU/ml",
            transform = transform_pseudo_log(),
            labels = label_log(),
            breaks = c(0, 10^3, 10^6, 10^9),
            expand = expansion(mult = c(0, 0.05))
        ) +
        scale_fill_okabe_ito() +
        theme_bw(base_size = 10) +
        theme(
            text = element_text(family = "Arial"),
            axis.text.x = element_markdown()
        ) +
        coord_cartesian(ylim = c(0, NA)) +
        suppressWarnings(geom_signif(
            data = data.frame(
                condition = "1",
                start = "PC UPEC 7958",
                end = "UPEC 7958 + Phage G10400 MOI 0.01",
                label = "**",
                y = 25
            ),
            aes(xmin = start, xmax = end, annotations = label, y_position = y),
            manual = TRUE,
            inherit.aes = FALSE,
            size = 0.25,
            textsize = 10 * 0.8 / .pt,
            family = "Arial"
        ))) /
        (ggplot(data9, aes(name, value, fill = name)) +
            stat_summary(fun = "mean", geom = "bar", show.legend = FALSE) +
            stat_summary(
                fun.data = "mean_cl_normal",
                geom = "errorbar",
                linewidth = 0.25,
                show.legend = FALSE
            ) +
            geom_jitter(
                position = position_jitter(height = 0, seed = 5),
                show.legend = FALSE
            ) +
            facet_grid(
                cols = vars(condition),
                labeller = as_labeller(c(`1` = "pH 7", `2` = "pH 4.5"))
            ) +
            scale_x_discrete(
                name = NULL,
                labels = c(
                    "UPEC 8923 + Phage MM02 MOI 1" = "UPEC 8923 +<br>&Phi; MM02 MOI 1",
                    "UPEC 8923 + Phage MM02 MOI 0.01" = "UPEC 8923 +<br>&Phi; MM02 MOI 0.01"
                )
            ) +
            scale_y_continuous(
                name = "CFU/ml",
                transform = transform_pseudo_log(),
                labels = label_log(),
                breaks = c(0, 10^3, 10^6, 10^9),
                expand = expansion(mult = c(0, 0.05))
            ) +
            scale_fill_okabe_ito() +
            theme_bw(base_size = 10) +
            theme(
                text = element_text(family = "Arial"),
                axis.text.x = element_markdown()
            ) +
            coord_cartesian(ylim = c(0, NA)) +
            suppressWarnings(geom_signif(
                data = data.frame(
                    condition = "1",
                    start = c("PC UPEC 8923", "PC UPEC 8923"),
                    end = c(
                        "UPEC 8923 + Phage MM02 MOI 1",
                        "UPEC 8923 + Phage MM02 MOI 0.01"
                    ),
                    label = c("*", "***"),
                    y = c(22.5, 25)
                ),
                aes(
                    xmin = start,
                    xmax = end,
                    annotations = label,
                    y_position = y
                ),
                manual = TRUE,
                inherit.aes = FALSE,
                size = 0.25,
                textsize = 10 * 0.8 / .pt,
                family = "Arial"
            ))) +
        plot_annotation(tag_levels = c("A", "B")),
    device = ragg::agg_png(),
    dpi = 1200,
    width = 7,
    height = 7,
    units = "in"
)

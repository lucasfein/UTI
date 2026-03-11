library(readxl)
library(dplyr)
library(tidyr)

library(forcats)

library(ggplot2)

library(ggokabeito)

library(scales)

library(ggsignif)

library(ragg)

library(ggtext)

data <- read_excel("data/Virucide excel.xlsx", range = "A1:D13")

data <- data %>%
    rename(phage = Phage) %>%
    mutate(
        phage = replace_values(
            phage,
            "G2062" ~ "G9062",
            "Ac-L3" ~ "AC-L3",
            "Ac-L3 + Virucide" ~ "AC-L3 + Virucide"
        )
    )

data <- data %>%
    mutate(
        group = replace_values(
            phage,
            "AC-L3" ~ "phage",
            "AC-L3 + Virucide" ~ "phage + virucide",
            "G10400" ~ "phage",
            "G10400 + Virucide" ~ "phage + virucide",
            "G2494" ~ "phage",
            "G2494 + Virucide" ~ "phage + virucide",
            "G9062" ~ "phage",
            "G9062 + Virucide" ~ "phage + virucide",
            "MM02" ~ "phage",
            "MM02 + Virucide" ~ "phage + virucide",
            "WFH" ~ "phage",
            "WFH + Virucide" ~ "phage + virucide"
        )
    ) %>%
    mutate(
        phage = replace_values(
            phage,
            "AC-L3" ~ "AC-L3",
            "AC-L3 + Virucide" ~ "AC-L3",
            "G10400" ~ "G10400",
            "G10400 + Virucide" ~ "G10400",
            "G2494" ~ "G2494",
            "G2494 + Virucide" ~ "G2494",
            "G9062" ~ "G9062",
            "G9062 + Virucide" ~ "G9062",
            "MM02" ~ "MM02",
            "MM02 + Virucide" ~ "MM02",
            "WFH" ~ "WFH",
            "WFH + Virucide" ~ "WFH"
        )
    )

data <- data %>%
    pivot_longer(V1:V3, names_to = "replicate", values_to = "value") %>%
    mutate(replicate = recode_values(replicate, "V1" ~ 1, "V2" ~ 2, "V3" ~ 3))

data <- data %>%
    mutate(
        phage = fct_inorder(phage),
        group = fct_relevel(group, "phage", "phage + virucide")
    )

p <- p.adjust(
    setNames(
        c(
            t.test(
                I(log10(value)) ~ group,
                data %>% filter(phage == "AC-L3")
            )$p.value,
            t.test(
                I(log10(value)) ~ group,
                data %>% filter(phage == "G10400")
            )$p.value,
            t.test(
                I(log10(value)) ~ group,
                data %>% filter(phage == "G2494")
            )$p.value,
            t.test(
                I(log10(value)) ~ group,
                data %>% filter(phage == "G9062")
            )$p.value,
            t.test(
                I(log10(value)) ~ group,
                data %>% filter(phage == "MM02")
            )$p.value,
            t.test(
                I(log10(value)) ~ group,
                data %>% filter(phage == "WFH")
            )$p.value
        ),
        c("AC-L3", "G10400", "G2494", "G9062", "MM02", "WFH")
    ),
    method = "holm"
)

p[p < 0.05]

ggsave(
    "Rplots-5.png",
    ggplot(data, aes(group, value, fill = group)) +
        stat_summary(fun = "mean", geom = "bar", show.legend = FALSE) +
        stat_summary(
            fun.data = "mean_cl_normal",
            geom = "errorbar",
            linewidth = 0.25,
            show.legend = FALSE
        ) +
        geom_jitter(position = position_jitter(height = 0, seed = 2)) +
        facet_grid(cols = vars(phage)) +
        scale_x_discrete(
            labels = c(
                "phage" = "phage",
                "phage + virucide" = "phage +<br>virucide"
            )
        ) +
        scale_y_continuous(
            transform = "log10",
            breaks = c(10^0, 10^3, 10^6, 10^9, 10^12),
            labels = label_log(base = 10),
            expand = expansion(mult = c(0, 0.05))
        ) +
        labs(y = "PFU/ml") +
        theme_bw(base_size = 8) +
        theme(
            text = element_text(family = "Arial"),
            axis.title.x = element_blank(),
            axis.text.x = element_markdown()
        ) +
        guides(fill = "none") +
        scale_fill_okabe_ito() +
        suppressWarnings(geom_signif(
            data = data.frame(
                phage = c("G10400", "G2494", "G9062", "MM02"),
                start = "phage",
                end = "phage + virucide",
                label = c("**", "**", "**", "**"),
                y = 13
            ),
            aes(xmin = start, xmax = end, annotations = label, y_position = y),
            manual = TRUE,
            inherit.aes = FALSE,
            size = 0.25,
            textsize = 8 * 0.8 / .pt,
            family = "Arial",
            tip_length = 0.025
        )),
    width = 7,
    height = 3.5,
    units = "in",
    dpi = 1200,
    device = ragg::agg_png()
)

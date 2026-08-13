library(dunn.test)

library(readxl)

library(dplyr)
library(stringr)
library(tidyr)

library(forcats)

library(ggplot2)

library(scales)

library(ggokabeito)
library(ggsignif)
library(ggtext)

library(patchwork)

library(ragg)

data1 <- bind_rows(
  read_excel(
    "44298_2026_222_MOESM2_ESM.xlsx",
    sheet = "MM02 invasion",
    range = "B33:G40",
    col_names = FALSE,
    .name_repair = "unique_quiet"
  ) %>%
    select(1, 6) %>%
    rename(name = `...1`, value = `...6`) %>%
    mutate(
      value = value -
        0.51 *
          read_excel(
            "44298_2026_222_MOESM2_ESM.xlsx",
            sheet = "MM02 invasion",
            range = "B45:G48",
            col_names = FALSE,
            .name_repair = "unique_quiet"
          ) %>%
            select(1, 6) %>%
            rename(name = `...1`, value = `...6`) %>%
            pull(value) %>%
            mean()
    ),
  read_excel(
    "44298_2026_222_MOESM2_ESM.xlsx",
    sheet = "MM02 invasion",
    range = "B50:N55",
    col_names = FALSE,
    .name_repair = "unique_quiet"
  ) %>%
    select(1, 6) %>%
    rename(name = `...1`, value = `...6`) %>%
    mutate(
      value = value -
        0.51 *
          read_excel(
            "44298_2026_222_MOESM2_ESM.xlsx",
            sheet = "MM02 invasion",
            range = "B56:G57",
            col_names = FALSE,
            .name_repair = "unique_quiet"
          ) %>%
            select(1, 6) %>%
            rename(name = `...1`, value = `...6`) %>%
            pull(value) %>%
            mean()
    )
) %>%
  mutate(
    name = fct_relevel(
      str_trim(str_remove(name, "^V[1-4]"), side = "left"),
      "MM02 + UPEC 8923 3h",
      "MM02 cells only 3h",
      "MM02 + UPEC 8923 24h",
      "MM02 cells only 24h"
    )
  )

data1 %>%
  filter(name == "MM02 + UPEC 8923 3h" | name == "MM02 cells only 3h") %>%
  count(name)

data1 %>%
  filter(name == "MM02 + UPEC 8923 3h" | name == "MM02 cells only 3h") %>%
  group_by(name) %>%
  summarise(SD = sd(value))

result1 <- t.test(
  log10(value) ~ fct_rev(name),
  data = data1 %>%
    filter(name == "MM02 + UPEC 8923 3h" | name == "MM02 cells only 3h")
)

round(unname(10^-diff(result1$estimate)), 2)

round(as.numeric(10^result1$conf.int), 2)

data1 %>%
  filter(name == "MM02 + UPEC 8923 24h" | name == "MM02 cells only 24h") %>%
  count(name)

data1 %>%
  filter(name == "MM02 + UPEC 8923 24h" | name == "MM02 cells only 24h") %>%
  group_by(name) %>%
  summarise(SD = sd(value))

result2 <- t.test(
  log10(value) ~ fct_rev(name),
  data = data1 %>%
    filter(name == "MM02 + UPEC 8923 24h" | name == "MM02 cells only 24h")
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

data2 <- bind_rows(
  read_excel(
    "44298_2026_222_MOESM2_ESM.xlsx",
    sheet = "P00 invasion",
    range = "B18:G25",
    col_names = FALSE,
    .name_repair = "unique_quiet"
  ) %>%
    select(1, 6) %>%
    rename(name = `...1`, value = `...6`) %>%
    mutate(
      value = value -
        0.51 *
          read_excel(
            "44298_2026_222_MOESM2_ESM.xlsx",
            sheet = "P00 invasion",
            range = "B29:G31",
            col_names = FALSE,
            .name_repair = "unique_quiet"
          ) %>%
            select(1, 6) %>%
            rename(name = `...1`, value = `...6`) %>%
            pull(value) %>%
            mean()
    ),
  read_excel(
    "44298_2026_222_MOESM2_ESM.xlsx",
    sheet = "P00 invasion",
    range = "B34:N39",
    col_names = FALSE,
    .name_repair = "unique_quiet"
  ) %>%
    select(1, 6) %>%
    rename(name = `...1`, value = `...6`) %>%
    mutate(
      value = value -
        0.51 *
          read_excel(
            "44298_2026_222_MOESM2_ESM.xlsx",
            sheet = "P00 invasion",
            range = "B40:G40",
            col_names = FALSE,
            .name_repair = "unique_quiet"
          ) %>%
            select(1, 6) %>%
            rename(name = `...1`, value = `...6`) %>%
            pull(value) %>%
            mean()
    )
) %>%
  mutate(
    name = fct_relevel(
      str_trim(
        str_remove(
          replace_when(
            name,
            name == "V3 G10400 + UPEC7958 3h" ~ "V3 G10400 + UPEC 7958 3h"
          ),
          "^V[1-4]"
        ),
        side = "left"
      ),
      "G10400 + UPEC 7958 3h",
      "G10400 cells only 3h",
      "G10400 + UPEC 7958 24h",
      "G10400 cells only 24h"
    )
  )

data2 %>%
  filter(
    name == "G10400 + UPEC 7958 3h" | name == "G10400 cells only 3h"
  ) %>%
  count(name)

data2 %>%
  filter(
    name == "G10400 + UPEC 7958 3h" | name == "G10400 cells only 3h"
  ) %>%
  group_by(name) %>%
  summarise(SD = sd(value))

result1 <- t.test(
  log10(value) ~ fct_rev(name),
  data = data2 %>%
    filter(
      name == "G10400 + UPEC 7958 3h" | name == "G10400 cells only 3h"
    )
)

round(unname(10^-diff(result1$estimate)), 2)

round(as.numeric(10^result1$conf.int), 2)

data2 %>%
  filter(
    name == "G10400 + UPEC 7958 24h" | name == "G10400 cells only 24h"
  ) %>%
  count(name)

data2 %>%
  filter(
    name == "G10400 + UPEC 7958 24h" | name == "G10400 cells only 24h"
  ) %>%
  group_by(name) %>%
  summarise(SD = sd(value))

result2 <- t.test(
  log10(value) ~ fct_rev(name),
  data = data2 %>%
    filter(
      name == "G10400 + UPEC 7958 24h" | name == "G10400 cells only 24h"
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

data1 <- data1 %>%
  mutate(
    group = case_when(
      grepl("3h", name) ~ "1",
      grepl("24h", name) ~ "2"
    )
  )

data2 <- data2 %>%
  mutate(
    group = case_when(
      grepl("3h", name) ~ "1",
      grepl("24h", name) ~ "2"
    )
  )

ggsave(
  "figure-1.png",
  (ggplot(
    data1 %>% filter(group == 1),
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
        "MM02 + UPEC 8923 3h" = "UPEC 8923 + &Phi; MM02",
        "MM02 cells only 3h" = "&Phi; MM02",
        "MM02 + UPEC 8923 24h" = "UPEC 8923 + &Phi; MM02",
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
        start = "MM02 + UPEC 8923 3h",
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
      data2 %>% filter(group == 1),
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
          "G10400 + UPEC 7958 3h" = "UPEC 7958 + &Phi; G10400",
          "G10400 cells only 3h" = "&Phi; G10400",
          "G10400 + UPEC 7958 24h" = "UPEC 7958 + &Phi; G10400",
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
          start = "G10400 + UPEC 7958 3h",
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
      data1 %>% filter(group == 2),
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
          start = "MM02 + UPEC 8923 24h",
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
        data2 %>% filter(group == 2),
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
            start = "G10400 + UPEC 7958 24h",
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

data3 <- read_excel(
  "44298_2026_222_MOESM2_ESM.xlsx",
  sheet = "MM02 invasion",
  range = "B2:G13",
  col_names = FALSE,
  .name_repair = "unique_quiet"
) %>%
  select(1, 6) %>%
  rename(name = `...1`, value = `...6`) %>%
  mutate(
    name = fct_relevel(
      sub("V[0-9] +(.+)", "\\1", name),
      "Cells + UPEC 8923 3h",
      "Cells + UPEC 8923 + Phage 3h",
      "Cells + UPEC 8923 24h",
      "Cells + UPEC 8923 + Phage 24/19h"
    )
  )

data3 %>%
  filter(
    name == "Cells + UPEC 8923 3h" | name == "Cells + UPEC 8923 + Phage 3h"
  ) %>%
  count(name)

data3 %>%
  filter(
    name == "Cells + UPEC 8923 3h" | name == "Cells + UPEC 8923 + Phage 3h"
  ) %>%
  group_by(name) %>%
  summarise(SD = sd(value))

data3 %>%
  filter(
    name == "Cells + UPEC 8923 24h" | name == "Cells + UPEC 8923 + Phage 24/19h"
  ) %>%
  count(name)

data3 %>%
  filter(
    name == "Cells + UPEC 8923 24h" | name == "Cells + UPEC 8923 + Phage 24/19h"
  ) %>%
  group_by(name) %>%
  summarise(SD = sd(value))

result2 <- t.test(
  log10(value) ~ fct_rev(name),
  data = data3 %>%
    filter(
      name == "Cells + UPEC 8923 24h" |
        name == "Cells + UPEC 8923 + Phage 24/19h"
    )
)

round(unname(10^-diff(result2$estimate)), 2)

round(as.numeric(10^result2$conf.int), 2)

result2$p.value

data4 <- read_excel(
  "44298_2026_222_MOESM2_ESM.xlsx",
  sheet = "P00 invasion",
  range = "B2:G13",
  col_names = FALSE,
  .name_repair = "unique_quiet"
) %>%
  select(1, 6) %>%
  rename(name = `...1`, value = `...6`) %>%
  mutate(
    name = fct_relevel(
      sub("V[0-9] +(.+)", "\\1", name),
      "Cells + UPEC 7958 3h",
      "Cells + UPEC 7958 + Phage 3h",
      "Cells + UPEC 7958 24 h",
      "Cells + UPEC 7958 + Phage 24/19 h"
    )
  )

data4 %>%
  filter(
    name == "Cells + UPEC 7958 3h" | name == "Cells + UPEC 7958 + Phage 3h"
  ) %>%
  count(name)

data4 %>%
  filter(
    name == "Cells + UPEC 7958 3h" | name == "Cells + UPEC 7958 + Phage 3h"
  ) %>%
  group_by(name) %>%
  summarise(SD = sd(value))

data4 %>%
  filter(
    name == "Cells + UPEC 7958 24 h" |
      name == "Cells + UPEC 7958 + Phage 24/19 h"
  ) %>%
  count(name)

data4 %>%
  filter(
    name == "Cells + UPEC 7958 24 h" |
      name == "Cells + UPEC 7958 + Phage 24/19 h"
  ) %>%
  group_by(name) %>%
  summarise(SD = sd(value))

result2 <- t.test(
  log10(value) ~ fct_rev(name),
  data = data4 %>%
    filter(
      name == "Cells + UPEC 7958 24 h" |
        name == "Cells + UPEC 7958 + Phage 24/19 h"
    )
)

round(unname(10^-diff(result2$estimate)), 2)

round(as.numeric(10^result2$conf.int), 2)

result2$p.value

data3 <- data3 %>%
  mutate(
    group = case_when(
      grepl("3h", name) ~ "1",
      grepl("24h|24/19h", name) ~ "2"
    )
  )

data4 <- data4 %>%
  mutate(
    group = case_when(
      grepl("3h", name) ~ "1",
      grepl("24 h|24/19 h", name) ~ "2"
    )
  )

ggsave(
  "figure-3.png",
  (ggplot(
    data3 %>% filter(group == 1),
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
        "Cells + UPEC 8923 3h" = "UPEC 8923",
        "Cells + UPEC 8923 + Phage 3h" = "UPEC 8923 + &Phi; MM02 SIM",
        "Cells + UPEC 8923 24h" = "UPEC 8923",
        "Cells + UPEC 8923 + Phage 24/19h" = "UPEC 8923 + &Phi; MM02 5h PI"
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
        start = "Cells + UPEC 8923 3h",
        end = "Cells + UPEC 8923 + Phage 3h",
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
      data4 %>% filter(group == 1),
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
          "Cells + UPEC 7958 3h" = "UPEC 7958",
          "Cells + UPEC 7958 + Phage 3h" = "UPEC 7958 + &Phi; G10400 SIM",
          "Cells + UPEC 7958 24 h" = "UPEC 7958",
          "Cells + UPEC 7958 + Phage 24/19 h" = "UPEC 7958 + &Phi; G10400 5h PI"
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
          start = "Cells + UPEC 7958 3h",
          end = "Cells + UPEC 7958 + Phage 3h",
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
      data3 %>% filter(group == 2),
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
          "Cells + UPEC 8923 3h" = "UPEC 8923",
          "Cells + UPEC 8923 + Phage 3h" = "UPEC 8923 + &Phi; MM02 SIM",
          "Cells + UPEC 8923 24h" = "UPEC 8923",
          "Cells + UPEC 8923 + Phage 24/19h" = "UPEC 8923 + &Phi; MM02 5h PI"
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
          start = "Cells + UPEC 8923 24h",
          end = "Cells + UPEC 8923 + Phage 24/19h",
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
        data4 %>% filter(group == 2),
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
            "Cells + UPEC 7958 3h" = "UPEC 7958",
            "Cells + UPEC 7958 + Phage 3h" = "UPEC 7958 + &Phi; G10400 SIM",
            "Cells + UPEC 7958 24 h" = "UPEC 7958",
            "Cells + UPEC 7958 + Phage 24/19 h" = "UPEC 7958 + &Phi; G10400 5h PI"
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
            start = "Cells + UPEC 7958 24 h",
            end = "Cells + UPEC 7958 + Phage 24/19 h",
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

data5 <- bind_rows(
  read_excel(
    "44298_2026_222_MOESM2_ESM.xlsx",
    sheet = "MM02 invasion",
    range = "B15:G22",
    col_names = FALSE,
    .name_repair = "unique_quiet"
  ) %>%
    select(1, 6) %>%
    rename(name = `...1`, value = `...6`),
  read_excel(
    "44298_2026_222_MOESM2_ESM.xlsx",
    sheet = "MM02 invasion",
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
      "UPEC 8923 1h",
      "MM02 after 1,5h",
      "UPEC 8923 0,5h",
      "MM02 after 1h"
    )
  )

data5 %>%
  filter(name == "UPEC 8923 1h" | name == "MM02 after 1,5h") %>%
  count(name)

data5 %>%
  filter(name == "UPEC 8923 1h" | name == "MM02 after 1,5h") %>%
  group_by(name) %>%
  summarise(SD = sd(value))

result1 <- t.test(
  log10(value) ~ fct_rev(name),
  data = data5 %>%
    filter(name == "UPEC 8923 1h" | name == "MM02 after 1,5h")
)

round(unname(10^-diff(result1$estimate)), 2)

round(as.numeric(10^result1$conf.int), 2)

data5 %>%
  filter(name == "UPEC 8923 0,5h" | name == "MM02 after 1h") %>%
  count(name)

data5 %>%
  filter(name == "UPEC 8923 0,5h" | name == "MM02 after 1h") %>%
  group_by(name) %>%
  summarise(SD = sd(value))

result2 <- t.test(
  log10(value) ~ fct_rev(name),
  data = data5 %>%
    filter(name == "UPEC 8923 0,5h" | name == "MM02 after 1h")
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
      grepl("UPEC 8923 1h|MM02 after 1,5h", name) ~ "1",
      grepl("UPEC 8923 0,5h|MM02 after 1h", name) ~ "2"
    )
  )

ggsave(
  "figure-4.png",
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
        "UPEC 8923 1h" = "UPEC 8923",
        "MM02 after 1,5h" = "UPEC 8923 + &Phi; MM02 1.5h PI",
        "UPEC 8923 0,5h" = "UPEC 8923",
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
        start = "UPEC 8923 1h",
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
          "UPEC 8923 1h" = "UPEC 8923",
          "MM02 after 1,5h" = "UPEC 8923 + &Phi; MM02 1.5h PI",
          "UPEC 8923 0,5h" = "UPEC 8923",
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
          start = "UPEC 8923 0,5h",
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

data6 <- bind_rows(
  read_excel(
    "44298_2026_222_MOESM2_ESM.xlsx",
    sheet = "Lysosome conditions",
    range = "B33:F35",
    col_names = FALSE,
    .name_repair = "unique_quiet"
  ) %>%
    select(1, 5) %>%
    rename(name = `...1`, value = `...5`) %>%
    mutate(condition = "1"),
  read_excel(
    "44298_2026_222_MOESM2_ESM.xlsx",
    sheet = "Lysosome conditions",
    range = "B39:F41",
    col_names = FALSE,
    .name_repair = "unique_quiet"
  ) %>%
    select(1, 5) %>%
    rename(name = `...1`, value = `...5`) %>%
    mutate(condition = "1"),
  read_excel(
    "44298_2026_222_MOESM2_ESM.xlsx",
    sheet = "Lysosome conditions",
    range = "B43:F45",
    col_names = FALSE,
    .name_repair = "unique_quiet"
  ) %>%
    select(1, 5) %>%
    rename(name = `...1`, value = `...5`) %>%
    mutate(condition = "1"),
  read_excel(
    "44298_2026_222_MOESM2_ESM.xlsx",
    sheet = "Lysosome conditions",
    range = "B48:F50",
    col_names = FALSE,
    .name_repair = "unique_quiet"
  ) %>%
    select(1, 5) %>%
    rename(name = `...1`, value = `...5`) %>%
    mutate(condition = "2"),
  read_excel(
    "44298_2026_222_MOESM2_ESM.xlsx",
    sheet = "Lysosome conditions",
    range = "B54:F56",
    col_names = FALSE,
    .name_repair = "unique_quiet"
  ) %>%
    select(1, 5) %>%
    rename(name = `...1`, value = `...5`) %>%
    mutate(condition = "2"),
  read_excel(
    "44298_2026_222_MOESM2_ESM.xlsx",
    sheet = "Lysosome conditions",
    range = "B58:F60",
    col_names = FALSE,
    .name_repair = "unique_quiet"
  ) %>%
    select(1, 5) %>%
    rename(name = `...1`, value = `...5`) %>%
    mutate(condition = "2")
) %>%
  mutate(
    name = fct_relevel(
      case_when(
        str_detect(
          str_squish(name),
          "PC [1-3] UPEC 7958 pH(7| 4,5)"
        ) ~ "PC UPEC 7958",
        str_detect(
          str_squish(name),
          "UPEC 7958 \\+ G10400 MOI 1 BR[1-3]"
        ) ~ "UPEC 7958 + Phage G10400 MOI 1",
        str_detect(
          str_squish(name),
          "UPEC 7958 \\+ G10400 MOI 0[\\.,]01 BR[1-3]"
        ) ~ "UPEC 7958 + Phage G10400 MOI 0.01"
      ),
      "PC UPEC 7958",
      "UPEC 7958 + Phage G10400 MOI 1",
      "UPEC 7958 + Phage G10400 MOI 0.01"
    ),
    condition = fct_relevel(condition, "1", "2")
  )

data7 <- bind_rows(
  read_excel(
    "44298_2026_222_MOESM2_ESM.xlsx",
    sheet = "Lysosome conditions",
    range = "B3:F5",
    col_names = FALSE,
    .name_repair = "unique_quiet"
  ) %>%
    select(1, 5) %>%
    rename(name = `...1`, value = `...5`) %>%
    mutate(condition = "1"),
  read_excel(
    "44298_2026_222_MOESM2_ESM.xlsx",
    sheet = "Lysosome conditions",
    range = "B9:F11",
    col_names = FALSE,
    .name_repair = "unique_quiet"
  ) %>%
    select(1, 5) %>%
    rename(name = `...1`, value = `...5`) %>%
    mutate(condition = "1"),
  read_excel(
    "44298_2026_222_MOESM2_ESM.xlsx",
    sheet = "Lysosome conditions",
    range = "B13:F15",
    col_names = FALSE,
    .name_repair = "unique_quiet"
  ) %>%
    select(1, 5) %>%
    rename(name = `...1`, value = `...5`) %>%
    mutate(condition = "1"),
  read_excel(
    "44298_2026_222_MOESM2_ESM.xlsx",
    sheet = "Lysosome conditions",
    range = "B18:F20",
    col_names = FALSE,
    .name_repair = "unique_quiet"
  ) %>%
    select(1, 5) %>%
    rename(name = `...1`, value = `...5`) %>%
    mutate(condition = "2"),
  read_excel(
    "44298_2026_222_MOESM2_ESM.xlsx",
    sheet = "Lysosome conditions",
    range = "B24:F26",
    col_names = FALSE,
    .name_repair = "unique_quiet"
  ) %>%
    select(1, 5) %>%
    rename(name = `...1`, value = `...5`) %>%
    mutate(condition = "2"),
  read_excel(
    "44298_2026_222_MOESM2_ESM.xlsx",
    sheet = "Lysosome conditions",
    range = "B28:F30",
    col_names = FALSE,
    .name_repair = "unique_quiet"
  ) %>%
    select(1, 5) %>%
    rename(name = `...1`, value = `...5`) %>%
    mutate(condition = "2")
) %>%
  mutate(
    name = fct_relevel(
      case_when(
        str_detect(
          str_squish(name),
          "PC [1-3]"
        ) ~ "PC UPEC 8923",
        str_detect(
          str_squish(name),
          "UPEC 8923 \\+ MM02 MOI 1 BR[1-3]"
        ) ~ "UPEC 8923 + Phage MM02 MOI 1",
        str_detect(
          str_squish(name),
          "UPEC 8923 \\+ MM02 MOI 0[\\.,]01 BR[1-3]"
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
      data = data6 %>%
        filter(
          name %in%
            c("PC UPEC 7958", "UPEC 7958 + Phage G10400 MOI 1") &
            condition == "1"
        )
    )$p.value,
    t.test(
      I(log10(value + 1)) ~ name,
      data = data6 %>%
        filter(
          name %in%
            c("PC UPEC 7958", "UPEC 7958 + Phage G10400 MOI 0.01") &
            condition == "1"
        )
    )$p.value,
    t.test(
      I(log10(value + 1)) ~ name,
      data = data6 %>%
        filter(
          name %in%
            c("PC UPEC 7958", "UPEC 7958 + Phage G10400 MOI 1") &
            condition == "2"
        )
    )$p.value,
    t.test(
      I(log10(value + 1)) ~ name,
      data = data6 %>%
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
  data = data6 %>%
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
      data = data7 %>%
        filter(
          name %in%
            c("PC UPEC 8923", "UPEC 8923 + Phage MM02 MOI 1") &
            condition == "1"
        )
    )$p.value,
    t.test(
      I(log10(value + 1)) ~ name,
      data = data7 %>%
        filter(
          name %in%
            c("PC UPEC 8923", "UPEC 8923 + Phage MM02 MOI 0.01") &
            condition == "1"
        )
    )$p.value,
    t.test(
      I(log10(value + 1)) ~ name,
      data = data7 %>%
        filter(
          name %in%
            c("PC UPEC 8923", "UPEC 8923 + Phage MM02 MOI 1") &
            condition == "2"
        )
    )$p.value,
    t.test(
      I(log10(value + 1)) ~ name,
      data = data7 %>%
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
  data = data7 %>%
    filter(
      name %in%
        c("PC UPEC 8923", "UPEC 8923 + Phage MM02 MOI 1") &
        condition == "1"
    )
)

t.test(
  I(log10(value + 1)) ~ name,
  data = data7 %>%
    filter(
      name %in%
        c("PC UPEC 8923", "UPEC 8923 + Phage MM02 MOI 0.01") &
        condition == "1"
    )
)

ggsave(
  "figure-5.png",
  (ggplot(data6, aes(name, value, fill = name)) +
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
    (ggplot(data7, aes(name, value, fill = name)) +
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

data8 <- read_excel(
  "44298_2026_222_MOESM2_ESM.xlsx",
  sheet = "LDH",
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

data8 %>% count(name)

data8 %>% group_by(name) %>% summarise(SD = sd(value))

dunn.test(
  data8 %>% pull(value),
  data8 %>% pull(name),
  method = "holm",
  list = TRUE,
  altp = TRUE
)

data9 <- read_excel(
  "44298_2026_222_MOESM2_ESM.xlsx",
  sheet = "LDH",
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
      "UPEC 7958 24h",
      "UPEC 7958 + G10400 24h",
      "G10400 24h",
      "Cells only 24h"
    )
  )

data9 %>% count(name)

data9 %>% group_by(name) %>% summarise(SD = sd(value))

dunn.test(
  data9 %>% pull(value),
  data9 %>% pull(name),
  method = "holm",
  list = TRUE,
  altp = TRUE
)

ggsave(
  "figure-6.png",
  (ggplot(data8, aes(name, value, fill = name)) +
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
    (ggplot(data9, aes(name, value, fill = name)) +
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
          "UPEC 7958 24h" = "UPEC 7958 24h",
          "UPEC 7958 + G10400 24h" = "UPEC 7958 + &Phi; G10400 24h",
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
        comparisons = list(c("UPEC 7958 24h", "UPEC 7958 + G10400 24h")),
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

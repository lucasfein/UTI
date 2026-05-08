library(readxl)
library(dplyr)
library(stringr)
library(ggplot2)
library(patchwork)
library(ragg)
library(ggtext)
library(forcats)
library(scales)
library(ggsignif)
library(ggokabeito)

data1 <- bind_rows(
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

data2 <- bind_rows(
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
      data = data1 %>%
        filter(
          name %in%
            c("PC UPEC 7958", "UPEC 7958 + Phage G10400 MOI 1") &
            condition == "1"
        )
    )$p.value,
    t.test(
      I(log10(value + 1)) ~ name,
      data = data1 %>%
        filter(
          name %in%
            c("PC UPEC 7958", "UPEC 7958 + Phage G10400 MOI 0.01") &
            condition == "1"
        )
    )$p.value,
    t.test(
      I(log10(value + 1)) ~ name,
      data = data1 %>%
        filter(
          name %in%
            c("PC UPEC 7958", "UPEC 7958 + Phage G10400 MOI 1") &
            condition == "2"
        )
    )$p.value,
    t.test(
      I(log10(value + 1)) ~ name,
      data = data1 %>%
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
  data = data1 %>%
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
      data = data2 %>%
        filter(
          name %in%
            c("PC UPEC 8923", "UPEC 8923 + Phage MM02 MOI 1") &
            condition == "1"
        )
    )$p.value,
    t.test(
      I(log10(value + 1)) ~ name,
      data = data2 %>%
        filter(
          name %in%
            c("PC UPEC 8923", "UPEC 8923 + Phage MM02 MOI 0.01") &
            condition == "1"
        )
    )$p.value,
    t.test(
      I(log10(value + 1)) ~ name,
      data = data2 %>%
        filter(
          name %in%
            c("PC UPEC 8923", "UPEC 8923 + Phage MM02 MOI 1") &
            condition == "2"
        )
    )$p.value,
    t.test(
      I(log10(value + 1)) ~ name,
      data = data2 %>%
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
  data = data2 %>%
    filter(
      name %in%
        c("PC UPEC 8923", "UPEC 8923 + Phage MM02 MOI 1") &
        condition == "1"
    )
)

t.test(
  I(log10(value + 1)) ~ name,
  data = data2 %>%
    filter(
      name %in%
        c("PC UPEC 8923", "UPEC 8923 + Phage MM02 MOI 0.01") &
        condition == "1"
    )
)

ggsave(
  "Rplots-6.png",
  (ggplot(data1, aes(name, value, fill = name)) +
    stat_summary(fun = "mean", geom = "bar", show.legend = FALSE) +
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
    suppressWarnings(geom_signif(
      data = data.frame(
        condition = "1",
        start = "PC UPEC 7958",
        end = "UPEC 7958 + Phage G10400 MOI 0.01",
        label = "**",
        y = 24,
        tip_length = 0.025
      ),
      aes(xmin = start, xmax = end, annotations = label, y_position = y),
      manual = TRUE,
      inherit.aes = FALSE,
      size = 0.25,
      textsize = 10 * 0.8 / .pt,
      family = "Arial"
    ))) /
    (ggplot(data2, aes(name, value, fill = name)) +
      stat_summary(fun = "mean", geom = "bar", show.legend = FALSE) +
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
      suppressWarnings(geom_signif(
        data = data.frame(
          condition = "1",
          start = c("PC UPEC 8923", "PC UPEC 8923"),
          end = c(
            "UPEC 8923 + Phage MM02 MOI 1",
            "UPEC 8923 + Phage MM02 MOI 0.01"
          ),
          label = c("*", "***"),
          y = c(22, 24),
          tip_length = 0.025
        ),
        aes(xmin = start, xmax = end, annotations = label, y_position = y),
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

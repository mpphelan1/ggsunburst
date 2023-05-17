move_rect <- function(df,stg){
  if(stg == 1){
    df0 <- filter(df, stage != stg)
    df1 <- filter(df, stage == stg) %>%
      mutate(
        x_min = pct_cum_min,
        x_max = pct_cum,
        multiplier = 1,
        m1 = 1)
    df <- bind_rows(df0, df1)
  }
  else{
    df0 <- filter(df, stage == stg -1) %>%
      mutate(mult1 = mult1/ct_prv,
             stage = stg,
             prv_seq = stage_seq,
             multiplier = pct*multiplier,
      ) %>%
      select(stage, prv_seq, multiplier,x_min, m1)
    df1 <- filter(df, stage == stg) %>% ## Want (stage, stage_seq, pct, pct_cum, pct_cum_min)
      select(-c(x_min, x_max, multiplier, m1)) %>%
      left_join(df0) %>%
      mutate(
        offset = multiplier*(m1 - mult1/ct_prv),
        across(starts_with('pct_'),~.x * multiplier),
        x1 = x_min,
        x_min = offset/2 + x1 + pct_cum_min,
        x_max = offset/2 + x1 + pct_cum) %>%
      select(-x1,-offset)
    df2 <- filter(df, stage != stg)
    df <- bind_rows(df1,df2)
  }
}

prep_for_plot <- function(df){
  max_levels = max(df$stage)
  for(i in 1:max_levels){
      df <- move_rect(df,i)
    }
  df$state <- substrRight(df$stage_seq)
  df
}

ggsunburst_basic <- function(df){
  df %>%
    ggplot() +
    geom_rect(aes(ymin = 0, ymax =1 , xmin = 0, xmax=1), fill = NA) +
    geom_rect(
      aes(
        ymin = stage + 1,
        ymax = stage + 2,
        xmin = x_min,
        xmax = x_max,
        fill = state),
      color = 'black',
    ) +
    labs(fill = "State\n") +
    coord_polar() +
    annotate(
      "text", label = paste0(max(df$mult1),"\nPatients"),
      x = 0, y = 0, size = 5, colour = "black"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      panel.grid = element_blank(),
      axis.title.x.bottom = element_blank(),
      axis.title.y.left = element_blank()
    )

}


## WITH TOOLTIP ##
# library(ggiraph)
ggsunburst_dynamic <- function(df){
  if(!('toolTip' %in% names(df))){
    df$toolTip = paste('Pathway:',df$stage_seq,
                       '\nCount:', df$ct)
  }

  df %>%
  mutate(pct_ring = ct/ringN
) %>%
  ggplot() +
  geom_rect(aes(ymin = 0, ymax =1 , xmin = 0, xmax=1), fill = NA) +
  ggiraph::geom_rect_interactive(
    aes(
      ymin = stage + 1,
      ymax = stage + 2,
      xmin = x_min,
      xmax = x_max,
      fill = state,
      tooltip = toolTip),
    color = 'black'
  ) +
  labs(fill = "State\n") +
  coord_polar() +
  annotate(
    "text", label = paste0(max(df$mult1),"\nPatients"),
    x = 0, y = 0, size = 5, colour = "black"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    panel.grid = element_blank(),
    axis.title.x.bottom = element_blank(),
    axis.title.y.left = element_blank()
  )
}
# girafe(ggobj = plt)


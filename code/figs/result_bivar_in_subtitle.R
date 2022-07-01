## result of the fit of a bivariant model in the subtitle
## it can be called inside the bivariant plot function
## ggtext package is required

result <- function (plot, fit, type) {
  
  if (type == "lm") {
    fit <- fit %>% 
      keep (names (.) %in% c("fstatistic", "adj.r.squared")) %>%
      flatten_dfc () %>%
      transmute (adjr = adj.r.squared,
                 pval = pf (q = value, df1 = numdf, df2 = dendf, lower.tail = F)) %>% 
      map2_dfc (c(2, 4),
                ~ round (.x, .y)) %>% 
      mutate (pval = if_else (pval == 0, "< 0.0001", paste ("=", pval)))
      
    plot <- plot +
      labs (subtitle = paste("*R*<sup>2</sup> =", fit$adjr, ", *p*", fit$pval)) +
      theme (plot.subtitle = element_markdown ())
  }
  
  if (type == "loess") {
    fit <- fit %>% 
      keep (names (.) == "s") %>% 
      flatten_dfc () %>% 
      mutate (s = round (s, 3))
    
    plot <- plot +
      labs (subtitle = paste0("RSE=", fit$s))
  }
  
  if (type == "gam") {
    fit <- fit %>% 
      keep (names (.) %in% c("s.table", "r.sq")) %>% 
      flatten_dfc () %>% 
      transmute (adjr = r.sq,
                 pval_s = V4) %>% 
      map2 (c(4, 2),
            ~ round (.x, .y))
    
    plot <- plot +
      labs (subtitle = paste0("*R*<sup>2</sup>=", fit$adjr, "<br>*p*=", fit$pval_s)) +
      theme (plot.subtitle = element_markdown ())
  }
  return (plot)
}

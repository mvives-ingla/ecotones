## boxplot + tukey test
## out --> outliers
## br --> breaks
## res --> show result of lm in the plot


boxplot <- function (data, x, y, fill = NULL, res = T, out = 0, br = 3, toupper = F,
                     outlier.shape = 19, ymax = NULL) {
  
  library (multcomp)
  library (tidyverse)
  library (emmeans)
  
  x <- ensym(x)
  y <- ensym(y)
  fill <- try(ensym(fill), silent = T)
  if(class(fill) == "try-error") {
    fill <- NULL
  }
  
  
  data <- data %>% 
    filter(!is.na(!!x),
          !is.na(!!y))
  
  if (nrow (data) == 0) {
    plot <- NULL
  }
    
  fo <- expr (!!y ~ !!x)
  
  fo.emm <- expr (~ !!x)
  if (is.null (ymax)){
    ymax <- data%>% 
      ungroup() %>% 
      summarize(max(!!y)*21/20) %>% 
      as.numeric()
  }
  
  ymin <- data %>% 
    ungroup() %>% 
    summarize(min(!!y)*19/20) %>% 
    as.numeric()
  
  if (is.null(fill)) {
    lm <- lm (formula = eval (fo), data = data) 
  
    lett <- lm %>%
      emmeans(eval(fo.emm)) %>% 
      multcomp::cld(method = "Tukey", Letters = letters)  %>%
      mutate(label = .group,
              .group = NULL,
              y = ymax)
  } else {
    lm <- data %>% 
      nest(data = -fill) %>% 
      mutate(lm = map(data, ~lm(formula = eval(fo), data = .)))
    
    lett <- lm %>%
      mutate(tuk = map(lm, ~emmeans(., eval(fo.emm))),
             tuk = map(tuk, ~multcomp::cld(., metod = "Tukey", Letters = letters)),
             tuk = map2(tuk,
                        list(ymax, ymin),
                        ~mutate(.x,
                                label = .group,
                                .group = NULL,
                                y = .y)), 
             tuk = map(tuk, dplyr::select, !!x, y, label)) %>%
      unnest(tuk)
    }
  
  
  
  if (toupper) {
    lett <- lett %>%
      mutate (label = toupper (label))
  }
  
  plot <- data %>% 
    ggplot (aes (x = !!x, y = !!y, fill = !!fill)) +
    geom_boxplot (outlier.shape = outlier.shape) +
    geom_text (data = lett, aes (x = !!x, y = y, label = label, color = !!fill)) +
    theme_classic ()
    
  
  if (res) {
  fit <- summary (lm)
  
  plot <- result (plot = plot, fit = fit, type = "lm")
  }
  
  return (plot)
}
    

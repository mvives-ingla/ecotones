## boxplot + tukey test
## out --> outliers
## br --> breaks
## res --> show result of lm in the plot

## repassar tidy evaluation, etc. per afegir arguments amb ... i acabar de solucionar problemes
## https://adv-r.hadley.nz/evaluation.html#quosure-dots

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
  
 #aconseguir fer-ho amb el procediment de multcomp
  # primer glht
  # segon cld
  # correccions de comparacions múltiples: https://cran.r-project.org/web/packages/afex/vignettes/afex_anova_example.html
  # https://cran.r-project.org/web/packages/multcomp/vignettes/multcomp-examples.pdf
  
  
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
    
  
  
  
  
  # if (out == 0) {
  #   df <- data %>%
  #     map (~ select (., patch_plot, yvar, species)) %>%
  #     map (~ mutate (., var = .[[yvar]]))
  #   
  #   lim <- df %>%
  #     map (~ group_by (., patch_plot)) %>%
  #     map (~ mutate (.,
  #                    Q75 = as.numeric (quantile (var, probs = 0.75, na.rm = T)),
  #                    Q25 = as.numeric (quantile (var, probs = 0.25, na.rm = T)),
  #                    IQR = IQR (var, na.rm = T),
  #                    group = case_when (var > Q75 + 1.5*IQR ~ 1,
  #                                       var >= Q75 ~ 2,
  #                                       var >= Q25 ~ 3,
  #                                       var >= Q25 - 1.5*IQR ~ 4,
  #                                       TRUE ~ 5))) %>%
  #     map (~ ungroup (.)) %>%
  #     map (~ filter (.,
  #                    group != 1,
  #                    group != 5)) %>%
  #     map (~ group_by (., patch_plot)) %>%
  #     map (~ summarise (.,
  #                       ymax = max (var, na.rm = T),
  #                       ymin = min (var, na.rm = T)))
  #   
  #   
  #   plot <- df %>%
  #     map (~ lm (formula = var ~ patch_plot, data = .)) %>%
  #     map (~ emmeans (., ~ patch_plot)) %>% 
  #     map (~ CLD (., method = "tukey", Letters = letters)) %>%
  #     map2 (lim,
  #           ~ data.frame (x = .x$patch_plot,
  #                         y = max (.y$ymax),
  #                         label = .x$.group)) %>%
  #     map2 (df,
  #           ~ ggplot (data = .y, aes (x = patch_plot, y = var)) +
  #             geom_boxplot (outlier.shape = NA,
  #                           fill = NA) +
  #             geom_text (data = .x, aes (x = x, y = y, label = label),
  #                        nudge_x = 0.1) +
  #             labs (title = .y$species,
  #                   y = yvar) +
  #             theme_classic ()) %>% 
  #     map2 (lim,
  #           ~ .x +
  #             scale_y_continuous (expand = expand_scale (mult = c(0.1, 0.2)),
  #                                 limits = c (min (.y$ymin),
  #                                             max (.y$ymax)),
  #                                 breaks = pretty_breaks (n = br + 1, min.n = br)))
  #   
  # }
    

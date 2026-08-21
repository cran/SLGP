## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----loadQuakes, warning=FALSE, message=FALSE---------------------------------
data("quakes")
library(tidyr)
library(dplyr)

# Bin data together for visualisation purpose
df <- quakes %>%
  mutate(long_bin = cut(long, breaks = seq(165, 190, by = 2.5), include.lowest = FALSE)) %>%
  group_by(long_bin) %>%
  mutate(long_bin = paste0(long_bin, "\nn=", n()))%>%
  ungroup()%>%
  mutate(long_bin = factor(long_bin, 
                           levels = sort(unique(long_bin), decreasing = FALSE))) %>%
  data.frame()

range_response <- c(4, 7) 
range_x <- c(165, 190)

## ----figureQuakesDiscrete, fig.cap = "A visual representation of the event magnitudes depending on the longitude in the `quakes` catalogue.", out.width = "0.99\\textwidth", fig.height=3.5, fig.width=10, fig.align='center',fig.pos="H"----
library(ggplot2)
library(ggpubr)
library(viridis)


scatter_plot <- ggplot(df, aes(x = long, y = mag)) +
  geom_point(alpha = 0.5, color = "navy") +
  labs(x = "Longitude (°)",
       y = "Magnitude",
       title = "Observed earthquake magnitudes")  +
  theme_bw()+
  coord_cartesian(xlim=range_x,
                  ylim=range_response)

# Compute normalized frequencies per long_bin
df_bar <- df %>%
  count(long_bin, mag) %>%
  group_by(long_bin) %>%
  mutate(prop = n / sum(n))

# Histogram: Distribution of mag by 'long' bin
hist_plot <- ggplot(df_bar, aes(x = mag)) +
  geom_bar(mapping=aes(y = prop), stat = "identity",
           fill = "darkgrey", color = "grey50", lwd = 0.18, alpha = 0.7)  +
  geom_rug(data = df, aes(x = mag), 
           sides = "b", color = "navy", alpha = 0.5) +
  facet_wrap(~ long_bin, scales = "free_y", nrow=2) +
  labs(x = "Magnitude", 
       y = "Probability density", 
       title = "Histogram of 'magnitude' by 'long' group") +
  theme_bw()+
  coord_cartesian(xlim=range_response,
                  ylim=c(0, 0.5)) 
ggarrange(scatter_plot, hist_plot, ncol = 2, nrow = 1,
          widths = c(0.3, 0.7))

## ----SLGPfitting--------------------------------------------------------------
library(SLGP)

modelMAP <- slgp(mag~long, # Use a formula with two indexing variables
                 data=df,
                 method="MAP", #Maximum a posteriori estimation scheme
                 basisFunctionsUsed = "RFF",
                 interpolateBasisFun="WNN", # Accelerate inference
                 hyperparams = list(lengthscale=c(0.1, 0.1), 
                                    sigma2=1), 
                 nIntegral = 31, 
                 sigmaEstimationMethod = "heuristic", 
                 # Set to heuristic for numerical stability                 
                 predictorsLower= c(range_x[1]),
                 predictorsUpper= c(range_x[2]),
                 responseRange= range_response,
                 opts_BasisFun = list(nFreq=150,
                                      MatParam=5/2),
                 discrete=TRUE)

## ----SLGPplottingPrior1, fig.cap = "Conditional magnitude probabilities across longitude under the MAP estimate of the SLGP.", fig.fullwidth=TRUE, fig.height=5, fig.width=10, fig.align='center', fig.pos="H"----
plot( modelMAP,
      newdata = data.frame(long = seq(range_x[1], range_x[2], length.out = 6)),
      draw = "mean",
      panels = TRUE,
      n_response = 31,
      discrete = TRUE)

## ----SLGPplottingMAPDiscrete, fig.cap = "Empirical magnitude distributions within longitude bins and SLGP MAP estimates at the corresponding bin centers.", fig.fullwidth = TRUE, fig.height = 4, fig.width = 10, fig.align = "center", fig.pos = "H"----

selected_values <- c(167, 180, 185)
gap <- 0.5

df_filtered <- df %>%
  mutate(interval=findInterval(long, c(0, 
                                       selected_values[1]-gap, 
                                       selected_values[1]+gap, 
                                       selected_values[2]-gap, 
                                       selected_values[2]+gap, 
                                       selected_values[3]-gap, 
                                       selected_values[3]+gap)))%>%
  filter(interval %in% c(2, 4, 6))%>%
  group_by(interval)%>%
  mutate(category = paste0("long close to ", c("", selected_values[1],
                                               "", selected_values[2],
                                               "", selected_values[3])[interval], 
                           "\nn=", n()))

names <- sort(unique(df_filtered$category))
dfGrid <- data.frame(expand.grid(selected_values, 
                                 seq(range_response[1], range_response[2],, 31)))
colnames(dfGrid) <- c("long", "mag")
predMAP <- predict(modelMAP, newdata = dfGrid, discrete=TRUE, nIntegral=31)

colnames(predMAP) <- c("long", "mag", "MAP estimator")
predMAP <- predMAP%>%
  pivot_longer(-c("long", "mag"))
predMAP$category <-ifelse(predMAP$long==selected_values[1], names[1],
                          ifelse(predMAP$long==selected_values[2], names[2], names[3]))


df_emp <- df_filtered %>%
  count(category, mag) %>%
  group_by(category) %>%
  mutate(prob = n / sum(n)) %>%
  ungroup()

ggplot(mapping=aes(x = mag)) +
  geom_col(data = df_emp, aes(y = prob), width = 0.09, fill = "darkgrey",
    color = "grey50", linewidth = 0.2, alpha = 0.7)+
  geom_step(data=predMAP, mapping=aes(y=value, group=name, col=name), 
            lwd=1.1, direction = "mid")+
  facet_wrap(~ category, scales = "free_y", nrow=1) +
  labs(x = "Magnitude",
       y = "Probability",    
       title = "Binned 'magnitude' histograms vs SLGP MAP estimates at bins centers") +
  theme_bw()+
  theme(legend.position="bottom",
        legend.direction = "horizontal",
        legend.title = element_blank())+
  coord_cartesian(xlim=range_response,
                  ylim=c(0, 0.2)) 


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

range_response <- c(40, 680) # Can use range(df$depth), or user defined range as we do here
range_x <- c(165, 190) # Can use range(df$long), or user defined range as we do here

## ----figureQuakes, fig.cap = "A visual representation of the dependency of hypocentre depth on longitude in the Fiji-Tonga `quakes` catalogue.", fig.fullwidth=TRUE, fig.height=4, fig.width=10, fig.align='center',fig.pos="H"----
library(ggplot2)
library(ggpubr)
library(viridis)
# Scatterplot: long vs depth
scatter_plot <- ggplot(df, aes(x = long, y = depth)) +
  geom_point(alpha = 0.5, color = "grey", pch=20) +
  labs(y = "Hypocentral depth (km)",
       x = "Longitude (°)",
       title = "Depth measured at various longitudes") +
  theme_bw()+
  coord_cartesian(xlim=range_x,
                  ylim=range_response)

# Histogram: Distribution of depth by 'long' bin
hist_plot <- ggplot(df, aes(x = depth)) +
  geom_histogram(mapping=aes(y=after_stat(density)),
                 position = "identity", breaks = seq(40, 680, 40),
                 fill="darkgrey", col="grey50", lwd=0.2, alpha=0.7) +
  geom_rug(sides = "b", color = "navy", alpha = 0.5)+
  facet_wrap(~ long_bin, scales = "free_y", nrow=2) +
  labs(x = "Longitude (°)",
       y = "Probability density", 
       title = "Histogram of 'depth' by 'longitude' group") +
  theme_bw()+
  coord_cartesian(xlim=range_response,
                  ylim=c(0, 0.015))
ggarrange(scatter_plot, hist_plot, ncol = 2, nrow = 1,
          widths = c(0.3, 0.7))

## ----SLGPfitting1-------------------------------------------------------------
library(SLGP)
modelMAP <- slgp(depth~long, # Use a formula to specify predictors VS response
                 # Can use depth~. for all variables,
                 # Or 'depth' ~ long + var2 + var3 for more variables
                 data=df,
                 method="MAP", # MAP estimation scheme
                 basisFunctionsUsed = "RFF",
                 interpolateBasisFun="WNN", # Will Accelerate inference
                 hyperparams = list(lengthscale=c(0.15, 0.15), 
                                    # Applied to normalised data
                                    # So 0.15 is 15% of the range of values
                                    sigma2=1), 
                 # Will be re-selected with sigmaEstimationMethod
                 sigmaEstimationMethod = "heuristic", 
                 # Set to heuristic for numerical stability                 
                 predictorsLower= c(range_x[1]),
                 predictorsUpper= c(range_x[2]),
                 responseRange= range_response,
                 opts_BasisFun = list(nFreq=200,
                                      MatParam=5/2),
                 seed=1)

## ----SLGPplottingPrior1, fig.cap = "Conditional depth densities across longitude under the MAP estimate of the SLGP.", fig.fullwidth=TRUE, fig.height=5, fig.width=10, fig.align='center', fig.pos="H"----
plot( modelMAP,
      newdata = data.frame(long = seq(range_x[1], range_x[2], length.out = 6)),
      draw = "mean",
      panels = TRUE,
      n_response = 101,
      discrete = FALSE)

## ----SLGPplotting2, fig.cap = "Predictive probability density of 'depth' at 'long', seen over slices.", fig.fullwidth=TRUE, fig.height=5, fig.width=10, fig.align='center',fig.pos="H"----
library(viridis)
dfGrid <- data.frame(expand.grid(seq(range_x[1], range_x[2], 1), 
                                 seq(range_response[1], range_response[2],, 101)))
colnames(dfGrid) <- c("long", "depth")
pred <- predict(modelMAP, newdata= dfGrid)

scale_factor <- 300
ggplot()  +
  labs(x = "Hypocentral depth (km)",
       y = "Longitude (°)")+
  theme_bw()+
  geom_ribbon(data=pred,
              mapping=aes(x=depth, ymax=scale_factor*pdf_1+long, 
                          ymin=long, group=-long, fill=long),
              col="grey", alpha=0.9)+
  geom_point(data=df,
             mapping=aes(x = depth, y = long), alpha = 0.5, 
             pch=20, color = "grey")+
  scale_fill_viridis(option = "plasma",
                     guide = guide_colorbar(title = "Indexing variable: Longitude",
                                            barheight = unit(2, units = "mm"),
                                            barwidth = unit(55, units = "mm"),
                                            title.position = 'top',
                                            label.position = "bottom",
                                            title.hjust = 0.5))+
  theme(legend.position = "bottom")+
  coord_flip()

## ----SLGPplottingMAP, fig.cap = "Binned 'depth' histograms by 'long' (width = 1) VS SLGP MAP estimators at bins centers", fig.fullwidth=TRUE, fig.height=4, fig.width=10, fig.align='center',fig.pos="H"----

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
                                 seq(range_response[1], range_response[2],, 101)))
colnames(dfGrid) <- c("long", "depth")
predMAP <- predict(modelMAP, newdata = dfGrid)
colnames(predMAP) <- c("long", "depth", "MAP estimator")
predMAP <- predMAP%>%
  pivot_longer(-c("long", "depth"))
predMAP$category <-ifelse(predMAP$long==selected_values[1], names[1],
                          ifelse(predMAP$long==selected_values[2], names[2], names[3]))

ggplot(mapping=aes(x = depth)) +
  geom_histogram(df_filtered,
                 mapping=aes(y=after_stat(density)),
                 position = "identity", breaks = seq(40, 680, 20),
                 fill="darkgrey", col="grey50", lwd=0.2, alpha=0.7) +
  geom_rug(data=df_filtered, sides = "b", color = "navy", alpha = 0.5)+
  geom_line(data=predMAP, mapping=aes(y=value, group=name), 
            color = "black", lwd=0.1, alpha=0.5)+
  geom_line(data=predMAP, mapping=aes(y=value, group=name, col=name), lwd=1.1)+
  facet_wrap(~ category, scales = "free_y", nrow=1) +
  labs(x = "Hypocentral depth (km)",
       y = "Probability density",
       title = "Binned 'depth' histograms by 'long' (width = 1) VS SLGP-MAP estimators at bins centers") +
  theme_bw()+
  theme(legend.position="bottom",
        legend.direction = "horizontal",
        legend.title = element_blank())+
  coord_cartesian(xlim=range_response,
                  ylim=c(0, 0.02))

## ----SLGPfitting2-------------------------------------------------------------
modelLaplace <- update(modelMAP, 
                       newdata = df, 
                       method="Laplace")

## ----SLGPfitting22, eval=FALSE------------------------------------------------
#  # Or equivalent, more explicit in the re-using of the elements
#  # From the SLGP prior
#  modelLaplace <- slgp(medv~age,
#                       data=df,
#                       method="Laplace", #Maximum a posteriori estimation scheme
#                       basisFunctionsUsed = "RFF",
#                       interpolateBasisFun="WNN", # Accelerate inference
#                       hyperparams = modelMAP@hyperparams,
#                       sigmaEstimationMethod = "none",# Already selected in the prior
#                       predictorsLower= c(range_x[1]),
#                       predictorsUpper= c(range_x[2]),
#                       responseRange= range_response,
#                       opts_BasisFun = modelMAP@opts_BasisFun,
#                       BasisFunParam = modelMAP@BasisFunParam,
#                       seed=1)

## ----SLGPLaplaceplo0, fig.cap = "Predictive probability density of 'depth' at 'long', as predicted by a SLGP with Laplace approximation", fig.fullwidth=TRUE, fig.height=5, fig.width=10, fig.align='center', fig.pos="H"----
plot( modelLaplace,
      newdata = data.frame(long = seq(range_x[1], range_x[2], length.out = 6)),
      draw = c("mean", 1:10),
      panels = TRUE,
      n_response = 101,
      discrete = FALSE)

## ----SLGPLaplaceplot, fig.cap = "Predictive probability density (and draws from a Laplace approximation) of 'depth' at 'long', seen over 3 slices.", fig.fullwidth=TRUE, fig.height=4, fig.width=8, fig.align='center',fig.pos="H"----

dfGrid <- data.frame(expand.grid(seq(3), 
                                 seq(range_response[1], range_response[2],, 101)))
colnames(dfGrid) <- c("ID", "depth")
dfGrid$long <- selected_values[dfGrid$ID]
pred <- predict(modelLaplace, newdata = dfGrid)
pred$meanpdf <- rowMeans(pred[, -c(1:3)])

library(tidyr)
# Filter the data: keep values within ±5 of the selected ones
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
pred$category <- names[pred$ID]

set.seed(1)
selected_cols <- sample(seq(1000), size=10, replace=FALSE)
df_plot <- pred %>%
  dplyr::select(c("long", "depth", "category", 
                  paste0("pdf_", selected_cols)))%>%
  pivot_longer(-c("long", "depth", "category"))


ggplot(mapping=aes(x = depth)) +
  geom_histogram(df_filtered,
                 mapping=aes(y=after_stat(density)),
                 position = "identity", breaks = seq(40, 680, 20),
                 fill="darkgrey", col="grey50", lwd=0.2, alpha=0.7) +
  geom_rug(data=df_filtered, sides = "b", color = "navy", alpha = 0.5)+
  geom_line(data=df_plot, mapping=aes(y=value, group=name), 
            color = "black", lwd=0.1, alpha=0.5)+
  geom_line(data=pred, mapping=aes(y=meanpdf, group=category), color = "red")+
  facet_wrap(~ category, scales = "free_y", nrow=1) +
  labs(x = "Hypocentral depth (km)",
       y = "Probability density",
       title = "Binned 'depth' histograms by 'long' (width = 1) VS SLGP-Laplace estimators at bins centers") +
  theme_bw()+
  coord_cartesian(xlim=range_response,
                  ylim=c(0, 0.02))


## ----SLGPfitting3, eval=FALSE-------------------------------------------------
#  modelMCMC <- slgp(medv~age, # Use a formula to specify predictors VS response
#                    # Can use medv~. for all variables,
#                    # Or medv ~ age + var2 + var3 for more variables
#                    data=df,
#                    method="MCMC", #MCMC
#                    basisFunctionsUsed = "RFF",
#                    interpolateBasisFun="WNN", # Accelerate inference
#                    hyperparams = list(lengthscale=c(0.15, 0.15),
#                                       # Applied to normalised data
#                                       # So 0.15 is 15% of the range of values
#                                       sigma2=1),
#                    # Will be re-selected with sigmaEstimationMethod
#                    sigmaEstimationMethod = "heuristic", # Set to heuristic for numerical stability
#                    predictorsLower= c(range_x[1]),
#                    predictorsUpper= c(range_x[2]),
#                    responseRange= range_response,
#                    opts_BasisFun = list(nFreq=100,
#                                         MatParam=5/2),
#                    opts = list(stan_chains=2, stan_iter=1000))

## ----SLGPfitting3b, eval=FALSE------------------------------------------------
#  modelMCMC <- update(modelMAP,
#                      newdata = df,
#                      method="MCMC")

## ----SLGPMCMCplotMoments, fig.cap = "Simultaneous prediction of the fields moments (and associated uncertainty) using a SLGP model", fig.fullwidth=TRUE, fig.height=4, fig.width=10, fig.align='center',fig.pos="H"----
dfX <- data.frame(long=seq(range_x[1], range_x[2], 1))

predMean <- predict(modelLaplace, type= "moments", newdata = dfX, 
                    power=c(1),
                    centered=FALSE) # Uncentered moments
# For the mean
predVar <- predict(modelLaplace, type= "moments", newdata = dfX, 
                   power=c(2, 3, 4),
                   centered=TRUE) # Centered moments
# For the variance, Kurtosis and Skewness

pred <- rbind(predMean, predVar)
pred <- pred %>%
  pivot_longer(-c("long", "power"))%>%
  mutate(value=ifelse(power==2, sqrt(value), value))%>% # Define std
  pivot_wider(values_from = value,
              names_from = power)%>%
  mutate(`3`=`3`/`2`^2,
         `4`=`4`/`2`^4)%>% # Kurtosis and Skewness
  pivot_longer(-c("long", "name"), names_to = "power")%>%
  data.frame()

pred$power <- factor(c("Expected value", 
                       "Standard deviation",
                       "Skewness","Kurtosis")[as.numeric(pred$power)],
                     levels=c("Expected value", "Standard deviation",
                              "Skewness", "Kurtosis"))
df_plot <- pred %>%
  group_by(long, power)%>%
  summarise(q10 = quantile(value, probs=c(0.1)),
            q50 = quantile(value, probs=c(0.5)),
            q90 = quantile(value, probs=c(0.9)),
            mean = mean(value), .groups="keep")%>%
  ungroup() # summarise uncertainty


ggplot(df_plot, mapping=aes(x = long, group=power)) +
  geom_ribbon(mapping = aes(ymin=q10, ymax=q90),
              alpha = 0.25, lty=2, col="black", fill="cornflowerblue")+
  geom_line(mapping=aes(y=q50))+
  facet_wrap(.~power,
             scales = "free", nrow=1)+
  labs(x =  "Longitude (°)", 
       y = "Moment value") +
  theme_bw()+
  coord_cartesian(xlim=range_x)


## ----SLGPMCMCplotQuantiles2, fig.cap = "Simultaneous quantile prediction (and associated uncertainty) using a SLGP model, with estimation performed by MCMC", fig.fullwidth=TRUE, fig.height=5, fig.width=10, fig.align='center',fig.pos="H", warning=FALSE, message=FALSE----
probsL <-  c(10, 25, 50, 75, 90)/100

# SLGP prediction
pred <- predict(modelLaplace, type = "quantiles",
                newdata = dfX,
                probs = probsL)

df_plot <- pred %>%
  pivot_longer(-c("long", "probs"))%>%
  group_by(long, probs)%>%
  summarise(q10 = quantile(value, probs=c(0.1)),
            q50 = quantile(value, probs=c(0.5)),
            q90 = quantile(value, probs=c(0.9)),
            mean = mean(value), .groups="keep")%>%
  ungroup()%>%
  mutate(probs=factor(paste0("Quantile: ", 100*probs, "%"),
                      levels=paste0("Quantile: ", 100*probsL, "%")))


plot1 <-ggplot(df_plot, mapping=aes(x = long)) +
  geom_ribbon(mapping = aes(ymin=q10, ymax=q90, 
                            col=probs, fill=probs, group=probs),
              alpha = 0.25, lty=3)+
  geom_line(mapping=aes(y=q50, lty="SLGP", 
                        col=probs,group=probs), lwd=0.75)+
  labs(x = "Longitude [°E]",
       y = "Hypocentre depth [km]",
       col = "Quantile levels",
       fill = "Quantile levels",
       lty = "Quantile estimation method") +
  theme_bw() +
  coord_cartesian(xlim = range_x, ylim = range_response) 

plot2 <- ggplot(df, mapping = aes(x = long)) +
  geom_histogram(col = "navy", fill = "grey", alpha = 0.3,
                 breaks = seq(range_x[1], range_x[2], 1)) +
  theme_minimal() +
  labs(x = NULL, y = "Sample\ncount")


p_with_marginal <- ggarrange(
  plot2, plot1,
  ncol = 1, heights = c(1, 3),  # adjust height ratio
  align = "v"
)

# Display the plot
print(p_with_marginal)


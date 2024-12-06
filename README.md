# Distance based permutational analysis of variance for R

This code implements a nested analysis of variance via permutations based on a distance measure between observations.

Essentially it implements the model `outcome ~ A + A:B` where `A` is a fixed main effect (e.g. a grouping) and `B`
is either a fixed or random effect nested within `A`.

Code to run is as follows:

```{r}
source('R/permanova_nested.R')

library(tidyverse)
library(palmerpenguins)
dat <- penguins |> na.omit()

A <- dat$species |> as.factor()
B <- dat$sex |> as.factor()

d <- dat[,c("flipper_length_mm", "body_mass_g")] |>
  dist()

dist2 <- as.matrix(d)^2

species_effect <- perm_main(dist2, A, B)
sex_in_species_effect <- perm_nested(dist2, A, B)              

ss <- get_ss_main(dist2, A, B)   

as.data.frame(ss) |>
  mutate(across(everything(), ~ . / SS_tot * 100))

species_effect
sex_in_species_effect
```


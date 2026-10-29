source("code/R/taxon_display.R")
library(ggplot2)

keys <- c("S. pulcher", "S. purpuratus", "P. interruptus")
stopifnot(identical(taxon_display(keys), c("B. pulcher", keys[-1])))
stopifnot(identical(taxon_display("Semicossyphus pulcher - Biomass"),
                    "Bodianus pulcher - Biomass"))
stopifnot(identical(taxon_display("S. pulcher lnRR"), "B. pulcher lnRR"))
stopifnot(identical(taxon_display(expression(italic("S. pulcher") ~ "lnRR")),
                    expression(italic("B. pulcher") ~ "lnRR")))
stopifnot(identical(taxon_display(expression(italic("S. purpuratus") ~ "lnRR")),
                    expression(italic("S. purpuratus") ~ "lnRR")))

d <- data.frame(taxon = factor(keys, levels = keys), value = c(1, 2, 3))
before <- d
p <- ggplot(d, aes(taxon, value)) + geom_point() +
  scale_x_discrete(labels = taxon_display)
b <- ggplot_build(p)
stopifnot(identical(b$layout$panel_params[[1]]$x$get_labels(),
                    c("B. pulcher", keys[-1])))
stopifnot(identical(d, before))
stopifnot(identical(levels(d$taxon), keys))
stopifnot(identical(b$data[[1]]$y, d$value))

labels <- labeller(taxon = taxon_display)(data.frame(taxon = keys))
stopifnot(identical(unname(labels[[1]]), c("B. pulcher", keys[-1])))
cat("PASS: taxon display labels, preserved source keys and unchanged plot values.\n")

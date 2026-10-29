source("code/R/temporal_uncertainty_utils.R")
input <- read.csv("data/harmonized/harmonized_response_ratios.csv")
primary <- read.csv("data/sumstats_final.csv")
d <- temporal_primary_data(input, primary)
stopifnot(nrow(d) == 1750L, sum(d$resp == "Bio") == 831L,
  length(unique(d$history[d$resp == "Bio"])) == 63L,
  all(d$resp %in% c("Bio", "Den")), all(d$BA == "After"), all(d$source != "Landsat"))
stopifnot(inherits(try(temporal_primary_data(rbind(input, input[1, ]), primary), silent = TRUE), "try-error"))
satellite <- primary[primary$Source == "Landsat", ][1, ]
fake <- input[1, ]
fake$CA_MPA_Name_Short <- satellite$MPA
fake$y <- "Macrocystis pyrifera"
fake$resp <- "Bio"; fake$source <- "Landsat"
fake$BA <- "After"; fake$time <- 1; fake$year <- 3000
stopifnot(nrow(temporal_primary_data(rbind(input, fake), primary)) == nrow(d))
display <- d[d$resp == "Bio" & d$time <= 15, ]
stopifnot(nrow(display) == 723L,
  setequal(unique(display$history), unique(d$history[d$resp == "Bio"])))
cache <- readRDS("data/cache/temporal_primary_models.rds")
for (response in c("bio", "den")) {
  m <- cache[[response]]
  stopifnot(!lme4::isSingular(m), m@optinfo$conv$opt == 0L)
  sp <- levels(m@frame$Species)[2]
  z <- temporal_contrast(m, sp, slope = TRUE)
  nd <- data.frame(Species = factor(sp, levels = levels(m@frame$Species)), time = c(0, 1))
  x <- model.matrix(~ Species * time, nd)
  g <- x[2, ] - x[1, ]
  b <- lme4::fixef(m)
  stopifnot(abs(z$estimate - sum(g * b)) < 1e-12,
    abs(z$se - as.numeric(sqrt(drop(t(g) %*% vcov(m) %*% g)))) < 1e-12,
    is.finite(z$df), z$df > 0)
  grid <- temporal_prediction_grid(m, ifelse(response == "bio", "Bio", "Den"), times = c(0, 11, 15))
  stopifnot(nrow(grid) == length(levels(m@frame$Species)) * 3L,
    all(is.finite(grid$ci_low)), all(grid$ci_low < grid$ci_high))
}
cat("PASS: exact temporal primary cohort, nonsingular fits and contrast-specific slopes/predictions.\n")

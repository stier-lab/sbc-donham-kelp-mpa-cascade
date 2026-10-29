siting_tracking <- function(x, y) {
  if (length(x) < 4L || sd(x) < 1e-12 || sd(y) < 1e-12) return(NA_real_)
  cor(x, y)
}

siting_window <- function(d, anchor, width) {
  pre <- d[d$year >= anchor-width & d$year < anchor, ]
  post <- d[d$year > anchor & d$year <= anchor+width, ]
  eligible <- nrow(pre) >= ceiling(.8*width) && nrow(post) >= ceiling(.8*width)
  means <- c(mean(pre$inside), mean(pre$reference), mean(post$inside), mean(post$reference))
  positive <- eligible && all(is.finite(means) & means > 0)
  data.frame(anchor, width, n_before=nrow(pre), n_after=nrow(post), eligible,
    before_years=paste(pre$year,collapse=";"), after_years=paste(post$year,collapse=";"),
    before_inside=means[1], before_reference=means[2],
    after_inside=means[3], after_reference=means[4],
    change_log=if(eligible) mean(post$contrast)-mean(pre$contrast) else NA_real_,
    change_index=if(eligible) mean(post$inside-post$reference)-mean(pre$inside-pre$reference) else NA_real_,
    positive_means=positive,
    change_period_log=if(positive) log(means[3]/means[4])-log(means[1]/means[2]) else NA_real_)
}

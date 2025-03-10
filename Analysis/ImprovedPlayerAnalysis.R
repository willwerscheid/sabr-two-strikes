library(tidyverse) 
theme_set(theme_minimal())
library(mgcv)
library(scam)
source("Code/est_bs_dist.R")

fig_path <- "Writing/images/"
dat_path <- "Writing/data/"
my_ggsave <- function(fname, width = 4, height = 3) {
  ggsave(paste0(fig_path, fname), width = width, height = height)
}

statcast_split <- readRDS(paste0(dat_path, "statcast_split.rds"))
all_bs_dist <- readRDS(paste0(dat_path, "all_bs_dist.rds"))
player_map <- map_chr(statcast_split, list("player_name", 1))

statcast_split <- map(
  statcast_split,
  \(df) df |> 
    mutate(pitch_cat_2s = fct_cross(pitch_cat, two_strikes),
           .after = pitch_cat)
)

statcast_split <- map2(
  statcast_split, all_bs_dist,
  \(df, dist) df |> 
    mutate(
      relative_bs = qnorm(
        my_pskt(bat_speed, dist$df1, dist$gamma1, dist$mu1, dist$s1)
      ),
      .after = bat_speed
    )
)

statcast <- bind_rows(statcast_split)
eg_idx <- str_which(player_map, "Yordan")
eg_df <- statcast_split[[eg_idx]] 

min_x <- -3
max_x <- 3
min_z <- -1.25
max_z <- 5.5
min_z_rel <- -1.75
max_z_rel <- 2.25

min_x_plot <- -1.5
max_x_plot <- 1.5
min_z_rel_plot <- -0.5
max_z_rel_plot <- 1.5

mean_sz_bot <- mean(statcast$sz_bot)
mean_sz_top <- mean(statcast$sz_top)

sz <- tibble(
  x = c(-0.71, 0.71, 0.71, -0.71),
  xend = c(0.71, 0.71, -0.71, -0.71),
  y = c(mean_sz_top, mean_sz_top, mean_sz_bot, mean_sz_bot),
  yend = c(mean_sz_top, mean_sz_bot, mean_sz_bot, mean_sz_top)
)
relative_sz <- tibble(
  x = c(-0.71, 0.71, 0.71, -0.71),
  xend = c(0.71, 0.71, -0.71, -0.71),
  y = c(1, 1, 0, 0),
  yend = c(1, 0, 0, 1)
)

eg_sz_bot <- mean(eg_df$sz_bot)
eg_sz_top <- mean(eg_df$sz_top)

eg_sz <- tibble(
  x = c(-0.71, 0.71, 0.71, -0.71),
  xend = c(0.71, 0.71, -0.71, -0.71),
  y = c(eg_sz_top, eg_sz_top, eg_sz_bot, eg_sz_bot),
  yend = c(eg_sz_top, eg_sz_bot, eg_sz_bot, eg_sz_top)
)

##### Contact ----------

all_contact <- gam(is_contact ~ s(plate_x_flipped, plate_z_relative, k = 100),
                   data = statcast,
                   weights = p_full,
                   family = binomial)
# Not different enough to be worth the trouble:
# all_contact2 <- gam(is_contact ~ s(plate_x_flipped, plate_z_relative, 
#                                    k = 100, by = pitch_cat),
#                     data = statcast,
#                     weights = p_full,
#                     family = binomial)

pred_df <- expand.grid(
  plate_x_flipped = seq(min_x, max_x, by = 0.05), 
  plate_z_relative = seq(min_z_rel, max_z_rel, by = 0.05)
)
pred_df <- pred_df |>
  mutate(
    plate_x = plate_x_flipped,
    plate_z = plate_z_relative * (mean_sz_top - mean_sz_bot) + mean_sz_bot,
    contact_pred = predict(all_contact, newdata = pred_df, type = "response")
  )

ggplot(pred_df) +
  geom_tile(aes(x = plate_x, y = plate_z, fill = contact_pred)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = mean(statcast$is_contact, na.rm = TRUE)) +
  geom_segment(data = sz, aes(x = x, y = y, xend = xend, yend = yend)) +
  labs(x = "plate_x (RHB)", y = "plate_z", fill = "Prob. contact")
my_ggsave("all_contact.pdf")


n_grid_pts <- 25
# bound_df <- tibble(
#   plate_x_flipped = c(rep(seq(min_x, max_x, length.out = n_grid_pts), times = 2), 
#                       rep(min_x, n_grid_pts), rep(max_x, n_grid_pts)),
#   plate_z_relative = c(rep(min_z_rel, n_grid_pts), rep(max_z_rel, n_grid_pts), 
#                        rep(seq(min_z_rel, max_z_rel, length.out = n_grid_pts), times = 2)),
#   is_contact = 0,
#   p_full = 1
# ) 
# bound_df <- bound_df |>
#   crossing(pitch_cat_2s = factor(levels(eg_df$pitch_cat_2s)))

even_prior_df <- expand.grid(
  plate_x_flipped = seq(min_x, max_x, length.out = n_grid_pts), 
  plate_z_relative = seq(min_z_rel, max_z_rel, length.out = n_grid_pts)
) 
quantile_prior_df <- expand.grid(
  plate_x_flipped = quantile(
    statcast$plate_x_flipped, seq(.005, .995, length.out = n_grid_pts)
  ), 
  plate_z_relative = quantile(
    statcast$plate_z_relative, seq(.005, .995, length.out = n_grid_pts)
  )
)
augment_prior_df <- function(prior_df) {
  prior_df <- prior_df |>
    mutate(contact_pred = predict(
      all_contact, newdata = prior_df, type = "response"
    ))
  prior_df <- prior_df |>
    crossing(is_contact = 0:1) |>
    mutate(wt = is_contact * contact_pred + (1 - is_contact) * (1 - contact_pred))
  return(prior_df)
}
even_prior_df <- augment_prior_df(even_prior_df)
quantile_prior_df <- augment_prior_df(quantile_prior_df)

test_contact_mod <- function(player_df, 
                             n_pseudopitches, 
                             pseudopitch_spacing = c("even", "quantile"),
                             pitch_cat_model = c("smooth_by", "fixed_only", "ignore"),
                             seed = 666) {
  cat(player_df$player_name[1], n_pseudopitches, "\n")
  
  pseudopitch_spacing <- match.arg(pseudopitch_spacing)
  pitch_cat_model <- match.arg(pitch_cat_model)
  
  set.seed(seed)
  train_idx <- sample(1:nrow(player_df), size = floor(0.8 * nrow(player_df)))
  test_idx <- setdiff(1:nrow(player_df), train_idx)
  
  train_df <- player_df |> slice(train_idx) 
  test_df <- player_df |> slice(test_idx)
  
  if (pseudopitch_spacing == "even") {
    prior_df <- even_prior_df
  } else {
    prior_df <- quantile_prior_df
  }
  
  if (is.infinite(n_pseudopitches)) {
    test_df <- test_df |>
      mutate(p_contact = predict(all_contact, newdata = test_df, type = "response")) 
  } else {
    # if (n_pseudopitches < 0) {
    #   train_df <- train_df |> bind_rows(bound_df)
    # } else 
    if (n_pseudopitches > 0) {
      prior_wt <- n_pseudopitches / (nrow(prior_df) / 2)
      prior_df <- prior_df |>
        mutate(p_full = prior_wt * wt) 
      prior_df <- prior_df |>
        crossing(two_strikes = factor(levels(player_df$two_strikes)))
      if (pitch_cat_model != "ignore") {
        prior_df <- prior_df |>
          crossing(pitch_cat = factor(levels(player_df$pitch_cat))) |>
          mutate(pitch_cat_2s = fct_cross(pitch_cat, two_strikes))
      }
      train_df <- train_df |>
        bind_rows(prior_df) |>
        filter(p_full > 0)
    }
    if (pitch_cat_model == "smooth_by") {
      gmod <- gam(
        is_contact ~ s(plate_x_flipped, plate_z_relative, by = pitch_cat_2s) +
          pitch_cat + two_strikes, 
        data = train_df, 
        family = binomial, 
        weights = p_full
      )
    } else if (pitch_cat_model == "fixed_only") {
      gmod <- gam(
        is_contact ~ s(plate_x_flipped, plate_z_relative, by = two_strikes) +
          pitch_cat + two_strikes, 
        data = train_df, 
        family = binomial, 
        weights = p_full
      )
    } else {
      gmod <- gam(
        is_contact ~ s(plate_x_flipped, plate_z_relative, by = two_strikes) +
          two_strikes, 
        data = train_df, 
        family = binomial, 
        weights = p_full
      )
    }
    test_df <- test_df |>
      mutate(p_contact = predict(gmod, newdata = test_df, type = "response")) 
  }
  
  test_df <- test_df |>
    mutate(llik = p_full * (log(p_contact) * is_contact + 
                               log(1 - p_contact) * (1 - is_contact)))
  
  return(sum(test_df$llik))
}

set.seed(123)
samp_idx <- sample(1:length(statcast_split), 50)
pseudopitch_spacing = c("even", "quantile")
pitch_cat_model = c("smooth_by", "fixed_only", "ignore")
n_pseudopitches <- c(0, 50, 100, 200, 500, 1000, 2000, Inf)

all_contact_tests <- tibble()
for (grid in pseudopitch_spacing) {
  for (mod in pitch_cat_model) {
    all_contact_tests <- all_contact_tests |>
      bind_rows(list_rbind(map(
        statcast_split[samp_idx],
        \(df) list_rbind(map(
          n_pseudopitches, 
          \(n) tibble(
            nPseudopitches = n, 
            grid = grid,
            mod = mod,
            llik = test_contact_mod(df, n, grid, mod),
            ntest = ceiling(0.2 * nrow(df))
          )
        ))
      ), names_to = "PlayerId"))
  }
}
saveRDS(all_contact_tests, paste0(dat_path, "contact_tests.rds"))

contact_test_summary <- all_contact_tests |>
  group_by(nPseudopitches, grid, mod) |>
  summarize(llik = sum(llik) / sum(ntest))
# TODO: Create reactable table


fit_contact_mod <- function(player_df, n_pseudopitches = 100) {
  cat(player_df$player_name[1], "\n")
  
  prior_df <- quantile_prior_df
  prior_wt <- n_pseudopitches / (nrow(prior_df) / 2)
  prior_df <- prior_df |>
    mutate(p_full = prior_wt * wt) 
  prior_df <- prior_df |>
    crossing(two_strikes = factor(levels(player_df$two_strikes))) |>
    crossing(pitch_cat = factor(levels(player_df$pitch_cat))) |>
    mutate(pitch_cat_2s = fct_cross(pitch_cat, two_strikes))
  
  player_df <- player_df |>
    bind_rows(prior_df) |>
    filter(p_full > 0)
  
  gmod <- gam(
    is_contact ~ s(plate_x_flipped, plate_z_relative, by = two_strikes) +
      pitch_cat + two_strikes, 
    data = player_df, 
    family = binomial, 
    weights = p_full
  )
  
  return(gmod)
}

contact_models <- map(statcast_split, fit_contact_mod)
statcast_split <- map2(
  statcast_split, contact_models,
  \(df, mod) df |> 
    mutate(
      p_contact_1s = predict.gam(
        mod, 
        newdata = df |> mutate(two_strikes = 0), 
        type = "response"
      ),
      p_contact_2s = predict.gam(
        mod, 
        newdata = df |> mutate(two_strikes = 1), 
        type = "response"
      )
    )
)

saveRDS(statcast_split, paste0(dat_path, "statcast_split_w_con.rds"))


##### Contact model example -----

eg_contact_df <- quantile_prior_df |>
  mutate(p_full = 100 * wt / 625) |> 
  crossing(two_strikes = factor(levels(eg_df$two_strikes))) |>
  crossing(pitch_cat = factor(levels(eg_df$pitch_cat))) |>
  mutate(pitch_cat_2s = fct_cross(pitch_cat, two_strikes)) |>
  mutate(
    plate_x = -plate_x_flipped,
    plate_z = plate_z_relative * (eg_sz_top - eg_sz_bot) + eg_sz_bot,
    type = "Pseudopitch"
  ) |>
  bind_rows(eg_df |> mutate(type = "Actual pitch")) |>
  filter(p_full > 0)
ggplot(eg_contact_df, aes(x = plate_x, y = plate_z)) +
  geom_point(aes(size = p_full, color = type)) +
  geom_segment(data = eg_sz, aes(x = x, y = y, xend = xend, yend = yend)) +
  scale_size_area(max_size = 0.5) +
  labs(color = "Type", size = "Weight")
my_ggsave("contact_prior_example.pdf")

eg_contact <- contact_models[[eg_idx]]
pred_df <- expand.grid(
  plate_x_flipped = seq(min_x_plot, max_x_plot, by = 0.02), 
  plate_z_relative = seq(min_z_rel_plot, max_z_rel_plot, by = 0.02),
  pitch_cat = "Fastball"
) |>
  crossing(two_strikes = factor(levels(eg_df$two_strikes)))
pred_df <- pred_df |>
  mutate(
    plate_x = -plate_x_flipped,
    plate_z = plate_z_relative * (eg_sz_top - eg_sz_bot) + eg_sz_bot,
    contact_pred = predict(eg_contact, newdata = pred_df, type = "response")
  ) 
pred_df <- pred_df |>
  mutate(two_strikes = fct_recode(
    two_strikes, "usual swing" = "0", "two-strike swing" = "1"
  ))
ggplot(pred_df) +
  geom_tile(aes(x = plate_x, y = plate_z, fill = contact_pred)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = mean(statcast$is_contact, na.rm = TRUE)) +
  geom_segment(data = eg_sz, aes(x = x, y = y, xend = xend, yend = yend)) +
  labs(fill = "Prob. contact") +
  facet_wrap(~two_strikes)
my_ggsave("contact_model_example.pdf")


##### Fouls ----------

all_fouls <- gam(is_foul ~ s(plate_x_flipped, plate_z_relative, by = pitch_cat) +
                   pitch_cat,
                 data = statcast |> filter(is_contact == 1),
                 family = binomial)
# The following fit is far worse:
# all_fouls2 <- gam(is_foul ~ s(plate_x_flipped, plate_z_relative) + pitch_cat,
#                   data = statcast |> filter(is_contact == 1),
#                   family = binomial)

pred_df <- expand.grid(
  plate_x_flipped = seq(min_x, max_x, by = 0.05), 
  plate_z_relative = seq(min_z_rel, max_z_rel, by = 0.05),
  pitch_cat = factor(levels(statcast$pitch_cat))
)
pred_df <- pred_df |>
  mutate(
    plate_x = plate_x_flipped,
    plate_z = plate_z_relative * (mean_sz_top - mean_sz_bot) + mean_sz_bot,
    foul_pred = predict(all_fouls, newdata = pred_df, type = "response")
  )

p_foul = sum(statcast$is_foul) / sum(statcast$is_contact)
ggplot(pred_df) +
  geom_tile(aes(x = plate_x, y = plate_z, fill = foul_pred)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = p_foul) +
  geom_segment(data = sz, aes(x = x, y = y, xend = xend, yend = yend)) +
  facet_wrap(~pitch_cat, nrow = 3) +
  labs(x = "plate_x (RHB)", y = "plate_z", fill = "Prob. foul")
my_ggsave("all_foul.pdf")


foul_prior_df <- expand.grid(
  plate_x_flipped = quantile(
    statcast$plate_x_flipped, seq(.005, .995, length.out = n_grid_pts)
  ), 
  plate_z_relative = quantile(
    statcast$plate_z_relative, seq(.005, .995, length.out = n_grid_pts)
  ),
  pitch_cat = factor(levels(statcast$pitch_cat))
)
foul_prior_df <- foul_prior_df |>
  mutate(foul_pred = predict(
    all_fouls, newdata = foul_prior_df, type = "response"
  ))
foul_prior_df <- foul_prior_df |>
  crossing(is_foul = 0:1) |>
  mutate(wt = is_foul * foul_pred + (1 - is_foul) * (1 - foul_pred))

test_foul_mod <- function(player_df, n_pseudopitches, seed = 666) {
  cat(player_df$player_name[1], n_pseudopitches, "\n")
  
  player_df <- player_df |>
    filter(is_contact == 1)
  
  set.seed(seed)
  train_idx <- sample(1:nrow(player_df), size = floor(0.8 * nrow(player_df)))
  test_idx <- setdiff(1:nrow(player_df), train_idx)

  test_df <- player_df |>
    slice(test_idx)
  if (is.infinite(n_pseudopitches)) {
    test_df <- test_df |>
      mutate(p_foul = predict(all_fouls, newdata = test_df, type = "response")) 
  } else {
    prior_wt <- n_pseudopitches / (nrow(foul_prior_df) / 2)
    prior_df <- foul_prior_df |>
      mutate(p_full = prior_wt * wt) 
    train_df <- player_df |>
      slice(train_idx) |>
      bind_rows(prior_df) |>
      filter(p_full > 0)
    gmod <- gam(is_foul ~ s(plate_x_flipped, plate_z_relative, by = pitch_cat) +
                  pitch_cat, 
                data = train_df, 
                family = binomial, 
                weights = p_full)
    test_df <- test_df |>
      mutate(p_foul = predict(gmod, newdata = test_df, type = "response")) 
  }
  
  test_df <- test_df |>
    mutate(llik = p_full * (log(p_foul) * (1 - is_foul) + log(1 - p_foul) * is_foul))
    
  return(sum(test_df$llik))
}

set.seed(123)
samp_idx <- sample(1:length(statcast_split), 50)
n_pseudopitches <- c(0, 100, 200, 500, 1000, 2000, 4000, Inf)

all_foul_tests <- list_rbind(map(
  statcast_split[samp_idx],
  \(df) list_rbind(map(
    n_pseudopitches, 
    \(n) tibble(
      nPseudopitches = n, 
      llik = test_foul_mod(df, n),
      ntest = ceiling(0.2 * nrow(df))
    )
  ))
), names_to = "PlayerId")
saveRDS(all_foul_tests, paste0(dat_path, "foul_tests.rds"))

foul_test_summary <- all_foul_tests |>
  group_by(nPseudopitches) |>
  summarize(llik = sum(llik) / sum(ntest))
# TODO: Create reactable table


fit_foul_mod <- function(player_df, n_pseudopitches = 1000) {
  cat(player_df$player_name[1], "\n")
  
  player_df <- player_df |>
    filter(is_contact == 1)
  
  prior_df <- foul_prior_df
  prior_wt <- n_pseudopitches / (nrow(prior_df) / 2)
  prior_df <- prior_df |>
    mutate(p_full = prior_wt * wt) 
  
  player_df <- player_df |>
    bind_rows(prior_df) |>
    filter(p_full > 0)
  
  gmod <- gam(is_foul ~ s(plate_x_flipped, plate_z_relative, by = pitch_cat) +
                pitch_cat, 
              data = player_df, 
              family = binomial, 
              weights = p_full
  )

  return(gmod)
}

foul_models <- map(statcast_split, fit_foul_mod)
statcast_split <- map2(
  statcast_split, foul_models,
  \(df, mod) df |> 
    mutate(p_foul = predict.gam(mod, newdata = df, type = "response"))
)

saveRDS(statcast_split, paste0(dat_path, "statcast_split_w_foul.rds"))


##### Foul model example -----

eg_foul_df <- foul_prior_df |>
  mutate(p_full = 1000 * wt / 1875) |> 
  mutate(
    plate_x = -plate_x_flipped,
    plate_z = plate_z_relative * (eg_sz_top - eg_sz_bot) + eg_sz_bot,
    type = "Pseudopitch"
  ) |>
  bind_rows(eg_df |> 
              filter(is_contact == 1) |> 
              mutate(type = "Actual pitch")) |>
  filter(p_full > 0)
ggplot(eg_foul_df, aes(x = plate_x, y = plate_z)) +
  geom_point(aes(size = p_full, color = type)) +
  geom_segment(data = eg_sz, aes(x = x, y = y, xend = xend, yend = yend)) +
  scale_size_area(max_size = 1.5) +
  labs(color = "Type", size = "Weight") +
  facet_wrap(~pitch_cat, nrow = 1)
my_ggsave("foul_prior_example.pdf")

eg_foul <- foul_models[[eg_idx]]
pred_df <- expand.grid(
  plate_x_flipped = seq(min_x_plot, max_x_plot, by = 0.02), 
  plate_z_relative = seq(min_z_rel_plot, max_z_rel_plot, by = 0.02),
  pitch_cat = factor(levels(eg_df$pitch_cat))
)
pred_df <- pred_df |>
  mutate(
    plate_x = -plate_x_flipped,
    plate_z = plate_z_relative * (eg_sz_top - eg_sz_bot) + eg_sz_bot,
    foul_pred = predict(eg_foul, newdata = pred_df, type = "response")
  ) 
ggplot(pred_df) +
  geom_tile(aes(x = plate_x, y = plate_z, fill = foul_pred)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = p_foul) +
  geom_segment(data = eg_sz, aes(x = x, y = y, xend = xend, yend = yend)) +
  labs(fill = "Prob. foul") +
  facet_wrap(~pitch_cat, nrow = 1)
my_ggsave("foul_model_example.pdf")


##### Bat speed -----

all_bs <- gam(relative_bs ~ s(plate_x_flipped, plate_z_relative),
              data = statcast,
              weights = p_full)
# The following does not improve the fit much (and looks worse at the edges):
# all_bs2 <- gam(relative_bs ~ s(plate_x_flipped, plate_z_relative, by = pitch_cat),
#               data = statcast,
#               weights = p_full)

pred_df <- expand.grid(
  plate_x_flipped = seq(min_x, max_x, by = 0.05), 
  plate_z_relative = seq(min_z_rel, max_z_rel, by = 0.05)
)
pred_df <- pred_df |>
  mutate(
    plate_x = plate_x_flipped,
    plate_z = plate_z_relative * (mean_sz_top - mean_sz_bot) + mean_sz_bot,
    bs_pred = predict(all_bs, newdata = pred_df, type = "response")
  )

ggplot(pred_df) +
  geom_tile(aes(x = plate_x, y = plate_z, fill = bs_pred)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0) +
  geom_segment(data = sz, aes(x = x, y = y, xend = xend, yend = yend)) +
  labs(x = "plate_x (RHB)", y = "plate_z", fill = "BS (z-scored)")
my_ggsave("all_bs.pdf")


bs_prior_df <- expand.grid(
  plate_x_flipped = quantile(
    statcast$plate_x_flipped, seq(.005, .995, length.out = n_grid_pts)
  ), 
  plate_z_relative = quantile(
    statcast$plate_z_relative, seq(.005, .995, length.out = n_grid_pts)
  )
)
bs_prior_df <- bs_prior_df |>
  mutate(relative_bs = predict(
    all_bs, newdata = bs_prior_df, type = "response"
  ))

test_bs_mod <- function(player_df, 
                        n_pseudopitches, 
                        include_pitch_cat = FALSE,
                        seed = 666) {
  cat(player_df$player_name[1], n_pseudopitches, "\n")
  
  set.seed(seed)
  train_idx <- sample(1:nrow(player_df), size = floor(0.8 * nrow(player_df)))
  test_idx <- setdiff(1:nrow(player_df), train_idx)
  
  train_df <- player_df |> slice(train_idx) 
  test_df <- player_df |> slice(test_idx)
  
  if (is.infinite(n_pseudopitches)) {
    test_df <- test_df |>
      mutate(bs_pred = predict(all_bs, newdata = test_df, type = "response")) 
  } else {
    if (n_pseudopitches > 0) {
      prior_df <- bs_prior_df |>
        mutate(p_full = n_pseudopitches / nrow(prior_df)) |>
        crossing(two_strikes = factor(levels(player_df$two_strikes)))
      if (include_pitch_cat) {
        prior_df <- prior_df |>
          crossing(pitch_cat = factor(levels(player_df$pitch_cat)))
      }
      train_df <- train_df |>
        bind_rows(prior_df)
    }
    
    if (include_pitch_cat) {
      gmod <- gam(
        relative_bs ~ s(plate_x_flipped, plate_z_relative, by = two_strikes) +
          pitch_cat + two_strikes, 
        data = train_df, 
        weights = p_full
      )
    } else {
      gmod <- gam(
        relative_bs ~ s(plate_x_flipped, plate_z_relative, by = two_strikes) + 
          two_strikes, 
        data = train_df, 
        weights = p_full
      )
    } 
    
    test_df <- test_df |>
      mutate(bs_pred = predict(gmod, newdata = test_df, type = "response")) 
  }
  
  test_df <- test_df |>
    mutate(SqError = (relative_bs - bs_pred)^2)
  
  return(sum(test_df$SqError))
}

set.seed(123)
samp_idx <- sample(1:length(statcast_split), 100)
n_pseudopitches <- c(0, 50, 100, 200, 400, 1000, Inf)

all_bs_tests <- tibble()
for (incl_pitch_cat in c(TRUE, FALSE)) {
  all_bs_tests <- all_bs_tests |>
    bind_rows(list_rbind(map(
      statcast_split[samp_idx],
      \(df) list_rbind(map(
        n_pseudopitches, 
        \(n) tibble(
          nPseudopitches = n, 
          IncludePitchCat = incl_pitch_cat,
          SqError = test_bs_mod(df, n, incl_pitch_cat),
          ntest = ceiling(0.2 * nrow(df))
        )
      ))
    ), names_to = "PlayerId"))
}
saveRDS(all_bs_tests, paste0(dat_path, "bs_tests.rds"))

bs_test_summary <- all_bs_tests |>
  group_by(nPseudopitches, IncludePitchCat) |>
  summarize(RMSE = sqrt(sum(SqError) / sum(ntest)))
# TODO: Create reactable table


fit_bs_mod <- function(player_df, n_pseudopitches = 50) {
  cat(player_df$player_name[1], "\n")
  
  prior_df <- bs_prior_df 
  prior_df <- prior_df |>
    mutate(p_full = n_pseudopitches / nrow(prior_df)) 
  
  player_df <- player_df |>
    bind_rows(prior_df)
  
  gmod <- gam(
    relative_bs ~ s(plate_x_flipped, plate_z_relative, by = two_strikes) +
      pitch_cat + two_strikes, 
    data = player_df, 
    weights = p_full
  )
  
  return(gmod)
}

bs_models <- map(statcast_split, fit_bs_mod)
statcast_split <- map2(
  statcast_split, bs_models,
  \(df, mod) df |> 
    mutate(rbs_pred_1s = predict.gam(mod, 
                                     newdata = df |> mutate(two_strikes = 0),
                                     type = "response"),
           rbs_pred_2s = predict.gam(mod,
                                     newdata = df |> mutate(two_strikes = 1), 
                                     type = "response")) 
)
statcast_split <- map2(
  statcast_split, all_bs_dist,
  \(df, dist) df |>
    mutate(
      bs_pred_1s = my_qskt(
        pnorm(rbs_pred_1s), dist$df1, dist$gamma1, dist$mu1, dist$s1
      ),
      bs_pred_2s = my_qskt(
        pnorm(rbs_pred_2s), dist$df1, dist$gamma1, dist$mu1, dist$s1
      )
    )
)

saveRDS(statcast_split, paste0(dat_path, "statcast_split_w_bs.rds"))


##### Bat speed model example -----

eg_bs_df <- bs_prior_df |>
  mutate(p_full = 50 / nrow(bs_prior_df)) |> 
  mutate(
    plate_x = -plate_x_flipped,
    plate_z = plate_z_relative * (eg_sz_top - eg_sz_bot) + eg_sz_bot,
    type = "Pseudopitch"
  ) |>
  bind_rows(eg_df |> mutate(type = "Actual pitch"))
ggplot(eg_bs_df, aes(x = plate_x, y = plate_z)) +
  geom_point(aes(size = p_full, color = type)) +
  geom_segment(data = eg_sz, aes(x = x, y = y, xend = xend, yend = yend)) +
  scale_size_area(max_size = 1.5) +
  labs(color = "Type", size = "Weight") 
my_ggsave("bs_prior_example.pdf")

eg_bs <- bs_models[[eg_idx]]
pred_df <- expand.grid(
  plate_x_flipped = seq(min_x_plot, max_x_plot, by = 0.02), 
  plate_z_relative = seq(min_z_rel_plot, max_z_rel_plot, by = 0.02),
  pitch_cat = "Fastball",
  two_strikes = factor(levels(eg_df$two_strikes))
)
pred_df <- pred_df |>
  mutate(
    plate_x = -plate_x_flipped,
    plate_z = plate_z_relative * (eg_sz_top - eg_sz_bot) + eg_sz_bot,
    relative_bs = predict(eg_bs, newdata = pred_df, type = "response")
  ) 
pred_df <- pred_df |>
  mutate(two_strikes = fct_recode(
    two_strikes, "usual swing" = "0", "two-strike swing" = "1"
  ))
dist <- all_bs_dist[[eg_idx]]
pred_df <- pred_df |>
  mutate(bs_pred = my_qskt(
    pnorm(relative_bs), dist$df1, dist$gamma1, dist$mu1, dist$s1
  ))
ggplot(pred_df) +
  geom_tile(aes(x = plate_x, y = plate_z, fill = bs_pred)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = median(eg_df$bat_speed)) +
  geom_segment(data = eg_sz, aes(x = x, y = y, xend = xend, yend = yend)) +
  labs(fill = "Bat speed") +
  facet_wrap(~two_strikes, nrow = 1)
my_ggsave("bs_model_example.pdf")

pred_df2 <- pred_df |>
  select(-relative_bs) |>
  pivot_wider(names_from = two_strikes, values_from = bs_pred) |>
  mutate(bs_diff = `two-strike swing` - `usual swing`)
ggplot(pred_df2) +
  geom_tile(aes(x = plate_x, y = plate_z, fill = bs_diff)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0) +
  geom_segment(data = eg_sz, aes(x = x, y = y, xend = xend, yend = yend)) +
  labs(fill = "Bat speed diff.") 
my_ggsave("bs_diff_example.pdf")


##### E(wOBA) -----

# Release speed goes unused:
ewoba_mod <- scam(woba_value ~ s(plate_x_flipped, plate_z_relative) + 
                    s(bat_speed, bs = "mpi"), 
                  data = statcast |> filter(is_contact == 1, is_foul == 0),
                  weights = p_full)
pred_df <- tibble(
  bat_speed = seq(50, 85, by = 0.5),
  plate_x_flipped = 0,
  plate_z_relative = 0.5
)
pred_df <- pred_df |>
  mutate(ewoba = predict.scam(ewoba_mod, newdata = pred_df))
median_bs <- statcast |>
  filter(p_full > 0.5) |>
  group_by(batter) |>
  summarize(median_bs = median(bat_speed))
ggplot(pred_df) +
  geom_line(aes(x = bat_speed, y = ewoba)) +
  geom_rug(aes(x = median_bs), data = median_bs) +
  labs(x = "Bat speed (mph)", y = "Expected wOBA")
my_ggsave("ewoba.pdf")


statcast_split <- map(
  statcast_split,
  \(df) df |> 
    mutate(ewobacon_1s = predict.scam(ewoba_mod, 
                                      newdata = df |> mutate(bat_speed = bs_pred_1s), 
                                      type = "response"),
           ewobacon_2s = predict.scam(ewoba_mod,
                                      newdata = df |> mutate(bat_speed = bs_pred_2s), 
                                      type = "response"))
)

statcast_split <- map(
  statcast_split,
  \(df) df |>
    mutate(ewoba_1s = p_contact_1s * (1 - p_foul) * ewobacon_1s,
           wt_1s = 1 - p_contact_1s * p_foul,
           ewoba_2s = p_contact_2s * (1 - p_foul) * ewobacon_2s,
           wt_2s = 1 - p_contact_2s * p_foul)
)

saveRDS(statcast_split, paste0(dat_path, "statcast_split_final.rds"))


##### Summaries -----

player_summs <- map(
  statcast_split,
  \(df) df |>
    filter(strikes == 2) |>
    summarize(player_name = player_name[1],
              ewoba_1s = sum(ewoba_1s) / sum(wt_1s),
              ewoba_2s = sum(ewoba_2s) / sum(wt_2s))
) |>
  list_rbind() |>
  mutate(batter = as.numeric(names(statcast_split)))

player_woba <- readRDS(paste0(dat_path, "player_woba.rds")) |>
  select(-player_name)

player_summs <- player_summs |>
  left_join(player_woba, by = "batter") |>
  mutate(ewoba_diff = ewoba_2s - ewoba_1s)

ggplot(player_summs |> filter(n > 200), aes(x = ewoba_diff)) + 
  geom_histogram(binwidth = 0.002) +
  labs(x = "Per-PA wOBA differential",
       y = "Number of hitters")
my_ggsave("results_woba_gain.pdf")

neg_diffs <- player_summs |>
  filter(ewoba_diff < 0) |>
  pull(ewoba_diff)
sd_est <- sd(c(neg_diffs, -neg_diffs))
ebnm_res <- ebnm::ebnm_unimodal_nonnegative(player_summs$ewoba_diff, s = sd_est)
pct_increasing_woba <- 1 - ebnm_res$fitted_g$pi[1]

leader_table <- player_summs |>
  filter(n > 500) |>
  select(player_name, n, woba, prop_2s, ewoba_diff) |>
  slice_max(n = 10, ewoba_diff)
saveRDS(leader_table, paste0(dat_path, "leader_table.rds"))


##### Examples -----

do_example <- function(player_name) {
  eg_idx <- str_which(player_map, player_name)
  eg_df <- statcast_split[[eg_idx]] 
  
  eg_sz_bot <- mean(eg_df$sz_bot)
  eg_sz_top <- mean(eg_df$sz_top)
  
  eg_sz <- tibble(
    x = c(-0.71, 0.71, 0.71, -0.71),
    xend = c(0.71, 0.71, -0.71, -0.71),
    y = c(eg_sz_top, eg_sz_top, eg_sz_bot, eg_sz_bot),
    yend = c(eg_sz_top, eg_sz_bot, eg_sz_bot, eg_sz_top)
  )
  
  eg_contact <- contact_models[[eg_idx]]
  pred_df <- expand.grid(
    plate_x_flipped = seq(min_x_plot, max_x_plot, by = 0.02), 
    plate_z_relative = seq(min_z_rel_plot, max_z_rel_plot, by = 0.02),
    pitch_cat = "Fastball"
  ) |>
    crossing(two_strikes = factor(levels(eg_df$two_strikes)))
  is_LHB <- mean(eg_df$stand == "L") > 0.5
  pred_df <- pred_df |>
    mutate(
      plate_x = plate_x_flipped * ifelse(is_LHB, -1, 1),
      plate_z = plate_z_relative * (eg_sz_top - eg_sz_bot) + eg_sz_bot,
      contact_pred = predict(eg_contact, newdata = pred_df, type = "response")
    ) 
  pred_df <- pred_df |>
    mutate(two_strikes = fct_recode(
      two_strikes, "usual swing" = "0", "two-strike swing" = "1"
    ))
  p1 <- ggplot(pred_df) +
    geom_tile(aes(x = plate_x, y = plate_z, fill = contact_pred)) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = mean(statcast$is_contact, na.rm = TRUE)) +
    geom_segment(data = eg_sz, aes(x = x, y = y, xend = xend, yend = yend)) +
    labs(fill = "p. contact") +
    facet_wrap(~two_strikes)
  
  eg_bs <- bs_models[[eg_idx]]
  pred_df <- expand.grid(
    plate_x_flipped = seq(min_x_plot, max_x_plot, by = 0.02), 
    plate_z_relative = seq(min_z_rel_plot, max_z_rel_plot, by = 0.02),
    pitch_cat = "Fastball",
    two_strikes = factor(levels(eg_df$two_strikes))
  )
  pred_df <- pred_df |>
    mutate(
      plate_x = plate_x_flipped * ifelse(is_LHB, -1, 1),
      plate_z = plate_z_relative * (eg_sz_top - eg_sz_bot) + eg_sz_bot,
      relative_bs = predict(eg_bs, newdata = pred_df, type = "response")
    ) 
  pred_df <- pred_df |>
    mutate(two_strikes = fct_recode(
      two_strikes, "usual swing" = "0", "two-strike swing" = "1"
    ))
  dist <- all_bs_dist[[eg_idx]]
  pred_df <- pred_df |>
    mutate(bs_pred = my_qskt(
      pnorm(relative_bs), dist$df1, dist$gamma1, dist$mu1, dist$s1
    ))
  pred_df <- pred_df |>
    select(-relative_bs) |>
    pivot_wider(names_from = two_strikes, values_from = bs_pred) |>
    mutate(bs_diff = `two-strike swing` - `usual swing`)
  p2 <- ggplot(pred_df) +
    geom_tile(aes(x = plate_x, y = plate_z, fill = bs_diff)) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = 0) +
    geom_segment(data = eg_sz, aes(x = x, y = y, xend = xend, yend = yend)) +
    labs(fill = "BS diff.") 
  
  return(cowplot::plot_grid(p1, p2, nrow = 1, rel_widths = c(.6, .4)))
}

do_example("Olson")
my_ggsave("example_olson.pdf", width = 8, height = 6)

do_example("Rafaela")
my_ggsave("example_rafaela.pdf", width = 8, height = 6)

do_example("McCutchen")
my_ggsave("example_mccutchen.pdf", width = 8, height = 6)

do_example("Bellinger")
my_ggsave("example_bellinger.pdf", width = 8, height = 6)

do_example("Gim")
my_ggsave("example_gimenez.pdf", width = 8, height = 6)

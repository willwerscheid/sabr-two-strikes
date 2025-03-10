library(tidyverse) 
theme_set(theme_minimal())

library(mgcv)
library(scam)

# Load and prep data: -----

min_swings <- 200

# Includes both regular season and postseason:
statcast <- readRDS("Data/statcast2024.rds")

bunt_descriptions <- c("bunt_foul_tip",
                       "foul_bunt",
                       "missed_bunt")
called_pitch_descriptions <- c("ball",
                               "blocked_ball",
                               "called_strike")
swing_descriptions <- c("foul",
                        "foul_tip",
                        "hit_into_play", 
                        "swinging_strike", 
                        "swinging_strike_blocked")
contact_descriptions <- c("foul", "hit_into_play")

statcast <- statcast |>
  select(
    game_date, game_pk, game_type, 
    at_bat_number, batter, player_name, stand, p_throws, 
    pitch_number, balls, strikes, 
    pitch_type, release_speed, release_extension,
    plate_x, plate_z, sz_top, sz_bot,
    type, events, description, woba_value, woba_denom,
    launch_speed, launch_angle, bat_speed, swing_length
  ) |>
  mutate(across(
    c(game_pk, at_bat_number, batter, pitch_number, balls),
    as.integer
  )) |>
  mutate(
    balls = ifelse(balls > 3, 3, balls),
    strikes = ifelse(strikes > 2, 2, strikes)
  ) |>
  mutate(
    the_count = factor(paste(balls, strikes, sep = "-")),
    .after = strikes
  ) |>
  mutate(
    stand = as.factor(stand),
    p_throws = as.factor(p_throws),
  ) |>
  mutate(
    pitch_cat = factor(case_when(
      pitch_type %in% c("FF", "SI", "FC") ~ "Fastball",
      pitch_type %in% c("CH", "FS", "FO", "SC") ~ "Offspeed",
      pitch_type %in% c("", "EP", "FA", "IN", "PO") ~ "Other",
      TRUE ~ "Breaking"
    )), 
    .before = pitch_type
  ) |>
  mutate(
    plate_x_flipped = ifelse(stand == "L", -plate_x, plate_x),
    plate_z_relative = (plate_z - sz_bot) / (sz_top - sz_bot),
    .after = plate_z
  ) |>
  mutate(
    la_diff_from_optimal = pmin(abs(launch_angle - 20), 30), 
    .after = launch_angle
  ) |>
  mutate(
    two_strikes = factor(1L * (strikes == 2)),
    is_strike = 1L * (type == "S"),
    is_swing = 1L * (description %in% swing_descriptions) * (events != "sac_bunt"),
    is_contact = 1L * (description %in% contact_descriptions) * (events != "sac_bunt"),
    is_foul = 1L * (description == "foul"),
    is_hbp = 1L * (description == "hit_by_pitch"),
    is_called = 1L * (description %in% called_pitch_descriptions)
  ) |>
  arrange(game_pk, at_bat_number, pitch_number)

statcast_split <- statcast |>
  filter(is_swing == 1, !is.na(bat_speed), pitch_cat != "Other") |>
  mutate(pitch_cat = fct_drop(pitch_cat)) |>
  group_by(player_name) |>
  filter(n() >= min_swings) |>
  ungroup()
statcast_split <- statcast_split |>
  split(statcast_split$batter)

player_map <- map_chr(statcast_split, list("player_name", 1))


# wOBA as function of location and bat speed:

# gam is from mcgv package:
ewoba_mod <- gam(woba_value ~ s(plate_x_flipped, plate_z) + 
                   s(bat_speed), 
                 data = statcast |> filter(is_contact == 1, is_foul == 0))
plot(ewoba_mod)

# overfitting; examples of 40 mph swings:
slow <- statcast |> filter(is_contact == 1, is_foul == 0, round(bat_speed) == 40)
# https://www.mlb.com/video/george-springer-grounds-out-softly-first-baseman-triston-casas-to-pitcher

fast <- statcast |> filter(is_contact == 1, is_foul == 0, round(bat_speed) == 88)
# https://www.mlb.com/video/aaron-judge-lines-out-to-right-fielder-jo-adell

# scam is from scam package:
ewoba_mod <- scam(woba_value ~ s(plate_x_flipped, plate_z) + 
                    s(bat_speed, bs = "mpi"), 
                  data = statcast |> filter(is_contact == 1, is_foul == 0))
plot(ewoba_mod)

eg_idx <- str_which(player_map, "Yordan")
eg_df <- statcast_split[[eg_idx]]

eg_contact <- gam(is_contact ~ s(plate_x, plate_z, by = two_strikes) + pitch_cat, 
                  data = eg_df, 
                  family = binomial)
eg_bs <- gam(bat_speed ~ s(plate_x, plate_z, by = two_strikes), 
             data = eg_df)

sz <- tibble(
  x = c(-0.71, 0.71, 0.71, -0.71),
  xend = c(0.71, 0.71, -0.71, -0.71),
  y = c(mean(eg_df$sz_top), mean(eg_df$sz_top), mean(eg_df$sz_bot), mean(eg_df$sz_bot)),
  yend = c(mean(eg_df$sz_top), mean(eg_df$sz_bot), mean(eg_df$sz_bot), mean(eg_df$sz_top))
)

pred_df <- expand.grid(
  plate_x = seq(-1.5, 1.5, by = 0.02), 
  plate_z = seq(1, 4.25, by = 0.02),
  two_strikes = 0:1,
  pitch_cat = "Fastball"
) 
pred_df <- pred_df |>
  mutate(contact_pred = predict.gam(eg_contact, newdata = pred_df, type = "response")) |>
  mutate(bs_pred = predict.gam(eg_bs, newdata = pred_df, type = "response"))
ggplot(pred_df)  +
  geom_tile(aes(x = plate_x, y = plate_z, fill = contact_pred)) +
  geom_segment(data = sz, aes(x = x, y = y, xend = xend, yend = yend)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = median(pred_df$contact_pred)) +
  facet_wrap(~two_strikes, 
             labeller = as_labeller(c(`0` = "Usual swing", `1` = "Two-strike swing"))) +
  labs(fill = "Prob. contact")
ggplot(pred_df)  +
  geom_tile(aes(x = plate_x, y = plate_z, fill = bs_pred)) +
  geom_segment(data = sz, aes(x = x, y = y, xend = xend, yend = yend)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = median(pred_df$bs_pred)) +
  facet_wrap(~two_strikes, 
             labeller = as_labeller(c(`0` = "Usual swing", `1` = "Two-strike swing"))) +
  labs(fill = "Bat speed")

pred_df <- pred_df |>
  pivot_wider(names_from = two_strikes, 
              values_from = c(contact_pred, bs_pred))
pred_df <- pred_df |>
  mutate(contact_diff = contact_pred_1 - contact_pred_0,
         bs_diff = bs_pred_1 - bs_pred_0)
ggplot(pred_df, aes(x = plate_x, y = plate_z)) + 
  lims(x = range(pred_df$plate_x), y = range(pred_df$plate_z)) +
  geom_tile(aes(fill = contact_diff)) +
  geom_segment(data = sz, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_point(data = eg_df |> filter(strikes == 2), aes(x = plate_x, y = plate_z), size = 0.5) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(fill = "Contact diff.")
ggplot(pred_df, aes(x = plate_x, y = plate_z)) + 
  lims(x = range(pred_df$plate_x), y = range(pred_df$plate_z)) +
  geom_tile(aes(fill = bs_diff)) +
  geom_segment(data = sz, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_point(data = eg_df |> filter(strikes == 2), aes(x = plate_x, y = plate_z), size = 0.5) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(fill = "BS diff.")

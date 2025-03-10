library(tidyverse) 
theme_set(theme_minimal())
library(mgcv)
library(scam)
source("Code/est_bs_dist.R")

fig_path <- "Writing/images/"
dat_path <- "Writing/data/"
my_ggsave <- function(fname) {
  ggsave(paste0(fig_path, fname), width = 4, height = 3)
}

min_swings <- 200

# Load and prep data: -----

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


# Part 0. How much does the two-strike approach matter? -----

player_summs <- statcast |>
  group_by(batter, player_name, game_pk, game_date, at_bat_number) |>
  filter(pitch_number == max(pitch_number)) |>
  group_by(batter, player_name) |>
  summarize(n = n(), 
            prop_2s = mean(strikes == 2 & is_swing == 1), 
            woba = mean(woba_value, na.rm = TRUE)) |>
  ungroup()
saveRDS(player_summs, paste0(dat_path, "player_woba.rds"))

# For most hitters, between 35% and 50% of all PAs end with a two-strike swing:
ggplot(player_summs |> filter(n >= 100), aes(x = prop_2s)) + 
  geom_histogram(binwidth = 0.01) +
  labs(x = "Proportion of PAs ending in two-strike swings",
       y = "Number of hitters")
my_ggsave("prop2s_hist.pdf")
ggplot(player_summs |> filter(n > 100), aes(x = prop_2s, y = woba)) + 
  geom_point(size = 0.5) + 
  geom_smooth(method = "lm") +
  labs(x = "Proportion of PAs ending in two-strike swings",
       y = "wOBA") +
  scale_y_continuous(labels = scales::label_number(0.001))
my_ggsave("prop2s_woba.pdf")


# Part I. Illustrate contact/power tradeoff. -----

### Part IA. Contact. -----

contact_models <- map(
  statcast_split,
  \(df) gam(is_contact ~ s(plate_x, plate_z, k = 4) + two_strikes + pitch_cat, 
            data = df, 
            family = binomial)
)

mod_summs <- tibble(
  Batter = names(statcast_split),
  Contact = map(contact_models, summary) |> map_dbl(list("p.coeff", 2)),
  Contact_p = map(contact_models, summary) |> map_dbl(list("p.pv", 2))
)

ggplot(mod_summs, aes(x = Contact)) + 
  geom_histogram(binwidth = 0.1, boundary = 0) +
  labs(x = "Two-strike coefficient estimate", y = "Number of hitters")
my_ggsave("contact_mod_coef.pdf")
ggplot(mod_summs, aes(x = Contact_p)) + 
  geom_histogram(binwidth = 0.05, boundary = 0) +
  labs(x = "p-value", y = "Number of hitters")
my_ggsave("contact_mod_p.pdf")

contact_pi0 <- qvalue::pi0est(mod_summs$Contact_p)$pi0

ggplot(mod_summs, aes(x = Contact_p)) + 
  geom_histogram(binwidth = 0.05, boundary = 0) +
  labs(x = "p-value", y = "Number of hitters") +
  geom_hline(yintercept = nrow(mod_summs) * contact_pi0 * 0.05,
             linetype = "dashed", color = "red")
my_ggsave("contact_mod_q.pdf")


### Part IB. Power. -----

##### Exit velocity approach: -----

# Why not EV?
#.  1. Don't have data for all swings; we lose statistical power.
#.  2. Bat speed is under the immediate control of the hitter; exit velocity depends on 
#.       external variables like release speed and release extension, as well as 
#.       variables that do not correlate well with two-strike counts ("quality of
#.       contact" variables like launch angle).

ev_models <- map(
  statcast_split,
  \(df) gam(launch_speed ~ s(plate_x, plate_z, k = 4) + two_strikes + pitch_cat, 
            data = df)
)
mod_summs <- mod_summs |> 
  mutate(
    EV_p = map(ev_models, summary) |> map_dbl(list("p.pv", 2))
  ) 

# Lacks statistical power:
ggplot(mod_summs, aes(x = EV_p)) + 
  geom_histogram(binwidth = 0.1, boundary = 0) +
  labs(x = "p-value", y = "Number of hitters")
my_ggsave("ev_mod_p.pdf")
# Estimated proportion of true nulls:
ev_pi0 <- qvalue::pi0est(mod_summs$EV_p)$pi0
ggplot(mod_summs, aes(x = EV_p)) + 
  geom_histogram(binwidth = 0.1, boundary = 0) +
  labs(x = "p-value", y = "Number of hitters") +
  geom_hline(yintercept = nrow(mod_summs) * ev_pi0 * 0.1,
             linetype = "dashed", color = "red")
my_ggsave("ev_mod_q.pdf")

rm(ev_models)


##### Bat speed, naive approach -----

bs_models_unweighted <- map(
  statcast_split,
  \(df) gam(bat_speed ~ s(plate_x, plate_z, k = 4) + two_strikes, 
            data = df)
)

mod_summs <- mod_summs |> 
  mutate(
    BSnowts = map(bs_models_unweighted, summary) |> 
      map_dbl(list("p.coeff", 2)),
    BSnowts_p = map(bs_models_unweighted, summary) |> 
      map_dbl(list("p.pv", 2)),
  ) 

ggplot(mod_summs, aes(x = BSnowts_p)) + 
  geom_histogram(binwidth = 0.05, boundary = 0) +
  labs(x = "p-value", y = "Number of hitters")
my_ggsave("bs_mod_naive_p.pdf")
bs_pi0 <- qvalue::pi0est(mod_summs$BSnowts_p)$pi0
ggplot(mod_summs, aes(x = BSnowts_p)) + 
  geom_histogram(binwidth = 0.05, boundary = 0) +
  labs(x = "p-value", y = "Number of hitters") +
  geom_hline(yintercept = nrow(mod_summs) * bs_pi0 * 0.05,
             linetype = "dashed", color = "red")
my_ggsave("bs_mod_naive_q.pdf")

rm(bs_models_unweighted)


##### What is a swing? -----

swings_bs <- statcast |>
  filter(is_swing == 1, !is.na(bat_speed))

# Assume that the overall distribution is a mixture of two skew-t distns ("weak"
#   and "full" swings) and perform a deconvolution. First fit a league-wide mixture:
lg_bs_dist <- est_bs_dist(swings_bs$bat_speed)
swings_bs <- swings_bs |>
  mutate(p_full = lg_bs_dist$lik_mat$p_full)

# The fit is excellent:
ggplot(swings_bs) +
  xlim(c(0, 100)) +
  geom_histogram(aes(x = bat_speed, after_stat(density)), binwidth = 1, fill = "light grey") +
  geom_function(fun = \(x) do.call(my_dmix, c(list(x), lg_bs_dist[1:9])), color = "black") +
  labs(x = "Bat speed (mph)", y = "Density")
my_ggsave("lg_bs_dist.pdf")

# But it is not appropriate for players who swing very hard:
swings_bs |> 
  filter(player_name == "Stanton, Giancarlo") |>
  group_by(bat_speed = round(bat_speed)) |>
  summarize(n = n(), p_full = mean(p_full)) |>
  ggplot() +
  geom_col(aes(x = bat_speed, y = n, fill = p_full)) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0.5) +
  labs(x = "Bat speed (mph)", y = "Number of swings", fill = "Prob. full swing",
       title = "Deconvolution using league-wide dist.")
my_ggsave("stanton_lg_deconv.pdf")

rm(swings_bs)
rm(lg_bs_dist)

##### Player-wise deconvolution of bat speeds: -----

# Remove very obvious non-swings:
statcast_split <- map(
  statcast_split,
  \(df) df |> filter(!is.na(bat_speed), bat_speed > 30)
)

# Fit a mixture to each player:
all_bs_dist <- map(statcast_split, \(df) est_bs_dist(df$bat_speed))
saveRDS(all_bs_dist, paste0(dat_path, "all_bs_dist.rds"))

# all_bs_dist_par <- map(1:9, \(i) map_dbl(all_bs_dist, i))
# all_bs_dist_par <- Reduce(cbind, all_bs_dist_par)
# colnames(all_bs_dist_par) <- names(all_bs_dist[[1]][1:9])
# all_bs_dist_par <- as_tibble(all_bs_dist_par) |>
#   pivot_longer(everything(), names_to = "Parameter", values_to = "Value")
# ggplot(all_bs_dist_par) +
#   geom_histogram(aes(x = Value), bins = 20) +
#   facet_wrap(~Parameter, scales = "free")

statcast_split <- map2(
  statcast_split, all_bs_dist, 
  \(x, y) x |> 
    left_join(y |> pluck("lik_mat") |> select(bat_speed, p_full) |> distinct(), 
              by = "bat_speed")
)

# Stanton's histogram looks much better:
statcast_split[[str_which(player_map, "Stanton")]] |> 
  group_by(bat_speed = round(bat_speed)) |>
  summarize(n = n(), p_full = mean(p_full)) |>
  ggplot() +
  geom_col(aes(x = bat_speed, y = n, fill = p_full)) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0.5) +
  labs(x = "Bat speed (mph)", y = "Number of swings", fill = "Prob. full swing")
       #title = "Deconvolution using player-specific dist.")
my_ggsave("stanton_pl_deconv.pdf")

# Look at overall trends; some points look problematic:
ggplot(list_rbind(statcast_split), aes(bat_speed, p_full)) +
  geom_bin_2d(bins = 100) + geom_smooth(method = "glm", method.args = list(family = binomial)) +
  scale_fill_gradient(trans = "log10") +
  labs(x = "Bat speed (mph)", y = "Probability of full swing")
my_ggsave("all_pl_deconv.pdf")

# These are often indicative of outlying bat speeds:
statcast_split[[str_which(player_map, "Vidal")]] |> 
  group_by(bat_speed = round(bat_speed)) |>
  summarize(n = n(), p_full = mean(p_full)) |>
  ggplot() +
  geom_col(aes(x = bat_speed, y = n, fill = p_full)) +
  scale_fill_gradient(low = "red", high = "blue") +
  labs(x = "Bat speed (mph)", y = "Number of swings", fill = "Prob. full swing")
# Swinging strike on 9/22/24
my_ggsave("vidal_pl_deconv.pdf")

saveRDS(statcast_split, paste0(dat_path, "statcast_split.rds"))
rm(all_bs_dist)


##### Bat speed, improved approach -----

bs_models <- map(
  statcast_split,
  \(df) gam(bat_speed ~ s(plate_x, plate_z, k = 4) + two_strikes, 
            data = df,
            weights = p_full)
)

mod_summs <- mod_summs |> 
  mutate(
    BSwts = map(bs_models, summary) |> map_dbl(list("p.coeff", 2)),
    BSwts_p = map(bs_models, summary) |> map_dbl(list("p.pv", 2)),
  ) 

ggplot(mod_summs, aes(x = BSwts_p)) + 
  geom_histogram(binwidth = 0.05, boundary = 0) +
  labs(x = "p-value", y = "Number of hitters")
my_ggsave("bs_mod_p.pdf")
qvalue::pi0est(mod_summs$BSwts_p)$pi0

ggplot(mod_summs, aes(BSwts)) +
  geom_histogram(binwidth = 0.2) +
  labs(x = "Two-strike coefficient estimate", y = "Number of hitters")
my_ggsave("bs_mod_coef.pdf")

ggplot(mod_summs, aes(BSwts, BSnowts)) + 
  geom_point() +
  geom_abline(linetype = "dashed") +
  labs(y = "Est. BS differential using all swings",
       x = "Est. BS differential using full swings only")


### Final results -----

mod_summs <- mod_summs |>
  mutate(BSwts_padj = qvalue::qvalue(BSwts_p, fdr.level = 0.05)$qvalues,
         Contact_padj = qvalue::qvalue(Contact_p, fdr.level = 0.05)$qvalues,
         BSSignif = (BSwts_padj < 0.05),
         ContactSignif = (Contact_padj < 0.05)) |>
  mutate(Signif = case_when(BSSignif & ContactSignif ~ "Both",
                            BSSignif ~ "BS only",
                            ContactSignif ~ "Contact only",
                            TRUE ~ "Neither"))
ggplot(mod_summs, aes(BSwts, Contact)) +
  geom_point(aes(color = Signif)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_smooth(method = "lm", linetype = "dashed", color = "black", se = FALSE) + 
  labs(x = "Two-strike bat speed diff. (mph)",
       y = "Two-strike contact rate diff. (log odds)",
       color = "Significant") +
  scale_color_manual(values = c("red", "orange", "violet", "blue"))
ggsave(paste0(fig_path, "all_coefs.pdf"), width = 6, height = 3)  

rm(bs_models)
rm(contact_models)
gc()


# Part II. Quantify yield of tradeoff. -----

# Since we are no longer interesting in estimating a lone coefficient, we can
#.  instead model contact "maps." There is no reason to assume that a two-strike
#.  swing will uniformly increase contact rates across all zones.

contact_models <- map(
  statcast_split,
  \(df) gam(is_contact ~ s(plate_x, plate_z, by = two_strikes) + pitch_cat, 
            data = df, 
            family = binomial, 
            weights = p_full)
)
statcast_split <- map2(
  statcast_split, contact_models,
  \(df, mod) df |> 
    mutate(p_contact_1s = predict.gam(mod, 
                                      newdata = df |> mutate(two_strikes = 0), 
                                      type = "response"),
           p_contact_2s = predict.gam(mod,
                                      newdata = df |> mutate(two_strikes = 1), 
                                      type = "response"))
)

foul_models <- map(
  statcast_split,
  \(df) gam(is_foul ~ s(plate_x, plate_z) + pitch_cat, 
            data = df |> filter(is_contact == 1), 
            family = binomial, 
            weights = p_full)
)
statcast_split <- map2(
  statcast_split, foul_models,
  \(df, mod) df |> 
    mutate(p_foul = predict.gam(mod, newdata = df, type = "response"))
)

bs_models <- map(
  statcast_split,
  \(df) gam(bat_speed ~ s(plate_x, plate_z, by = two_strikes) + pitch_cat, 
            data = df, 
            weights = p_full)
)
statcast_split <- map2(
  statcast_split, bs_models,
  \(df, mod) df |> 
    mutate(bs_1s = predict.gam(mod, 
                               newdata = df |> mutate(two_strikes = 0), 
                               type = "response"),
           bs_2s = predict.gam(mod,
                               newdata = df |> mutate(two_strikes = 1), 
                               type = "response"))
)

# Release speed goes unused:
ewoba_mod <- scam(woba_value ~ s(plate_x_flipped, plate_z) + 
                    s(bat_speed, bs = "mpi"), 
                  data = statcast |> filter(is_contact == 1, is_foul == 0))
pred_df <- tibble(
  bat_speed = seq(0, 90, by = 0.5),
  plate_x_flipped = 0,
  plate_z = mean(statcast$sz_top, na.rm = TRUE) / 2 + 
    mean(statcast$sz_bot, na.rm = TRUE) / 2
)
pred_df <- pred_df |>
  mutate(preds = predict.scam(ewoba_mod, newdata = pred_df))
ggplot(pred_df, aes(x = bat_speed, y = preds)) +
  geom_line() +
  labs(x = "Bat speed (mph)", y = "Expected wOBA")
my_ggsave("bs_effect.pdf")

statcast_split <- map(
  statcast_split,
  \(df) df |> 
    mutate(ewobacon_1s = predict.scam(ewoba_mod, 
                                      newdata = df |> mutate(bat_speed = bs_1s), 
                                      type = "response"),
           ewobacon_2s = predict.scam(ewoba_mod,
                                      newdata = df |> mutate(bat_speed = bs_2s), 
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

statcast_split_summ <- map(
  statcast_split,
  \(df) df |>
    filter(strikes == 2) |>
    summarize(ewoba_1s = sum(ewoba_1s) / sum(wt_1s),
              ewoba_2s = sum(ewoba_2s) / sum(wt_2s))
) |>
  list_rbind() |>
  mutate(batter = as.numeric(names(statcast_split)))

player_summs <- player_summs |>
  left_join(statcast_split_summ, by = "batter")
player_summs <- player_summs |>
  mutate(ewoba_diff = ewoba_2s - ewoba_1s,
         diff_2s = prop_2s * ewoba_diff)

ggplot(player_summs |> filter(n > 200), aes(x = ewoba_diff)) + 
  geom_histogram(binwidth = 0.002) +
  labs(x = "Per-PA wOBA differential",
       y = "Number of hitters")
my_ggsave("results_woba_gain.pdf")
ggplot(player_summs |> filter(n > 200), aes(x = diff_2s)) + 
  geom_histogram(binwidth = 0.001) +
  labs(x = "Share of overall wOBA",
       y = "Number of hitters")
my_ggsave("results_woba_share.pdf")
ggplot(player_summs |> filter(n > 200), aes(x = prop_2s, y = ewoba_diff)) + 
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "Proportion of PAs ending in two-strike swings",
       y = "Two-strike wOBA differential due to swing adjustments")

player_summs |>
  filter(n > 200) |>
  summarize(ewoba_diff = mean(ewoba_diff), diff_2s = mean(diff_2s))

leader_table <- player_summs |>
  filter(n > 500) |>
  select(player_name, n, woba, ewoba_diff, diff_2s) |>
  slice_max(n = 10, ewoba_diff)
saveRDS(leader_table, paste0(dat_path, "leader_table.rds"))


# Part III. Analyze "adjustments." -----

eg_idx <- str_which(player_map, "Yordan")
mod_summs[eg_idx, ]

eg_df <- statcast_split[[eg_idx]] 
eg_contact <- contact_models[[eg_idx]]
summary(eg_contact)
eg_bs <- bs_models[[eg_idx]]
summary(eg_bs)
eg_foul <- foul_models[[eg_idx]]

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
my_ggsave("yordan_contact.pdf")
ggplot(pred_df)  +
  geom_tile(aes(x = plate_x, y = plate_z, fill = bs_pred)) +
  geom_segment(data = sz, aes(x = x, y = y, xend = xend, yend = yend)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = median(pred_df$bs_pred)) +
  facet_wrap(~two_strikes, 
             labeller = as_labeller(c(`0` = "Usual swing", `1` = "Two-strike swing"))) +
  labs(fill = "Bat speed")
my_ggsave("yordan_batspeed.pdf")

pred_df <- pred_df |>
  pivot_wider(names_from = two_strikes, 
              values_from = c(contact_pred, bs_pred))
pred_df <- pred_df |>
  mutate(contact_diff = contact_pred_1 - contact_pred_0,
         bs_diff = bs_pred_1 - bs_pred_0)
pred_df <- pred_df |> 
  mutate(foul_pred = predict.gam(eg_foul, newdata = pred_df, type = "response"))
ggplot(pred_df, aes(x = plate_x, y = plate_z)) + 
  lims(x = range(pred_df$plate_x), y = range(pred_df$plate_z)) +
  geom_tile(aes(fill = contact_diff)) +
  geom_segment(data = sz, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_point(data = eg_df |> filter(strikes == 2), aes(x = plate_x, y = plate_z), size = 0.5) +
  # geom_contour(aes(z = foul_pred), binwidth = 0.1, color = "white") +
  # metR::geom_text_contour(aes(z = foul_pred), size = 2.5) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(fill = "Contact diff.")
my_ggsave("yordan_contact_diff.pdf")
ggplot(pred_df, aes(x = plate_x, y = plate_z)) + 
  lims(x = range(pred_df$plate_x), y = range(pred_df$plate_z)) +
  geom_tile(aes(fill = bs_diff)) +
  geom_segment(data = sz, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_point(data = eg_df |> filter(strikes == 2), aes(x = plate_x, y = plate_z), size = 0.5) +
  # geom_contour(aes(z = foul_pred), binwidth = 0.1, color = "white") +
  # metR::geom_text_contour(aes(z = foul_pred), size = 2.5) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(fill = "BS diff.")
my_ggsave("yordan_bs_diff.pdf")

swing_len_2s_mod <- gam(swing_length ~ s(plate_x, plate_z, k = 4) + two_strikes, 
                           data = eg_df)
summary(swing_len_2s_mod)

bs_mod_no2s <- gam(bat_speed ~ s(plate_x, plate_z) + pitch_cat, 
                   data = eg_df, 
                   weights = p_full)
bs_swing_len_mod <- gam(bs_resid ~ s(plate_x, plate_z, k = 4) + swing_length, 
                        data = eg_df |> mutate(bs_resid = resid(bs_mod_no2s)),
                        weights = p_full)
summary(bs_swing_len_mod)

# Alvarez
focus <- eg_df |> 
  filter(between(plate_x, -1, -0.5), between(plate_z, 3.25, 3.75))
# Donovan
focus <- eg_df |> 
  filter(between(plate_x, -0.75, -0.25), between(plate_z, 1.25, 1.75))
focus <- focus |> 
  filter(pitch_cat == names(which.max(table(pitch_cat)))) |>
  arrange(two_strikes, desc(ewoba_2s - ewoba_1s))

# 1-0: 8/4/24, Manuel Rodriguez swinging strike to Yordan Alvarez
# 2-2: 7/30/24, Yordan Alvarez singles on a line drive

# 1-1: 7/10/24, James McArthur swinging strike to Brendan Donovan
# 1-2: 7/10/24, James McArthur foul to Brendan Donovan


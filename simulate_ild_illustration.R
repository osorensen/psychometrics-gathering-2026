# simulate_ild_illustration.R --------------------------------------------------
#
# Illustration of intensive longitudinal data.
#
# Three latent traits (pain, stress, sleep) follow a VAR(1) process over 100
# timepoints, and each trait is measured by three items through a factor
# measurement model. Latent traits are shown in the left column and the nine
# observed items in the right column. Colour identifies the latent trait; within
# the item column, the three shades of a trait's colour identify its three
# indicators.
#
# Two figures are written:
#
#   ild_illustration.png    a single participant i
#   ild_illustration_n.png  N = 20 participants overlaid, with participant
#                           i = 1, the one of the first figure, highlighted
#
# Run from the project root.
# ------------------------------------------------------------------------------

library(ggplot2)
library(patchwork)

# Simulation settings ----------------------------------------------------------

n_timepoints <- 100
n_participants <- 20
trait_names <- c("Pain", "Stress", "Sleep")
n_traits <- length(trait_names)
n_items <- 3L # items per latent trait

# Population VAR(1) transition matrix: moderate autoregression, some cross-lagged
# effects.
Phi <- matrix(
  c(
    0.55, 0.25, -0.20,
    0.10, 0.60, -0.25,
    -0.05, -0.15, 0.50
  ),
  nrow = n_traits, byrow = TRUE,
  dimnames = list(trait_names, trait_names)
)

# Innovation covariance of the latent traits.
Psi <- matrix(
  c(
    1.00, 0.30, -0.20,
    0.30, 1.00, -0.25,
    -0.20, -0.25, 1.00
  ),
  nrow = n_traits, byrow = TRUE,
  dimnames = list(trait_names, trait_names)
) * 0.5^2

# Measurement model: loadings, intercepts and residual standard deviations for
# the three items of each trait. Items are on a 1-7 scale centred at 4.
lambda <- c(1.00, 0.85, 0.70)
nu <- 4
sigma <- 0.35

# Between-participant variation, used for the second figure.
sd_trait_mean <- 0.55 # SD of participant-specific trait levels
sd_autoregression <- 0.10 # SD of participant-specific autoregressions
sd_log_innovation <- 0.25 # SD of log innovation scale

item_names <- paste(rep(trait_names, each = n_items), rep(seq_len(n_items), n_traits))

# Simulation -------------------------------------------------------------------

# Simulate one participant: a VAR(1) latent process and its indicators. Draws
# from N(0, Psi) are made with the Cholesky factor, so that no extra package is
# needed.
simulate_participant <- function(Phi_i = Phi, Psi_i = Psi, mu_i = rep(0, n_traits)) {
  innovations <- matrix(rnorm(n_timepoints * n_traits), nrow = n_timepoints) %*%
    chol(Psi_i)

  eta <- matrix(NA_real_, nrow = n_timepoints, ncol = n_traits,
                dimnames = list(NULL, trait_names))
  eta[1, ] <- innovations[1, ]
  for (t in seq(2, n_timepoints)) {
    eta[t, ] <- Phi_i %*% eta[t - 1, ] + innovations[t, ]
  }
  # Participant-specific trait level, zero for the participant of figure 1.
  eta <- eta + rep(mu_i, each = n_timepoints)

  y <- matrix(NA_real_, nrow = n_timepoints, ncol = n_traits * n_items,
              dimnames = list(NULL, item_names))
  for (j in seq_len(n_traits)) {
    for (k in seq_len(n_items)) {
      y[, (j - 1) * n_items + k] <-
        nu + lambda[k] * eta[, j] + rnorm(n_timepoints, sd = sigma)
    }
  }

  list(eta = eta, y = y)
}

# Long data frames for one participant.
latent_frame <- function(sim, participant = 1L) {
  data.frame(
    participant = participant,
    time = rep(seq_len(n_timepoints), times = n_traits),
    trait = factor(rep(trait_names, each = n_timepoints), levels = trait_names),
    value = as.vector(sim$eta)
  )
}

item_frame <- function(sim, participant = 1L) {
  data.frame(
    participant = participant,
    time = rep(seq_len(n_timepoints), times = n_traits * n_items),
    trait = factor(rep(rep(trait_names, each = n_items), each = n_timepoints),
                   levels = trait_names),
    item = factor(rep(item_names, each = n_timepoints), levels = item_names),
    value = as.vector(sim$y)
  )
}

set.seed(2026)

# Participant 1 uses the population parameters and has no deviation from the
# population trait levels. It is simulated first, so that it is the same
# participant in both figures.
participants <- vector("list", n_participants)
participants[[1]] <- simulate_participant()

for (i in seq(2, n_participants)) {
  # Participant-specific autoregressions, innovation scale and trait levels.
  Phi_i <- Phi
  diag(Phi_i) <- diag(Phi) + rnorm(n_traits, sd = sd_autoregression)
  # Keep every participant's process stationary by shrinking the transition
  # matrix whenever its spectral radius gets too close to one.
  radius <- max(Mod(eigen(Phi_i, only.values = TRUE)$values))
  if (radius > 0.9) Phi_i <- Phi_i * (0.9 / radius)
  Psi_i <- Psi * exp(rnorm(1, sd = sd_log_innovation))^2
  mu_i <- rnorm(n_traits, sd = sd_trait_mean)

  participants[[i]] <- simulate_participant(Phi_i, Psi_i, mu_i)
}

latent_one <- latent_frame(participants[[1]])
item_one <- item_frame(participants[[1]])

latent_all <- do.call(rbind, Map(latent_frame, participants, seq_len(n_participants)))
item_all <- do.call(rbind, Map(item_frame, participants, seq_len(n_participants)))

# Colours ----------------------------------------------------------------------

# Colourblind-friendly base colour per latent trait (Okabe-Ito).
trait_colours <- c(Pain = "#0072B2", Stress = "#D55E00", Sleep = "#009E73")

mix_colour <- function(colour, target, amount) {
  mixed <- (1 - amount) * col2rgb(colour) + amount * col2rgb(target)
  rgb(mixed[1], mixed[2], mixed[3], maxColorValue = 255)
}

# Three shades per trait, so that each item is identifiable but clearly belongs
# to its latent trait.
item_colours <- unlist(lapply(trait_colours, function(colour) {
  c(
    mix_colour(colour, "black", 0.30),
    colour,
    mix_colour(colour, "white", 0.40)
  )
}))
names(item_colours) <- item_names

# Shared theme -----------------------------------------------------------------

theme_ild <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey92"),
    axis.line.x = element_line(colour = "grey70"),
    axis.ticks.x = element_line(colour = "grey70"),
    strip.text.y.right = element_text(face = "bold", angle = 0, hjust = 0),
    panel.spacing.y = unit(0.6, "lines"),
    plot.title = element_text(face = "bold", size = 13),
    plot.title.position = "plot",
    legend.position = "none"
  )

# Figure 1: a single participant -----------------------------------------------

latent_plot <- ggplot(latent_one, aes(time, value, colour = trait)) +
  geom_line(linewidth = 0.6) +
  facet_wrap(vars(trait), ncol = 1, strip.position = "right") +
  scale_colour_manual(values = trait_colours) +
  scale_x_continuous(expand = expansion(mult = 0.01)) +
  labs(
    title = expression(bold("Latent traits") ~ ~ eta[list(1, t)]),
    x = "Day", y = NULL
  ) +
  theme_ild

item_plot <- ggplot(item_one, aes(time, value, colour = item)) +
  geom_line(linewidth = 0.4) +
  facet_wrap(vars(item), ncol = 1, strip.position = "right") +
  scale_colour_manual(values = item_colours) +
  scale_y_continuous(breaks = c(2, 4, 6)) +
  scale_x_continuous(expand = expansion(mult = 0.01)) +
  labs(
    title = expression(bold("Measured items") ~ ~ y[list(1, t)]),
    x = "Day", y = NULL
  ) +
  theme_ild

ggsave("ild_illustration.png", latent_plot + item_plot,
       width = 11, height = 7.5, dpi = 300, bg = "white")

# Figure 2: N participants overlaid ---------------------------------------------

# The participant of the first figure is drawn on top of the others, at full
# opacity. Set to NA to draw all participants alike.
highlighted <- 1L

latent_plot_n <- ggplot(
  latent_all,
  aes(time, value, colour = trait, group = interaction(participant, trait))
) +
  geom_line(
    data = subset(latent_all, !participant %in% highlighted),
    linewidth = 0.3, alpha = 0.3
  ) +
  geom_line(
    data = subset(latent_all, participant %in% highlighted),
    linewidth = 0.6, alpha = 1
  ) +
  facet_wrap(vars(trait), ncol = 1, strip.position = "right") +
  scale_colour_manual(values = trait_colours) +
  scale_x_continuous(expand = expansion(mult = 0.01)) +
  labs(
    title = expression(bold("Latent traits") ~ ~ eta[list(1, it)]),
    x = "Day", y = NULL
  ) +
  theme_ild

item_plot_n <- ggplot(
  item_all,
  aes(time, value, colour = item, group = interaction(participant, item))
) +
  geom_line(
    data = subset(item_all, !participant %in% highlighted),
    linewidth = 0.2, alpha = 0.25
  ) +
  geom_line(
    data = subset(item_all, participant %in% highlighted),
    linewidth = 0.5, alpha = 1
  ) +
  facet_wrap(vars(item), ncol = 1, strip.position = "right") +
  scale_colour_manual(values = item_colours) +
  scale_x_continuous(expand = expansion(mult = 0.01)) +
  labs(
    title = expression(bold("Measured items") ~ ~ y[list(1, it)]),
    x = "Day", y = NULL
  ) +
  theme_ild

ggsave("ild_illustration_n.png", latent_plot_n + item_plot_n,
       width = 11, height = 7.5, dpi = 300, bg = "white")

# simulate_ild_illustration.R --------------------------------------------------
#
# Illustration of intensive longitudinal data for a single participant.
#
# Three latent traits (pain, stress, sleep) follow a VAR(1) process over 100
# timepoints, and each trait is measured by three items through a factor
# measurement model. The figure shows the latent traits in the left column and
# the nine observed items in the right column. Colour identifies the latent
# trait; within the item column, the three shades of a trait's colour identify
# its three indicators.
#
# Run from the project root. Writes ild_illustration.png.
# ------------------------------------------------------------------------------

library(ggplot2)
library(patchwork)

set.seed(2026)

# Simulation settings ----------------------------------------------------------

n_timepoints <- 100
trait_names <- c("Pain", "Stress", "Sleep")
n_traits <- length(trait_names)
n_items <- 3L # items per latent trait

# VAR(1) transition matrix: moderate autoregression, some cross-lagged effects.
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

# Simulate the latent traits ---------------------------------------------------

# Draws from N(0, Psi) via the Cholesky factor, so that no extra package is
# needed.
chol_psi <- chol(Psi)
innovations <- matrix(rnorm(n_timepoints * n_traits), nrow = n_timepoints) %*% chol_psi

eta <- matrix(NA_real_, nrow = n_timepoints, ncol = n_traits,
              dimnames = list(NULL, trait_names))
eta[1, ] <- innovations[1, ]
for (t in seq(2, n_timepoints)) {
  eta[t, ] <- Phi %*% eta[t - 1, ] + innovations[t, ]
}

# Simulate the observed items --------------------------------------------------

item_names <- paste(rep(trait_names, each = n_items), rep(seq_len(n_items), n_traits))

y <- matrix(NA_real_, nrow = n_timepoints, ncol = n_traits * n_items,
            dimnames = list(NULL, item_names))
for (j in seq_len(n_traits)) {
  for (k in seq_len(n_items)) {
    y[, (j - 1) * n_items + k] <-
      nu + lambda[k] * eta[, j] + rnorm(n_timepoints, sd = sigma)
  }
}

# Data frames for plotting -----------------------------------------------------

latent_df <- data.frame(
  time = rep(seq_len(n_timepoints), times = n_traits),
  trait = factor(rep(trait_names, each = n_timepoints), levels = trait_names),
  value = as.vector(eta)
)

item_df <- data.frame(
  time = rep(seq_len(n_timepoints), times = n_traits * n_items),
  trait = factor(rep(rep(trait_names, each = n_items), each = n_timepoints),
                 levels = trait_names),
  item = factor(rep(item_names, each = n_timepoints), levels = item_names),
  value = as.vector(y)
)

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

# Plot -------------------------------------------------------------------------

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

latent_plot <- ggplot(latent_df, aes(time, value, colour = trait)) +
  geom_line(linewidth = 0.6) +
  facet_wrap(vars(trait), ncol = 1, strip.position = "right") +
  scale_colour_manual(values = trait_colours) +
  scale_x_continuous(expand = expansion(mult = 0.01)) +
  labs(
    title = expression(bold("Latent traits") ~ ~ eta[list(1, t)]),
    x = "Day", y = NULL
  ) +
  theme_ild

item_plot <- ggplot(item_df, aes(time, value, colour = item)) +
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

ild_plot <- latent_plot + item_plot

ggsave("ild_illustration.png", ild_plot,
       width = 11, height = 7.5, dpi = 300, bg = "white")

# Cellucid's own category colors: the color the viewer gives category i of a
# field that brings no colors of its own (no Scanpy <key>_colors, no
# palette), as "#rrggbb". The web app's getCategoryColor
# (cellucid/assets/js/data/palettes.js): its 72 designed colors (the category
# palette, then the color picker's), then hues fanned out around the color
# wheel by the golden angle, each channel the byte it is drawn in
# (round(channel x 255), halves up: utils/hex-color.js hexOfRgb01).
# cellucid-python keeps the same rule (src/cellucid/_category_colors.py); the
# shared fixture category_default_colors.json (written by the web app's
# scripts/generate-category-default-colors.mjs) pins the first 1,024 colors
# byte for byte in all three repositories. An apply writes them as a user
# field's colors when the session gives none (R/session-apply.R
# .apply_shown_colors), and as a recolored dataset column's when it holds no
# colors of its own (.apply_column_colors): what the viewer showed.

.category_color_seeds <- c(
  "#ff9900", "#0077c8", "#ad50e0", "#ffe500", "#00b783", "#e61977", "#4cafee", "#d15f00",
  "#009e9e", "#7d57c2", "#fa8032", "#5c6cd1", "#008b6e", "#dd2e5b", "#74ade1", "#ecb100",
  "#9339a7", "#3497b7", "#cb4b2e", "#7080c0", "#00a891", "#9294e3", "#c29411", "#347789",
  "#1f77b4", "#4169e1", "#00bfff", "#87ceeb", "#0000cd", "#191970", "#2ca02c", "#008000",
  "#32cd32", "#90ee90", "#006400", "#2e8b57", "#d62728", "#ff0000", "#ff6347", "#ff7f50",
  "#dc143c", "#b22222", "#ff7f0e", "#ffa500", "#ffd700", "#ffff00", "#bdb76b", "#ffe4b5",
  "#9467bd", "#800080", "#8a2be2", "#ee82ee", "#e377c2", "#ff69b4", "#8c564b", "#d2691e",
  "#a0522d", "#8b4513", "#deb887", "#f5deb3", "#7f7f7f", "#a9a9a9", "#d3d3d3", "#696969",
  "#2f4f4f", "#000000", "#17becf", "#00ced1", "#40e0d0", "#00ffff", "#48d1cc", "#5f9ea0"
)

.category_golden_angle <- 137.50776405003785

# The colors Cellucid gives categories 0, 1, ... count - 1 of a field that
# brings none, as "#rrggbb".
.default_category_colors <- function(count) {
  if (!is.numeric(count) || length(count) != 1L || is.na(count) || count < 0 || count != floor(count)) {
    .abort("A count of categories is a non-negative whole number.")
  }
  if (count == 0) return(character())
  index <- seq_len(count) - 1
  out <- character(count)
  seeded <- index < length(.category_color_seeds)
  out[seeded] <- .category_color_seeds[index[seeded] + 1]
  offset <- index[!seeded] - length(.category_color_seeds)
  if (length(offset) > 0L) {
    clamp <- function(value, low, high) pmin(high, pmax(low, value))
    hue <- ((offset * .category_golden_angle) %% 360) / 360
    saturation <- clamp(0.62 + 0.1 * sin(offset * 0.8), 0.45, 0.82)
    lightness <- clamp(0.54 + 0.12 * cos(offset * 0.6), 0.38, 0.72)
    chroma <- (1 - abs(2 * lightness - 1)) * saturation
    sector <- (hue %% 1) * 6
    second <- chroma * (1 - abs((sector %% 2) - 1))
    zero <- rep(0, length(offset))
    red <- ifelse(sector < 1, chroma, ifelse(sector < 2, second, ifelse(sector < 4, zero, ifelse(sector < 5, second, chroma))))
    green <- ifelse(sector < 1, second, ifelse(sector < 3, chroma, ifelse(sector < 4, second, zero)))
    blue <- ifelse(sector < 2, zero, ifelse(sector < 3, second, ifelse(sector < 5, chroma, second)))
    light <- lightness - chroma / 2
    byte <- function(channel) as.integer(floor((channel + light) * 255 + 0.5))
    out[!seeded] <- sprintf("#%02x%02x%02x", byte(red), byte(green), byte(blue))
  }
  out
}

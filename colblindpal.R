# From : https://github.com/jrnold/ggthemes/blob/master/R/colorblind.R#L1

#' Colorblind Color Palette (Discrete) and Scales

colorblind <- c(black="#000000",
orange="#E69F00",
sky_blue="#56B4E9",
bluish_green="#009E73",
yellow="#F0E442",
blue="#0072B2",
vermillion="#D55E00",
reddish_purple="#CC79A7")

#' An 8-color colorblind safe qualitative discrete palette.
#'
#' @rdname colorblind
#' @references
#' Chang, W. "\href{http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette}{Cookbook for R}"
#'
#' \url{http://jfly.iam.u-tokyo.ac.jp/color}
#'
#' @export
#' @inheritParams ggplot2::scale_colour_hue
#' @family colour
#' @seealso The \pkg{dichromat} package, \code{\link[scales]{dichromat_pal}},
#'  and \code{\link{scale_color_tableau}} for other colorblind palettes.

colorblind_pal <- function() {
  manual_pal(colorblind)
}

#' @rdname colorblind
#' @export
scale_colour_colorblind <- function(...) {
  discrete_scale("colour", "colorblind", colorblind_pal(), ...)
}

#' @rdname colorblind
#' @export
scale_color_colorblind <- scale_colour_colorblind

#' @rdname colorblind
#' @export
scale_fill_colorblind <- function(...) {
  discrete_scale("fill", "colorblind", colorblind_pal(), ...)
}


#' @examples
#' library(ggplot2)
#' library(scales)
#' show_col(colorblind_pal()(8))
#' p <- ggplot(mtcars) + geom_point(aes(x = wt, y = mpg,
#'      colour=factor(gear))) + facet_wrap(~am)
#' p + theme_gray + scale_colour_colorblind()
#' 
#' 
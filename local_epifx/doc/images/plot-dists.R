#!/usr/bin/Rscript --vanilla
##
## Plot example probability density functions for the observation models
## provided by epifx.
##

## Probability mass function for the Beta-binomial distribution.
dbetabinom <- function(x, prob, size, theta) {
    v <- lfactorial(size) - lfactorial(x) -
        lfactorial(size - x) - lbeta(theta * (1 - prob), theta * prob) +
        lbeta(size - x + theta * (1 - prob), x + theta * prob)
    exp(v)
}

plot_sample_counts <- function() {
    xs <- seq(1, 500)
    mu <- 200
    denom <- 10000
    pr <- mu / denom
    ns <- c(50, 100, 500, 1000, 5000)

    df <- NULL
    for (n in ns) {
        ys <- dbetabinom(xs, prob = pr, size = denom, theta = n)
        df <- rbind(df, data.frame(x = xs, y = ys, n = n))
    }

    df$n <- factor(df$n)

    p <- ggplot(df, aes(x = x, y = y, colour = n)) +
        geom_line(size = 1) +
        geom_vline(xintercept = mu, linetype = 'dashed') +
        scale_colour_brewer(expression(k), palette = 'Set1') +
        scale_x_continuous(breaks = pretty_breaks(5)) +
        xlab(expression(c[t])) +
        ylab(bquote(PMF ~~ (E*group("[",c[t],"]") == .(mu)))) +
        theme_rgm(key.label = 1, hide.title = FALSE,
                  legend.position = c(1, 1))

    ## Plot dimensions in inches and plot resultion in dots per inch.
    w_inch <- 8
    h_inch <- 6
    dpi <- 150

    ## Print the plot
    CairoPNG('sample-counts-dist.png', bg = 'transparent',
             width = w_inch * dpi, height = h_inch * dpi, dpi = dpi)
    print(p)
    invisible(dev.off())
}

plot_popn_counts <- function() {
    xs <- seq(1, 500)
    mu <- 200
    sizes <- c(1000, 100, 25, 5, 1)

    df <- NULL
    for (s in sizes) {
        ys <- dnbinom(xs, mu = mu, size = s)
        df <- rbind(df, data.frame(x = xs, y = ys, k = s))
    }

    df$k <- factor(df$k)

    p <- ggplot(df, aes(x = x, y = y, colour = k)) +
        geom_line(size = 1) +
        geom_vline(xintercept = mu, linetype = 'dashed') +
        scale_colour_brewer(expression(k), palette = 'Set1') +
        xlab(expression(y[t])) +
        ylab(bquote(PMF ~~ (E*group("[",y[t],"]") == .(mu)))) +
        theme_rgm(key.label = 1, hide.title = FALSE,
                  legend.position = c(1, 1))

    ## Plot dimensions in inches and plot resultion in dots per inch.
    w_inch <- 8
    h_inch <- 6
    dpi <- 150

    ## Print the plot
    CairoPNG('popn-counts-dist.png', bg = 'transparent',
             width = w_inch * dpi, height = h_inch * dpi, dpi = dpi)
    print(p)
    invisible(dev.off())
}

main <- function() {
    req_pkgs <- c('ggplot2', 'scales', 'themergm', 'Cairo')
    for (pkg in req_pkgs) {
        if (! require(pkg, character.only = TRUE, quietly = TRUE)) {
            cat('ERROR: this script requires "', pkg, '" to be installed.\n',
                sep = '')
            quit(status = 2)
        }
    }

    ## Produce a plot for each observation model.
    plot_sample_counts()
    plot_popn_counts()
}

main()

set.seed(11)
xx <- rnorm(n = 2500, mean = 1, sd = 1)
zz <- rnorm(2500, 0, 1)
logit_y <- .25*xx + .25*zz + .05*xx*zz
prob_y <- exp(logit_y) / (1 + exp(logit_y))
yy <- rbinom(length(prob_y), 1, prob_y)

test_that("errors occur as intended", {
  expect_error(litl(x = 'a', y = yy, z = zz,
                    ntiles_x = 8, ntiles_z = 3),
               'x must be double, integer, or logical')
  expect_error(litl(x = xx, y = 'a', z = zz,
                    ntiles_x = 8, ntiles_z = 3),
               'y must be logical \\(T\\/F\\/NA\\) or binary numeric \\(1\\/0\\/NA\\)')
  expect_error(litl(x = xx, y = yy, z = 'a',
                    ntiles_x = 8, ntiles_z = 3),
               'z, x, and y must be the same length')
  expect_error(litl(x = xx, y = yy, z = zz,
                    ntiles_x = 'a', ntiles_z = 3),
               'ntiles_x must be double or integer')
  expect_error(litl(x = xx, y = yy, z = zz,
                    ntiles_x = 8, ntiles_z = 3, plot_line = 'a'),
               'plot_line must be "lm", "loess", or NULL')
})

test_that("warnings occur as intended", {
  expect_warning(litl(x = xx, y = yy, z = zz,
                      ntiles_x = 8, ntiles_z = 'a'),
                 'ntiles_z is NULL or non-numeric\\: dividing z by terciles')
  expect_warning(litl(x = c(rep(NA, 100), xx),
                      y = c(rep(NA, 100), yy),
                      z = c(rep(NA, 100), zz),
                      ntiles_x = 8, ntiles_z = 3),
                 'input variables contained NAs: dropping \\d+ incomplete cases')
  expect_warning(litl(x = xx[1:60], y = yy[1:60], ntiles_x = 4),
                 '\\d+ ntile.+')
})

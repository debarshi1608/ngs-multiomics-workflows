# Source PDF: Debarshi_011_Lab1.pdf
x <- runif(1000, 50, 100)
z <- runif(1000, 50, 10000)

log_x <- log10(x)

log_z <- log10(z)

d1 <- density(log_x)
d2 <- density(log_z)

plot(
  d1,
  xlim = range(c(0,5)),
  ylim = range(c(0,4)),
  col = "blue",
  lwd = 2,
  main = "Density Plot after Log10 Transformation",
  xlab = "Log10(Value)",
  ylab = "Density"
)

lines(d2, col = "red", lwd = 2)

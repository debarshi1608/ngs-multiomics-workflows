# Source PDF: Debarshi_011_Lab1.pdf
# Generate 100 random numbers from Uniform Distribution
x <- runif(100, min = 50, max = 100)

# Generate 100 random numbers from Normal Distribution
y <- rnorm(100, mean = 70, sd = 70)

print(x)
print(y)

cat("Uniform Distribution\n")
cat("Mean:", mean(x), "\n")
cat("Median:", median(x), "\n\n")

cat("Normal Distribution\n")
cat("Mean:", mean(y), "\n")
cat("Median:", median(y), "\n")

table_x <- table(x)
names(table_x)[which.max(table_x)]

plot(
  density(x),
  main = "Density Plot - Uniform Distribution",
  xlab = "Values",
  ylab = "Density",
  col = "blue",
  lwd = 2
)

plot(density(y),
  main = "Density Plot - Normal Distribution",
  xlab = "Values",
  ylab = "Density",
  col = "red",
  lwd = 2
)

plot(
  density(x),
  col = "blue",
  lwd = 2,
  xlim = c(5,1500),
  main = "Overlay of Uniform and Normal Density Plots",
  xlab = "Values",
  ylab = "Density"
)

lines( density(y),col = "red",lwd = 2)

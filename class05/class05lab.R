#' ---
#' title: "Week 3 Data Visualization Lab"
#' author: "Yane Lee (PID A17670350"
#' ---

# Week 3 Data Visualization Lab

# Install the package ggplot
#install.packages("ggplot2")
library(ggplot2)
#View(cars)

#A quick base R plot - this is not gglot
plot(cars)

# Our first ggplot, we need data + aes + geoms
ggplot(data=cars) +
  aes(x=speed, y=dist) +
  geom_point()

p <- ggplot(data=cars) + 
  aes(x=speed, y=dist) +
  geom_point()

# Add a line with geom with geom_line()
p + geom_line()

# Add a trend line close to the data
p + geom_smooth()

p + geom_smooth(method="lm") 

#---------------------------------------------------

# Read in our drug expression data
url <- "https://bioboot.github.io/bimm143_S20/class-material/up_down_expression.txt"
genes <- read.delim(url)
head(genes)

# Q. How many genes are in this dataset?
nrow(genes)

# Q. How many 'up' regulated genes?
table( genes$State )

# Q. What fraction of total genes is up-regulated?
round((table(genes$State) / nrow(genes)) * 100, 2)

# Let's make a first plot attempt
g <- ggplot(data=genes) + aes(x=Condition1, y=Condition2, col=State) + geom_point()

g
# Add some color
g + scale_color_manual(values=c("pink", "lightyellow", "lightblue")) +
  labs(title="Gene expression changes", x= "Control (no drug)", y= "Drug Treatment") +
  theme_bw()


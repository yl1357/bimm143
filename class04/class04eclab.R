head(cdc$height)
tail(cdc$height, 20)
plot(cdc$height, cdc$weight)
cor(cdc$height, cdc$weight)
hist(cdc$height)
weight_kg <- cdc$weight * 0.454
height_m <- cdc$height * 0.0254
bmi <- (weight_kg)/(height_m^2)
plot(cdc$height, bmi)
cor(cdc$height, bmi)
head(bmi >= 30, 100)
plot(cdc[1:100, 5], cdc[1:100, 6])
cdc[bmi >= 30, 9]
table(cdc$gender[bmi >= 30])
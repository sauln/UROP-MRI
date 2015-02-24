# comment

data(container.df)
fit = hotelling.test(.~gp, data = container.df)
fit

subs.df = container.df[1:10,]
subs.df$gp = rep(1:2, c(5,5))
fitPerm = hotelling.test(Al+Fe~gp, data = subs.df, perm = TRUE)
fitPerm
plot(fitPerm)

data(bottle.df)
fit12 = hotelling.test(.~Number, data = bottle.df)
fit12

fit23 = hotelling.test(.~Number, data = bottle.df, pair = c(2,3))
print(fit23)

print("hello world")


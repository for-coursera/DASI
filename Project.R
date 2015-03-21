load(url("http://bit.ly/dasi_gss_data"))
names(gss)
system("open https://d396qusza40orc.cloudfront.net/statistics%2Fproject%2Fgss1.html")
cats <- c("consci", "degree", "year")
my_na <- gss[cats]
my <- na.omit(my_na)
head(my)
my <- subset(my, my$year > 2005)
plot(my[,1:2])

my.mat <- as.matrix(table(my[1:2])) ; my.mat
barplot(prop.table(my.mat, margin=2))

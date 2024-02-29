
library(dplyr)
library(ggplot2)

x <- seq(0.5, 100, 0.5 )

y1 <- x^(-2)
y2 <- x^(-3)
y3 <- x^(-1)
y4 <- x^(-0.5)

data <- data.frame(x, y1, y2, y3, y4)

the_plotg <- ggplot(data) + 
  geom_line(aes(x=x, y=y2), color = 'red') + 
  geom_line(aes(x=x, y=y1), color = 'orange') +
  geom_line(aes(x=x, y=y3), color = 'blue') +
  geom_line(aes(x=x, y=y4), color = 'green') + 
  ylim(0,2) +
  xlim(0, 50)
the_plotg

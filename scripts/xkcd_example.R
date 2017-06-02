 library(extrafont)
 library(ggplot2)

if ('xkcd' %in% fonts()) {

  xrange <- range(mtcars$mpg)

  yrange <- range(mtcars$wt)

   ggplot(mtcars) +
    geom_line(aes(mpg, wt)) +
    xkcdaxis(xrange, yrange) +
    theme(text = element_text(size = 16, family = 'xkcd'))

  }else {
    +   warning("Not xkcd fonts installed!")
    + p <- ggplot() + geom_point(aes(x=mpg, y=wt), data=mtcars) +}
p

download.file("http://simonsoftware.se/other/xkcd.ttf", dest="xkcd.ttf", mode="wb")

system("mkdir ~/.fonts")

system("cp xkcd.ttf  ~/.fonts")

font_import(pattern = "[X/x]kcd", prompt=FALSE)

fonts()

fonttable()

if(.Platform$OS.type != "unix") {
  loadfonts(device="win")
} else {
  loadfonts()
}


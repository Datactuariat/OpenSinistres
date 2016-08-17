# devtools::install_github("hrbrmstr/waffle")
library(waffle)

# data
nombre <- c(
  'Responsabilité civile corporels'=2,
  'Responsabilité civile matériels'=18,
  'Vol et incendie'=3,
  'Bris de glaces'=33,
  'Dommages aux véhicules'=30,
  'Tempête, Grêle, Neige'=3,
  'Autres'=11
)

charge <- c(
  'Bris de glaces'=9,
  'Dommages aux véhicules'=27,
  'Responsabilité civile matériels'=16,
  'Autres'=4,
  'Vol et incendie'=7,
  'Tempête, Grêle, Neige'=4,
  'Responsabilité civile corporels'=33
)

col_sia=rgb(99,24,66,maxColorValue=255)
waffle_pal2 <- rev(c(col_sia,rev(c("#587C92","#ADC6AB","#E6ECC5","#ECDFBC","#CB9C87","#985960"))))
waffle_pal2 <- rev(waffle_pal)

waffle(sort(nombre,decreasing=TRUE), rows=5, color=waffle_pal,
       title="Nombre de sinistres  automobile en 2010", 
       xlab="un carré représente 77k sinistres")

waffle(charge, rows=5, color=waffle_pal,
       title="Charge de sinistres  automobile en 2010", 
       xlab="un carré représente 110 m€")



rij <- c("R", "W",	"87.5/12.5",	"75/25",	"62.5/37.5",	"50/50",	"37.5/62.5",	"25/75", "12.5/87.5",	"D")
opzet <- sample(rij)
for (i in 1:11) {
  print(i)
  opzet <- cbind(opzet, sample(rij))
}

rij <- c("R", "W",	"75/25",	"50/50",	"25/75",	"D")
opzet <- sample(rij)
for (i in 1:16) {
  print(i)
  opzet <- cbind(opzet, sample(rij))
}



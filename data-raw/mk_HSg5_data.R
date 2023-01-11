
setwd('D:/Sequoia/Docs/Data/2022-12')
setwd('../data')


Geno_HSg5 <- SimGeno(Ped = Ped_HSg5, nSnp = 200, ParMis=0.4,
                     CallRate = 0.9, SnpError = 0.005)
save(Geno_HSg5, file="Geno_HSg5.rda")


SeqOUT_HSg5 <- sequoia(GenoM = Geno_HSg5, LifeHistData = LH_HSg5, Module = "ped",
                       Err = 0.005)
save(SeqOUT_HSg5, file="SeqOUT_HSg5.rda")


PC <- PedCompare(Ped_HSg5, SeqOUT_HSg5$Pedigree)

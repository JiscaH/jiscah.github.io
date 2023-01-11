# - create mock pedigree & life history data for 'griffins'
# - simulate genotype data from this pedigree
# - create full set of output: sequoia() + GetMaybeRel() + EstConf()

# assuming wd is /sequoia


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LH_griffin ====

# Generate N=200 individuals with specified birth years & random sex
BY <- rep(c(2001:2010), each=20)
Sex <- sample.int(n=2, size=200, replace=TRUE)
ID <- paste0("i", formatC(1:200, width=3, flag="0"), "_", BY, "_", ifelse(Sex==1, "F", "M"))

LH_griffin <- data.frame(ID, Sex, BirthYear = BY)
write.table(LH_griffin[, 1:3],
            'data/LH_griffin.txt',
            row.names=F, sep="\t", quote=FALSE)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Ped_griffin ====

# For each individual, sample the age of its dam & sire with specified probabilities
Age.dam <- sample(1:3, size=200, prob=c(0.6, 0.3, 0.1), replace=TRUE)
Age.sire <- sample(1:3, size=200, prob=c(0.3, 0.5, 0.2), replace=TRUE)

# calculate parental birth years
BY.dam <- BY - Age.dam
BY.sire <- BY - Age.sire

# For each individual, sample a parent from the calculated parental birth year
# If the parental birth year is before the earliest birth year of 'known'/'sampled'
# individuals, the parent is NA.
Ped_griffin <- data.frame(id = ID, dam = NA, sire = NA,
                          stringsAsFactors=FALSE)
for (i in 1:200) {
  Ped_griffin$dam[i] <- ifelse(BY.dam[i] < 2001, NA,
                               sample(ID[BY == BY.dam[i] & Sex==1], size=1))
  Ped_griffin$sire[i] <- ifelse(BY.sire[i] < 2001, NA,
                                sample(ID[BY == BY.sire[i] & Sex==2], size=1))
}


Ped_griffin <- cbind(Ped_griffin, birthyear = BY)
write.table(Ped_griffin[, 1:4], "data/Ped_griffin.txt",
            row.names=FALSE, sep="\t", quote=FALSE)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Geno_griffin ====

# Genotypes accidentally not saved for version 2.0.7 example

# ensure same individuals are genotyped vs not in 2.4 as in 2.0.7
# write.table(sequoia::SeqOUT_griffin$PedigreePar[,'id'],
#             'genotyped_individuals_griffin.txt',
#             row.names=FALSE, sep="\t", quote=FALSE, col.names=FALSE)
genotyped <- sequoia::SeqOUT_griffin$PedigreePar[['id']]

Geno_griffin <- SimGeno(Ped_griffin, nSnp=400,
                        CallRate = c(setNames(rep(0.99, length(genotyped)), genotyped),
                                     setNames(rep(0, 200-length(genotyped)),
                                              setdiff(Ped_griffin$id, genotyped))))
# default values: MAF=0.3, CallRate=0.99, SnpError=5e-4, ErrorFM='version2.0'

save(Geno_griffin, file="data/Geno_griffin.rda")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SeqOUT_griffin ====

SeqOUT_A <- sequoia(GenoM = Geno_griffin,
                          LifeHistData = LH_griffin,
                          Module = 'ped')


chk <- PedCompare(sequoia::SeqOUT_griffin$Pedigree,
                  SeqOUT_A$Pedigree)
# different (sequoia version 2.0.7 vs 2.4)



SeqOUT_griffin <- SeqOUT_A
save(SeqOUT_griffin, file="data/SeqOUT_griffin.rda")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MaybeRel_griffin ====

# MaybeRel_griffin <- GetMaybeRel(GenoM = Geno_griffin,
#                                 SeqList = SeqOUT_griffin,
#                                 Module = 'ped')
# very quick --> not necessary to include in pkg.

MaybeRel_griffin <- GetMaybeRel(GenoM = Geno_griffin, Err=0.001, Module = 'par')

save(MaybeRel_griffin, file="data/MaybeRel_griffin.rda")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Conf_griffin ====

Conf_griffin <- EstConf(Pedigree = SeqOUT_griffin$Pedigree,
                        LifeHistData = LH_griffin,
                        args.sim = list(nSnp = 400, SnpError = 0.001, ParMis=0.4),
                        args.seq = list(Module = 'ped', Err=0.001),
                        nSim = 20,
                        nCores = 5,
                        quiet = TRUE)

save(Conf_griffin, file="data/Conf_griffin.rda")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FieldMums_griffin ====

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# base function 'replace' with match only replaces first match.
Replace <- function(V, old, new) {
  if (length(old) != length(new))  stop("'old' and 'new' must have same length")
  if (!all(old %in% V))  stop("all 'old' must be in V")
  these <- lapply(seq_along(old), function(x, y=V) which(y == old[x]))
  newr <- rep(new, sapply(these, length))
  replace(V, unlist(these), newr)
}


# version 2.0.7
# Replacements <- matrix(c("F0001", "GreenYellow",
#                          "F0002", "YellowBlue",
#                          "F0003", "OrangeGreen",
#                          "F0004", "RedOrange",
#                          "F0006", "PinkBlue",
#                          "F0007", "YellowPink",
#                          "i081_2005_F", "GreenBlue"),7,2,byrow=TRUE)

# version 2.4 (changes based on pedcompare() $DummyMatch)
Replacements <- matrix(c("F0003", "GreenYellow",
                         "F0005", "YellowBlue",
                         "F0001", "OrangeGreen",
                         "F0006", "RedOrange",
                         "F0002", "PinkBlue",
                         "F0007", "YellowPink",
                         "i081_2005_F", "GreenBlue"),7,2,byrow=TRUE)


FieldMums_griffin <- SeqOUT_griffin$Pedigree[, c("id", "dam")]
FieldMums_griffin$mum <- Replace(FieldMums_griffin$dam,
                                 old = Replacements[,1],
                                 new = Replacements[,2])
FieldMums_griffin$mum[substr(FieldMums_griffin$mum, 1, 2) == "F0"] <- NA
FieldMums_griffin <- FieldMums_griffin[!substr(FieldMums_griffin$id, 1, 2) %in% c("F0", "M0"), ]


# 1) split sibship F0001
FieldMums_griffin$mum[FieldMums_griffin$id %in% c("i158_2008_M",
                                                  "i160_2008_F",
                                                  "i170_2009_F",
                                                  "i173_2009_F")] <- "BlueRed"

# 2) mum swap
MumsX <- FieldMums_griffin[!is.na(FieldMums_griffin$mum), ]
View(MumsX[order(MumsX$mum), ])
View(FieldMums_griffin[!is.na(FieldMums_griffin$mum), ])

# FieldMums_griffin[FieldMums_griffin$id %in% c("i177_2009_M", "i175_2009_M"), ]
# not: use as mother-offspring example (Ped2-only)

FieldMums_griffin[FieldMums_griffin$id %in% c("i147_2008_F", "i148_2008_F"), ]
#              id   dam         mum
# 101 i147_2008_F F0003 GreenYellow
# 102 i148_2008_F F0005  YellowBlue

FieldMums_griffin$mum[FieldMums_griffin$id == "i147_2008_F"] <- "YellowBlue"
FieldMums_griffin$mum[FieldMums_griffin$id == "i148_2008_F"] <- "GreenYellow"

# perfect: colour misread, handwriting misread, or lab swap?


# 4) P1only: non-SNPd sib(s) & non-SNPd mum
PCX <- PedCompare(Ped_griffin, SeqOUT_griffin$Pedigree)
View(PCX$MergedPed[, c("id", "id.r", "dam.1", "dam.2", "id.dam.cat")])

Ped_griffin[which(Ped_griffin$dam == "i038_2002_F"), ]
#             id         dam        sire birthyear
# 53 i053_2003_M i038_2002_F i008_2001_M      2003  # not SNPd
# 57 i057_2003_M i038_2002_F i018_2001_M      2003  # SNPd

FieldMums_griffin$mum[FieldMums_griffin$id == "i057_2003_M"] <- "GreenRed"
FieldMums_griffin <- merge(FieldMums_griffin,
                           data.frame(id = "i053_2003_M",
                                      dam = NA,
                                      mum = "GreenRed",
                                      stringsAsFactors=FALSE),
                           all=TRUE)


write.table(FieldMums_griffin[, c("id", "mum")],
            "data/FieldMums_griffin.txt", sep="\t", row.names=FALSE)






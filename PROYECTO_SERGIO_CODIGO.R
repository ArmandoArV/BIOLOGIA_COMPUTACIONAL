cat("\014")

# Librerías

# library(seqinr)
# library(ggmsa)


# origina = read.fasta("")
#delta = read.fasta("deltaS.fasta")
#omicron = read.fasta("omicronS.fasta")

# Las length
# length(delta)
# length(delta)
# length(omicron)
trad =    c(UUU="F", UUC="F", UUA="L", UUG="L",
            UCU="S", UCC="S", UCA="S", UCG="S",
            UAU="Y", UAC="Y", UAA="STOP", UAG="STOP",
            UGU="C", UGC="C", UGA="STOP", UGG="W",
            CUU="L", CUC="L", CUA="L", CUG="L",
            CCU="P", CCC="P", CCA="P", CCG="P",
            CAU="H", CAC="H", CAA="Q", CAG="Q",
            CGU="R", CGC="R", CGA="R", CGG="R",
            AUU="I", AUC="I", AUA="I", AUG="M",
            ACU="T", ACC="T", ACA="T", ACG="T",
            AAU="N", AAC="N", AAA="K", AAG="K",
            AGU="S", AGC="S", AGA="R", AGG="R",
            GUU="V", GUC="V", GUA="V", GUG="V",
            GCU="A", GCC="A", GCA="A", GCG="A",
            GAU="D", GAC="D", GAA="E", GAG="E",
            GGU="G", GGC="G", GGA="G", GGG="G")

df = data.frame(
  Mutation = character(),
  Nucleotide = numeric(),
  Codon = character(),
  Protein = character(),
  index = numeric()
)


genOmi = solA
genDelta = solB
# genOmi = Omicron
# genDelta = Delta
genOmi
genDelta

genOmi[genOmi=="G"] = "C" # T por U
# attr1 = attr(genOmi, "Annot")
genDelta[genDelta=="T"] = "U" # T por U
diferentes = which(genOmi != genDelta);


# vec = unlist(strsplit(attr1,"\\[|\\]|:|=|\\."));
# vec[which(vec=="gene")+1]
# vec[which(vec=="location")+1]


if (length(diferentes) > 0){
  # print(paste(length(diferentes),"mutaciones de nucleótido único"))
  for (k in diferentes){
    mutation = paste(genOmi[k], "to", genDelta[k], sep="")
    inicio = k - (k-1) %% 3
    global = 21563+inicio # 21563 valor original
    index = as.integer(k/3+1)  # %/%
    codonOmi = paste(genOmi[inicio], genOmi[inicio+1], genOmi[inicio+2],sep="")
    codonDelta = paste(genDelta[inicio], genDelta[inicio+1], genDelta[inicio+2],sep="")
    codonChange = paste(codonOmi,global,codonDelta, sep="")
    aminoChange = paste(trad[codonOmi],index,trad[codonDelta], sep="")
    print(paste(mutation, global, codonChange, aminoChange, index))
    newRow = list(mutation, global, codonChange, aminoChange, index)
    df[nrow(df)+1, ] = newRow
  }
}
# codonChange
k
mutation
inicio
index
str(df)
df
View(df)


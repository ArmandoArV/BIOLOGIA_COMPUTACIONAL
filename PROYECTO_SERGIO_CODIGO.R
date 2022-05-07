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
trad =    c("GAC"="D","GAU"="D",
            "GAA"="E","GAG"="E",
            "CGA"="R","CGC"="R","CGG"="R","CGU"="R","AGA"="R","AGG"="R",
            "AAA"="K","AAG"="K",
            "AAC"="N","AAU"="N",
            "CAC"="H","CAU"="H",
            "CAA"="Q","CAG"="Q",
            "UCA"="S","UCC"="S","UCG"="S","UCU"="S","AGC"="S","AGU"="S",
            "ACA"="T","ACC"="T","ACG"="T","ACU"="T",
            "GCA"="A","GCC"="A","GCG"="A","GCU"="A",
            "GGA"="G","GGC"="G","GGG"="G","GGU"="G",
            "GUA"="V","GUC"="V","GUG"="V","GUU"="V",
            "CCA"="P","CCC"="P","CCG"="P","CCU"="P",
            "CUA"="L","CUC"="L","CUG"="L","CUU"="L","UUA"="L","UUG"="L",
            "UUC"="F","UUU"="F",
            "UAC"="Y","UAU"="Y",
            "AUA"="I","AUC"="I","AUU"="I",
            "AUG"="M",
            "UGG"="W",
            "UGC"="C","UGU"="C")



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

# genOmi[genOmi=="T"] = "U" # T por U
genOmi[genOmi=="T"] = "U" # T por U
# attr1 = attr(genOmi, "Annot")
# genDelta[genDelta=="T"] = "U" # T por U
genDelta[genDelta=="T"] = "U" # T por U

diferentes = which(genOmi != genDelta);


# vec = unlist(strsplit(attr1,"\\[|\\]|:|=|\\."));
# vec[which(vec=="gene")+1]
# vec[which(vec=="location")+1]


if (length(diferentes) > 0){
  # print(paste(length(diferentes),"mutaciones de nucleótido único"))
  for (k in diferentes){
    mutation = paste(genOmi[k], " to ", genDelta[k], sep="")
    inicio = k - (k-1) %% 3
    global = 0+inicio # 21563 valor original
    index = as.integer(k/3+1)  # %/%
    codonOmi = paste(genOmi[inicio], genOmi[inicio+1], genOmi[inicio+2],sep="")
    codonDelta = paste(genDelta[inicio], genDelta[inicio+1], genDelta[inicio+2],sep="")
    codonChange = paste(codonOmi," to ",codonDelta, sep="")
    aminoChange = paste(trad[codonOmi]," to ",trad[codonDelta], sep="")
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

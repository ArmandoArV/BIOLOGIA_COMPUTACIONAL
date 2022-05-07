cat("\f")
library(seqinr)

Delta = read.fasta("deltaS.fasta")
Omicron = read.fasta("omicronS.fasta")
Delta
Omicron




A = toupper(Delta[[1]])[1:100]
B = toupper(Omicron[[1]])[1:100]
ADN_to_ARNm = function(nucleotido){
  # return (switch(nucleotido,"C"="C","G"="G","T"="U","A"="A"))
  #return (switch(nucleotido,"C"="G","G"="C","T"="A","A"="U"))
  return (switch(nucleotido,"C"="C","G"="G","T"="U","A"="A"))
}

A = as.vector(sapply(A, ADN_to_ARNm))
B = as.vector(sapply(B, ADN_to_ARNm))
# A
# B

m = matrix(data=0,nrow=length(B)+1, ncol=length(A)+1)
m[1, ] = seq(0,length(A))*-2
m[ ,1] = seq(0,length(B))*-2
m

CalcularPeso = function(fila, col){
  if (A[col-1]==B[fila-1]) diagonal = m[fila-1,col-1] + 1
  else diagonal = m[fila-1,col-1] - 1
  up = m[fila-1,col] - 2
  left = m[fila,col-1] - 2
  peso = max (diagonal, up, left)
  return (peso)
}

for (fila in seq(2,length(B)+1)){
  for (col in seq(2,length(A)+1)){
    m[fila, col] = CalcularPeso(fila, col) 
  }
}
m

fila = length(B)+1
col = length(A)+1
solA = solB = c()
while(fila>1 || col>1){
  if (col>1 && fila>1 && A[col-1]==B[fila-1]){
    solA = c(A[col-1], solA)
    solB = c(B[fila-1], solB)
    col = col - 1
    fila = fila -1
  }else{
    if (col>1 && fila>1 && m[fila-1, col-1]>m[fila, col-1] && m[fila-1, col-1]>m[fila-1, col]){
      solA = c(A[col-1], solA)
      solB = c(B[fila-1], solB)
      col = col - 1
      fila = fila -1
    }else if (fila==1 || (col>1 && m[fila, col-1] > m[fila-1, col])){
      solA = c(A[col-1], solA)
      solB = c("_", solB)
      col = col - 1
    }else if (col==1 || (fila>1 && m[fila, col-1] <= m[fila-1, col])){
      solA = c("_", solA)
      solB = c(B[fila-1], solB)
      fila = fila - 1 
    }
  }
}
sol = matrix(data=c(solA,solB),nrow = 2, ncol = length(solA), byrow = TRUE)
sol


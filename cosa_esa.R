AMINOACIDS = c(UUU="FPhe", UUC="FPhe", UUA="LLeu", UUG="LLeu", UCU="SSer", UCC="SSer", UCA="SSer", UCG="SSer", UAU="YTyr", UAC="YTyr", UAA='-', UAG='-', UGU="CCys", UGC="CCys", UGA='-', UGG="WTrp", CUU="LLeu", CUC="LLeu", CUA="LLeu", CUG="LLeu", CCU="PPro", CCC="PPro", CCA="PPro", CCG="PPro", CAU="HHis", CAC="HHis", CAA="QGln", CAG="QGln", CGU="RArg", CGC="RArg", CGA="RArg", CGG="RArg", AUU="IIle", AUC="IIle", AUA="IIle", AUG="MMet", ACU="TThr", ACC="TThr", ACA="TThr", ACG="TThr", AAU="NAsn", AAC="NAsn", AAA="KLys", AAG="KLys", AGU="SSer", AGC="SSer", AGA="RArg", AGG="RArg", GUU="VVal", GUC="VVal", GUA="VVal", GUG="VVal", GCU="AAla", GCC="AAla", GCA="AAla", GCG="AAla", GAU="DAsp", GAC="DAsp", GAA="EGlu", GAG="EGlu", GGU="GGly", GGC="GGly", GGA="GGly", GGG="GGly");
globalIndex = 0;
mutations = data.frame(mutation = character(), nucleo = numeric(), codon = character(), protein = character(), index = numeric(), structure = character());
install.packages("seqinr");
library("seqinr");
for(structure in 1:length(read.fasta("original.txt"))){
  original = toupper(read.fasta("deltaS.fasta")[[structure]]);
  example = toupper(read.fasta("omicronS.fasta")[[structure]]);
  original[original=='T'] = 'U';
  example[example=='T'] = 'U';
  metadata = unlist(strsplit(attr(original,"Annot"), "\\[|\\]|:|=|\\.|join|\\(|\\)"));
  metadata = metadata[metadata!="" & metadata!=" "];
  globalIndex = as.integer(metadata[which(metadata=="location") + 1]);
  for(index in seq_along(example)){
    if(original[index]!=example[index]){
      codon = ceiling(index/3);
      originalCodon = original[(3*codon-2):(3*codon)];
      originalCodonDenomination = AMINOACIDS[paste(originalCodon, sep="", collapse = "")];
      exampleCodon = example[(3*codon-2):(3*codon)];
      exampleCodonDenomination = AMINOACIDS[paste(exampleCodon, sep="", collapse = "")];
      if(substr(originalCodonDenomination, 1, 1) != substr(exampleCodonDenomination, 1, 1))
        mutations[nrow(mutations)+1, ] = list(paste(original[index], "to", example[index], sep=""), globalIndex, paste(paste(originalCodon, collapse = ""), "to", paste(exampleCodon, collapse = ""), sep=""), paste(substr(originalCodonDenomination, 1, 1), "to", substr(exampleCodonDenomination, 1, 1), sep=""), codon, metadata[which(metadata=="protein") + 1]);
    }
    globalIndex = globalIndex + 1;
  }
}

View(mutations);


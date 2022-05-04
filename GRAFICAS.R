library(ggplot2)
p = ggplot(df)
p = p + aes(x=Mutation,fill=Mutation, label=..count..)
p = p + ggtitle("Mutaciones Delta-Omicrón")
p = p + labs(x="Mutación", y="Frecuencia", fill="Mutaciones")
p = p + geom_bar(stat = "count")
p = p + geom_text(stat = "count", vjust=1)
p = p + facet_grid(~index)
p
  

#Filter Singletons/Doubletons 

# A) Prevalence filtering

# Define prevalence of each taxa - (in how many samples did each taxa appear at least once)


# Find and remove singletons

prev0 = apply(X=otu_table(psq), MARGIN=ifelse(taxa_are_rows(psq), yes=1, no=2), FUN=function(x){sum(x > 0)})


# Define prevalence threshold as 1% of total samples

nsamples(psq)

0.01 * 100 = 1

psq

# Apply prevalence filter, using `prune_taxa()` function
# trippleton-species removed:

psq3 = prune_taxa((prev0 > 3), psq)

psq3
# mapping file: gg-map-adults.txt
# OTU table (classic): GG_100nt_even10k-adults.txt
# OTU table (biom): GG_100nt_even10k-adults.biom
# genus-level taxonomy: taxa/GG_100nt_even10k-adults_L6.txt
# multi-level taxonomy: taxa/merged-taxa.txt
# distance matrix: beta/unweighted_unifrac_GG_100nt_even10k-adults.txt
# pcoa file: beta/unweighted_unifrac_GG_100nt_even10k-adults_pcoa.txt
# Classification column: COUNTRY
# Regression column: AGE
# feature stats table (qvalues): stats/taxon-stats-table.txt


# GG input data were obtained from Antonion Gonzales Pena for 
# CHM "Rethinking Enterotypes" paper.
# They were picked againts greengenes 13_8 and rarefied at 10k.

# commands used to generate downstream tables:
filter_samples_from_otu_table.py -m gg-map.txt -s 'Age:*,!None' -i GG_100nt_even10k.biom -o GG_100nt_even10k-adults.biom
filter_otus_from_otu_table.py -s 20 -i GG_100nt_even10k-adults.biom -o GG_100nt_even10k-adults-s20.biom
summarize_taxa.py -i GG_100nt_even10k-adults.biom -L '2,3,4,5,6,7' -o taxa
sh ../../bin/util/merge_taxa taxa/GG_100nt_even10k-adults_L*.txt > taxa/merged-taxa.txt
beta_diversity.py -i GG_100nt_even10k-adults.biom -o beta -t 97_otus.tree 
principal_coordinates.py -i beta/unweighted_unifrac_GG_100nt_even10k-adults.txt -o beta/unweighted_unifrac_GG_100nt_even10k-adults_pcoa.txt

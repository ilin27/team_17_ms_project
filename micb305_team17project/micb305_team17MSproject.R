#load in the libraries
library(phyloseq)
library(tidyverse)
library(vegan)
library(ape)

#load in the files
ms_meta <- read_delim("MS_files_export/corrected_ms_metadata.tsv", delim="\t")
ms_otu <- read_delim("MS_files_export/feature-table.txt", delim="\t", skip=1)
ms_tax <- read_delim("MS_files_export/taxonomy.tsv", delim="\t")
ms_tree <- read.tree("MS_files_export/tree.nwk")

#formatting files
ms_otu_mat <- as.matrix(ms_otu[,-1])
rownames(ms_otu_mat) <- ms_otu$`#OTU ID`
MS_OTU <- otu_table(ms_otu_mat, taxa_are_rows = TRUE)
class(MS_OTU)

# add allergy status (TRUE if allergy_specific is filled in)
ms_meta <- mutate(ms_meta, allergy_status = !is.na(allergy_specific))

ms_meta_df <- as.data.frame(ms_meta[,-1])
rownames(ms_meta_df) <- ms_meta$`sample-id`
MS_META <- sample_data(ms_meta_df)
class(MS_META)

ms_tax_mat <- ms_tax |>
  select(-Confidence) |>
  separate(col=Taxon, sep="; ",
           into = c("Domain",
                    "Phylum",
                    "Class",
                    "Order",
                    "Family",
                    "Genus",
                    "Species"),
           fill = "right"
  ) |>
  as.matrix()
ms_tax_mat <- ms_tax_mat[,-1]
rownames(ms_tax_mat) <- ms_tax$`Feature ID`
MS_TAX <- tax_table(ms_tax_mat)
class(MS_TAX)

# combine OTU, META, TAX, tree into a single object
ms_ps <- phyloseq(MS_OTU, MS_META, MS_TAX, ms_tree)
ms_ps

# filtering out mitochondria and chloroplasts
ms_ps_filt <- ms_ps %>%
  subset_taxa(Domain == "d__Bacteria" &
                Class != "c__Chloroplast" &
                Family != "f__Mitochondria") %>%
  filter_taxa(function(x) sum(x) > 5, prune = TRUE) %>%
  prune_samples(sample_sums(.) > 100, .)

sample_sums(ms_ps_filt)
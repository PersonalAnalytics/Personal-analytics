"""

OVERVIEW: 

Python module for extracting taxonomic information from closed-reference OTU tables in classic dense format.

"""

import os, sys
import numpy as np


def convert_GGIDs_to_latin_names(OTU_table, new_OTU_table, taxonomy_file=0):
    # Takes as input a closed-reference OTU table with GreenGenes IDs and returns the same table but with latin names instead of GG IDs.
    if(taxonomy_file == 0):
        taxonomy_file = '/home/ubuntu/gg_13_5_otus/taxonomy/97_otu_taxonomy.txt'

    # Extract OTU IDs
    with open(OTU_table, 'r') as fid:
        OTU_table_lines = fid.readlines()
        OTU_IDs = [line.split()[0] for line in OTU_table_lines]
        OTU_IDs = OTU_IDs[1:]

    # Extract latin names for OTU IDs
    with open(taxonomy_file, 'r') as fid:
        all_lines = fid.readlines()
        gg_OTU_IDs = [line.split()[0] for line in all_lines]
        indices_of_OTUs_in_gg = [gg_OTU_IDs.index(OTU) for OTU in OTU_IDs]
        OTU_latin_names = ["".join(all_lines[i].split()[1:]) for i in indices_of_OTUs_in_gg]

    # Rename OTUs by latin names in OTU table
    for i in range(1,len(OTU_table_lines)):
        line = OTU_table_lines[i].split()
        line[0] = OTU_latin_names[i-1]
        OTU_table_lines[i] = '\t'.join(line)

    # Rewrite OTU table
    with open(new_OTU_table, 'w') as fid:
        fid.write(OTU_table_lines[0])
        for i in range(1,len(OTU_table_lines)):
            fid.write(OTU_table_lines[i] + '\n')

def collapse_taxonomic_contents(OTU_table, taxonomic_level, parent_node=0):
    # Takes as input an OTU table with latin names as OTU IDs,
    # the taxonomic level to which all taxonomies should be collapsed.
    # Accepts: kingdom, phylum, class, order, family, genus, species
    # e.g.  collapse_taxonomic_contents(OTU_table, class, Firmicutes)

    with open(OTU_table, 'r') as fid:
        OTU_table_lines = fid.readlines()
        OTU_IDs = [line.split()[0] for line in OTU_table_lines]
        OTU_IDs = OTU_IDs[1:]

    # Collapse to the right level        
    if(taxonomic_level == "kingdom"):
        OTU_taxa = [OTU_ID.split(';')[0] for OTU_ID in OTU_IDs]
    if(taxonomic_level == "phylum"):
        OTU_taxa = [';'.join(OTU_ID.split(';')[:2]) for OTU_ID in OTU_IDs]
        OTU_taxa_parents = [OTU_ID.split(';')[0][3:] for OTU_ID in OTU_IDs]
    if(taxonomic_level == "class"):
        OTU_taxa = [';'.join(OTU_ID.split(';')[:3]) for OTU_ID in OTU_IDs]
        OTU_taxa_parents = [OTU_ID.split(';')[1][3:] for OTU_ID in OTU_IDs]
    if(taxonomic_level == "order"):
        OTU_taxa = [';'.join(OTU_ID.split(';')[:4]) for OTU_ID in OTU_IDs]
        OTU_taxa_parents = [OTU_ID.split(';')[2][3:] for OTU_ID in OTU_IDs]
    if(taxonomic_level == "family"):
        OTU_taxa = [';'.join(OTU_ID.split(';')[:5]) for OTU_ID in OTU_IDs]
        OTU_taxa_parents = [OTU_ID.split(';')[3][3:] for OTU_ID in OTU_IDs]
    if(taxonomic_level == "genus"):
        OTU_taxa = [';'.join(OTU_ID.split(';')[:6]) for OTU_ID in OTU_IDs]
        OTU_taxa_parents = [OTU_ID.split(';')[4][3:] for OTU_ID in OTU_IDs]
    if(taxonomic_level == "species"):
        OTU_taxa = [';'.join(OTU_ID.split(';')[:7]) for OTU_ID in OTU_IDs]
        OTU_taxa_parents = [OTU_ID.split(';')[5][3:] for OTU_ID in OTU_IDs]

    # Get indices of each unique taxon
    taxa_indices = {}
    for i in range(len(OTU_taxa)):
        if (OTU_taxa[i] not in taxa_indices) and (OTU_taxa_parents[i] == parent_node) and (OTU_taxa[i][(len(OTU_taxa[i])-2):] != "__"):
            taxa_indices[OTU_taxa[i]] = []
        if (OTU_taxa_parents[i] == parent_node) and (OTU_taxa[i][(len(OTU_taxa[i])-2):] != "__"):
            taxa_indices[OTU_taxa[i]].append(i)            
    
    # Read in OTU abundances
    sample_labels = OTU_table_lines[0].split()[1:]
    nOTUs = len(OTU_IDs)
    nSamples = len(sample_labels)

    # Load OTU table row major order
    OTU_table_counts = np.zeros((nOTUs, nSamples))
    for i in range(1,nOTUs+1):
        OTU_table_counts[i-1,:] = OTU_table_lines[i].split('\t')[1:]
    
    # Normalize counts by sample (column)
    for i in range(nSamples):
        OTU_table_counts[:,i] = np.divide(OTU_table_counts[:,i], np.sum(OTU_table_counts[:,i]))
    
    # Get sample contents for each taxa of the chosen level
    taxa_dict = {}
    for key in taxa_indices:
        taxa_dict[key] = {}
        indices = taxa_indices[key]
        for i in range(nSamples):
            taxa_dict[key][sample_labels[i]] = np.sum(OTU_table_counts[indices,i])

    return taxa_dict





# imports
import numpy as np
from Bio import Phylo
import pandas as pd

# Create dictionary of amino acids and corresponding codons
aas = {'Phe': ['UUU', 'UUC'], 'Leu': ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'], 'Ile': ['AUU', 'AUC', 'AUA'], 'Met': ['AUG'], 'Val': ['GUU', 'GUC', 'GUA', 'GUG'],
       'Ser': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'], 'Pro': ['CCU', 'CCC', 'CCA', 'CCG'], 'Thr': ['ACU', 'ACC', 'ACA', 'ACG'], 'Ala': ['GCU', 'GCC', 'GCA', 'GCG'],
       'Tyr': ['UAU', 'UAC'], '***': ['UAA', 'UAG', 'UGA'], 'Trp': ['UGG'], 'His': ['CAU', 'CAC'], 'Gln': ['CAA', 'CAG'], 'Asn': ['AAU', 'AAC'], 'Lys': ['AAA', 'AAG'],
       'Cys': ['UGU', 'UGC'], 'Arg': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 'Gly': ['GGU', 'GGC', 'GGA', 'GGG'], 'Asp': ['GAU', 'GAC'], 'Glu': ['GAA', 'GAG']}


def normalize_codom_bias(df):
    '''Normalizes the codon usage by amino acid for a dataframe'''
    for animal in df['SpeciesName']:
        species_bool = df['SpeciesName'] == animal
        species = df[species_bool]
        for aa in aas:
            total = 0

            for codon in aas[aa]:
                total += species[codon]

            # loop through codons, normalize by amino acid, and replace in data frame
            for codon in aas[aa]:
                df.loc[species_bool, [codon]] = species[codon]/total


def species_diff(df, species1, species2):
    '''Finds the difference between the codon bias of two species'''
    diff = 0

    animal1 = df[df['SpeciesName'] == species1]
    animal2 = df[df['SpeciesName'] == species2]

    # Find the squared difference between each species' total codon bias
    for aa in aas:
        for codon in aas[aa]:
            diff += (animal1[codon].values[0] - animal2[codon].values[0])**2

    return diff


def create_diff_codon_dict(df):
    '''Creates a matrix of differences between the species in a 
    dataframe. Stored in a dictionary of dictionaries based on
    the names of the provided species.'''
    diff_dict = {}

    # Loop over all the animal pairs to find their differences
    for species1 in df['SpeciesName'].unique():
        for species2 in df['SpeciesName'].unique():
            if(species1 == species2):
                continue

            diff = species_diff(df, species1, species2)

            if(not species1 in diff_dict):
                diff_dict[species1] = {}

            if(not species2 in diff_dict):
                diff_dict[species2] = {}

            diff_dict[species1][species2] = diff
            diff_dict[species2][species1] = diff

    return diff_dict


def hemo_differences(seq1, seq2):
    '''Counts the number of pairwise differences between
    two sequences'''
    count = 0

    # Find differences between animo acid sequences
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            count += 1

    return count


def create_hemo_diff_dict(animals, nodes):
    '''Creates a matrix of differences between the given species.
    Stored in a dictionary of dictionary based on the indices of 
    the provided list.'''
    diff_dict = {}

    # Loop over all the animal pairs to find their differences
    for i, animal1 in enumerate(nodes):
        for j, animal2 in enumerate(nodes):

            # Don't check animals against themselves
            if(animal1 == animal2):
                continue

            diff = hemo_differences(animals[i], animals[j])

            # Update difference matrix with distance calculated
            if(not animal1 in diff_dict):
                diff_dict[animal1] = {}

            if(not animal2 in diff_dict):
                diff_dict[animal2] = {}

            diff_dict[animal1][animal2] = diff
            diff_dict[animal2][animal1] = diff

    return diff_dict

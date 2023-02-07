from Bio import Phylo
import pandas as pd

import species
import upgma


def main():
    # Read file into a data frame
    f = open('../data/codon_usage.csv', 'r')
    df = pd.read_csv(f)

    # Create three phylogenetic trees
    create_rodents_tree(df)
    create_hemoglobin_animal_tree()
    create_codon_animal_tree(df)


def create_rodents_tree(df):
    # Create the dataframe of rodents
    rodents_with_mit = df[df['Kingdom'] == 'rod']
    rodents = rodents_with_mit[rodents_with_mit['DNAtype'] == 0]
    rodents = rodents[8:13]
    species.normalize_codom_bias(rodents)

    # Create the phylogenetic tree
    diff_dict = species.create_diff_codon_dict(rodents)
    index_tree = upgma.generate_phylogenetic_tree(diff_dict)

    # Generate and visualize phylogenetic tree
    upgma.write_tree_xml(index_tree, 'rodent_tree')
    tree = Phylo.read('../data/rodent_tree.xml', 'phyloxml')
    tree.ladderize()
    Phylo.draw(tree)


def create_hemoglobin_animal_tree():
    # Create lists of all hemoglobin sequences
    human = "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHF.DLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR"
    gorilla = ".VLSPADKTNVKAAWGKVGAHAGDYGAEALERMFLSFPTTKTYFPHF.DLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR"
    cow = "MVLSAADKGNVKAAWGKVGGHAAEYGAEALERMFLSFPTTKTYFPHF.DLSHGSAQVKGHGAKVAAALTKAVEHLDDLPGALSELSDLHAHKLRVDPVNFKLLSHSLLVTLASHLPSDFTPAVHASLDKFLANVSTVLTSKYR"
    horse = "MVLSAADKTNVKAAWSKVGGHAGEYGAEALERMFLGFPTTKTYFPHF.DLSHGSAQVKAHGKKVGDALTLAVGHLDDLPGALSNLSDLHAHKLRVDPVNFKLLSHCLLSTLAVHLPNDFTPAVHASLDKFLSSVSTVLTSKYR"
    donkey = "MVLSAADKTNVKAAWSKVGGNAGEFGAEALERMFLGFPTTKTYFPHF.DLSHGSAQVKAHGKKVGDALTLAVGHLDDLPGALSNLSDLHAHKLRVDPVNFKLLSHCLLSTLAVHLPNDFTPAVHASLDKFLSTVSTVLTSKYR"
    rabbit = ".VLSPADKTNIKTAWEKIGSHGGEYGAEAVERMFLGFPTTKTYFPHF.DFTHGSZQIKAHGKKVSEALTKAVGHLDDLPGALSTLSDLHAHKLRVDPVNFKLLSHCLLVTLANHHPSEFTPAVHASLDKFLANVSTVLTSKYR"
    carp = "MSLSDKDKAAVKGLWAKISPKADDIGAEALGRMLTVYPQTKTYFAHWADLSPGSGPVKKHGKVIMGAVGDAVSKIDDLVGGLAALSELHAFKLRVDPANFKILAHNVIVVIGMLYPGDFPPEVHMSVDKFFQNLALALSEKYR"

    animals = [carp, cow, donkey, horse, human, gorilla, rabbit]
    nodes = ["carp", "cow", "donkey", "horse", "human", "gorilla", "rabbit"]

    # Create the phylogenetic tree
    diff_dict = species.create_hemo_diff_dict(animals, nodes)
    index_tree = upgma.generate_phylogenetic_tree(diff_dict)

    # Generate and visualize phylogenetic tree
    upgma.write_tree_xml(index_tree, 'test_tree')
    tree = Phylo.read('../data/test_tree.xml', 'phyloxml')
    tree.ladderize()
    Phylo.draw(tree)


def create_codon_animal_tree(df):
    # Create the dataframe of animals (carp has been replaced with trout)
    animals = ['Homo sapiens', 'Gorilla gorilla', 'Oryctolagus cuniculus',
               'Bos taurus', 'Equus asinus', 'Equus caballus', 'Oncorhynchus mykiss']

    for animal in animals:
        if animal == animals[0]:
            class_animals = df[df['SpeciesName'] == animal]
        else:
            class_animals = pd.concat(
                [class_animals, df[df['SpeciesName'] == animal]])

    species.normalize_codom_bias(class_animals)

    # Create the phylogenetic tree
    diff_dict = species.create_diff_codon_dict(class_animals)
    index_tree = upgma.generate_phylogenetic_tree(diff_dict)

    # Generate and visualize phylogenetic tree
    upgma.write_tree_xml(index_tree, 'animal_tree')
    tree = Phylo.read('../data/animal_tree.xml', 'phyloxml')
    tree.ladderize()
    Phylo.draw(tree)


if __name__ == '__main__':
    main()

# Codon Bias Phylogeny
This project was completed by <a href=https://github.com/mayamoweryevans>@mayamoweryevans</a> and <a href=https://github.com/aidanccraft>@aidanccraft</a> as the final project in our Quantitative Biology class. The purpose of this project is to investigate if <a href=https://en.wikipedia.org/wiki/Codon_usage_bias>codon bias</a> is an indicator for phylogenetic relationships. If codon bias does prove to be an indicator, this could be used to increase understanding of the relationship between different organisms. A csv file containing the codon bias for many different organisms was taken from <a href=www.kaggle.com/datasets/salikhussaini49/codon-usage>Kaggle</a>. To make the project more manageable, rodents were selected to be the focus of the project. The data was analyzed by creating a data frame and using pandas techniques. To create a phylogenetic tree, the algorithm UPGMA was used.

# Usage
The project can be run either using the Jupyter Notebook found in the ```notebooks``` folder or by running the ```main.py``` file from the ```src``` folder. In order to run the project <a href=https://pandas.pydata.org/>pandas</a>, <a href=https://numpy.org/>NumPy</a>, and <a href=https://biopython.org/>BioPython</a> must be installed.

# Results and Discussion
<p align="center">
  <img src="https://user-images.githubusercontent.com/29548845/217389829-d693eee4-c715-48b8-9981-2893018121c2.png"></img>
  <p align="center">A phylogenetic tree created from a simple amino acid sequence alignment on hemoglobin sequences</p>
</p>
<p align="center">
  <img src="https://user-images.githubusercontent.com/29548845/217389827-689c48e4-12a7-4633-81e3-557bd7344bc8.png"></img>
  <p align="center">A phylogentic tree of the same animals in the previous figure using codon bias data</p>
</p>
<p align="center">
  <img src="https://user-images.githubusercontent.com/29548845/217389830-2e3b37a0-6f3a-4c84-87d6-a47da9e05f77.svg"></img>
  <p align="center">A phylogentic tree created from a random assortment of rodents using codon bias data</p>
</p>
Overall, this method seemed to work on some animals. The tree produced for rodents was similar to accepted phylogenetic trees and the tree produced for in-class animals also had some resemblance to the tree created in class. For the rodents, Peromyscus leucopus and Peromyscus maniculatus should have been closely related, but they were instead placed at opposite ends of the tree. The most abnormal result for the in-class animals was the donkey. The donkey appeared to be distantly related to the rest of the mammals--more distantly than the trout. It also put cows and horses to be closely related unlike in class. However, some sources suggests that this may need corrected. The rodents saw similar results: most of the tree matched accepted relationships with a few vastly more distantly related than the should be. One reason for these discrepancies could be incomplete data. For many animals, the codon bias in the original csv file only accounts for a small fraction of the genome. Should more of the genome be analyzed, it is possible that this would improve the results. Another concern is that the distances for the rodents is much greater than the distances for in class animals. Theoretically, the rodents should be more closely related to each other than the in class animals are. This could also be caused by the incomplete data. Another potential source of error is an issue in the distance function that has not been accounted for. Overall, codon bias seems to be correlated to phylogenetics however more testing is required to determine the other factors at play.

# Next Steps
Future steps in this project could explore this method with a larger data set to confirm if the results hold true. The mitochondrial data could be used to see if that would be able to create a more accurate phylogenetic tree. Mitochondrial DNA could also be compared to bacteria to track the relationship between mitochondria and modern bacteria. Finally, the discrepancies in the tree created could be explored to determine the cause of the discrepancies. A different data set could be used which accounts for a larger percentage of the genome. This would determine if these discrepancies are due to incomplete data or another factor. Once this is determined, it would direct the next steps of this project as the relationship between codon bias and phylogenetics are explored.

# Sources
Codon Bias Dataset: www.kaggle.com/datasets/salikhussaini49/codon-usage

Rodent Phylogeny: www.ncbi.nlm.nih.gov/pmc/articles/PMC3532383/

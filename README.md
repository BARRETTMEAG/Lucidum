This project aims to identify which species among primates (50 taxa) possess or lack tapetum lucidum. The genes that we are using are: IRBP, MC1R, RAG2, Mitochondria, and G6PD.

The folders are broken down as follows:

0. Reading includes 15 papers and 3 peer reviews that have been found during researching this project. 

1. Datasets contains: raw nucleotide data for different species that our group has gathered. The R studio code has been included, showing how these nucleotide data was taken, trimmed, and then became Aligned. 

2. Alignments contains the genetics after trimming and aligning in R. 
SuperMatrix contains the alignment of all of the taxa. 

3. Draft trees shows trees pertaining to different genes, including the supermatrix (SM). 

   To create the supermatrix, this website and program was super helpful: https://github.com/andrewbudge/Liger. Thank you, Andrew!

4. Final Project contains the final trees found in Final SM and tree data containing the fasta file and taxa list.
   
   - Final SM contains trees from Beauti-Beast, NGPhylogeny, Iqtree, and the pictures from all of these files.
    
      - NGPhylogeny are the beginning trees. They were created in NGPhlogeny and although they do not     contain bootstrap values, their structure is still appliciable.
      
      - Beauti-Beast contains the SM files from Beauti and Beast, both the regular and the long version that we learned in Taming of the Beast (https://taming-the-beast.org/tutorials/Introduction-to-BEAST2/). 
      
      - Iqtree is the computer program that allowed for getting the Bayesian posterior values that helped us know the confidence of the trees.  



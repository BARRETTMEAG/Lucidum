Lucidum Overview:
    
    - This project aims to identify which species among primates (50 taxa) possess or lack tapetum lucidum. The genes that we are using are: IRBP, MC1R, RAG2, Mitochondria, and G6PD.

The folders are broken down as follows:

   0. Reading includes 15 papers and 3 peer reviews that have been found during researching this project. 

   1. Datasets contains: raw nucleotide data for different species that our group has gathered. The R studio code has been included, showing how these nucleotide data was taken, trimmed, and then became Aligned. 

   2. Alignments contains the genetics after trimming and aligning in R. 
      - SuperMatrix contains the alignment of all of the taxa. 

   3. Draft trees shows trees pertaining to different genes, including the supermatrix (SM). 
       - To create the supermatrix, this website and program was super helpful: https://github.com/andrewbudge/Liger. Thank you, Andrew!

   4. Final Project contains the final trees found in Best_Tree, SM_R_Code (code used to create the final tree) and SM_tree_data containing the fasta file and taxa list.
   
      - Final SM contains trees from NGPhylogeny R Code, Iqtree, and the pictures from all of these files.
    
      - NGPhylogeny are the beginning trees. They were created in NGPhlogeny and although they do not     contain bootstrap values, their structure is still appliciable.
       
      - SM_ Iqtree is the computer program that allowed for getting the Bayesian posterior values that helped us know the confidence of the trees.  



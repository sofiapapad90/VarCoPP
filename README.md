# VarCoPP

## Contents

1. Summary - What is this code about?
1. Prerequisites
1. Code usage - Predictions on the validation sets with **VarCoPP.py**
1. Code usage - Train the model with **VarCoPP_train.py**
1. Code usage - Assess performance with **VarCoPP_cross_validation.py**
1. Flexibility and Hydrophobicity annotation with **calc_aa.diff.py**
1. Contributors
1. License 



## 1. Summary - What is this code about?

This is the source code for **VarCoPP**, the Variant Combination Pathogenicity Predictor, which predicts the pathogenicity of
bi-locus variant combinations (*i.e.* variants between two genes). It is based on the in-house developed code that
was used to produce the performance results presented in the paper: Papadimitriou, S. *et al.* Predicting disease-causing variant combinations, *Proceedings of the National Academy of Sciences*. May 2019, DOI: https://doi.org/10.1073/pnas.1815601116.

With the current code you can:
1. re-train VarCoPP using the same annotated training datasets that were used in the paper
(present in the **./training_data** directory) or
2. try the current VarCoPP model on the annotated validation data used in the paper
(Datasets S1 - S5, an example of Dataset S1 is available in the **./example_data** directory). 

VarCoPP provides as output a .txt file containing a list of predictions for each digenic combination, including a final
class label (“**disease-causing**” or “**neutral**”) for each instance. This list is ranked based on the **Support Score** of 
each prediction (*i.e.* how many RFs agree for the disease-causing class) and the **Classification Score** (*i.e.* the median 
probability for the disease-causing class among all RFs), with the most probable disease-causing bi-locus combinations being 
at the top of the list.

**NOTE: this code requires already annotated data. An online version of VarCoPP that automatically annotates and predicts variant combinations is available at: https://varcopp.ibsquare.be. There, you can predict bi-locus combinations from a given list of variants (SNPs and indels) belonging to one individual.**


## 2. Prerequisites
* **Python (>=3.6)**. Installation information on: https://www.python.org/downloads/

* **NumPy (>=1.8.2) and SciPy (>=0.13.3)**. These can be installed using the SciPy pack: https://scipy.org/install.html

* **scikit-learn** python package (>=0.18.1). Installation:`$ pip install -U scikit-learn` or `$ conda install scikit-learn`.
Further information is found on: http://scikit-learn.org/stable/install.html

**NOTE:** The results of this work were produced using scikit-learn version 0.18.1. If you would like to replicate the 
results presented in the paper, you are encouraged to use this version to avoid potential compatibility errors.
	   
For scikit-learn version of >= 0.19 you may get the following warning when running VarCoPP.py:  

>UserWarning: Trying to unpickle estimator DecisionTreeClassifier from
>version 0.18.1 when using version 0.19.0. This might lead to breaking code or
>invalid results. Use at your own risk.




## 3. Code usage - Predictions on the validation sets with VarCoPP.py

It is advised that all input files are in the same directory as VarCoPP.py. Otherwise, the path of the files should also
be provided in the arguments. All output is by default stored in the same directory as VarCoPP.py.

#### INPUT: 
To run VarCoPP.py on a validation set (for e.g. the Dataset S1 present in the **./example_data** directory), 
you need: 
1. a .txt tab-delimited file with the features for each combination. This file should contain, first, a column with the combination ID and then the columns for the 
11 annotated features in the order they appear in the dataset (Flex1, Hydr1, CADD1, CADD2, CADD3, CADD4, HI_A, Rec_A, HI_B, 
Rec_B, Biol_Dist).
2. the trained VarCoPP model "**VarCoPP_model.p.gz**”

#### OUTPUT:
The “predictions.txt” file with the Classification and Support Score of each combination, and its predicted
class. The combinations in the file are ranked based on their **Support Score** (*i.e.* how many RFs agree for the disease-causing
class) and their **Classification Score** (*i.e.* the median probability for the disease-causing class among all RFs), with the 
most probable disease-causing combinations being at the top of the list.

#### USAGE EXAMPLES :
1. See help and input parameters of the script:
`$ python VarCoPP.py -h`

2. Run the script with default settings (the ones used in the paper) on the example validation combinations (Dataset S1 in the paper)
`$ python VarCoPP.py -v ./example_data/validation_combinations.txt`




## 4. Code usage - Train the model with VarCoPP_train.py

For validation purposes, we also provide the source code (VarCoPP_train.py) that can train a new model of the predictor,
in the **./training_data** directory. The training data sets are also present in this folder.

It is advised that all required input files are in the same directory as VarCoPP_train.py. Otherwise, the path of the files
should also be provided in the arguments. All output is by default stored in the same directory as VarCoPP_train.py.


#### INPUT :
To run VarCoPP_train.py you need: 
1. the "**1KGP_training_features.txt**" file with the 1KGP training sets
(200 combinations per set) present in the **./training_data** directory and the
2. "**DIDA_training_features.txt**" file with the DIDA training set (200 combinations) present in the same directory


#### OUTPUT:
A cPickle (.p) file with the trained model in the form of a LIST of 500 Random Forests. The output 
file is ready to be used in VarCoPP.py.


#### USAGE EXAMPLES :
1. See help and input parameters of the script:
`$ python VarCoPP_train.py -h`

2. Re-train the model using the same parameters (and training data sets) as the ones used in the paper:
`$ python VarCoPP_train.py`

2. In case you want to use a different output filename, you can specify it in the option -f.
`$ python VarCoPP_train.py -f my_filename`


**NOTE**: The user should keep in mind that due to the randomness of the RF procedure when training the predictor, they may
not obtain the same predictor model as the one present in the “VarCoPP_model.p.gz” with this script.
Therefore, if they attempt to run VarCoPP.py with the new trained model obtained by this script, the CS and SS scores of the
example data may very slightly differ.



## 5. Code usage - Assess performance with VarCoPP_cross_validation.py
	
	
This script can replicate the cross-validation procedure that was performed in the paper for the training set present in 
the **./training_data** directory. The validation procedure for each balanced set is the LeaveOneGenePairOut procedure,
as described in the paper.

The amount of iterations for each validation is around 310 - 311, based on the number of unique gene pairs present in each
balanced set.

This script does **NOT** store a VarCoPP model, but provides performance statistics among the different individual predictors.

It is advised that all required input files are in the same directory as VarCoPP_cross_validation.py. Otherwise, the path of the files
should also be provided in the arguments. All output is by default stored in the same directory as VarCoPP_cross_validation.py.

#### INPUT :
1. the "**1KGP_training_features.txt**" file with the 1KGP training sets (200 combinations per set)
2. the "**DIDA_training_features.txt**" file with the DIDA training set (200 combinations) and the 
3. "**training_data_pairs.txt** " file with the pair information for all combination IDs present in the training data
All files are present in the **./training_data** directory. 


#### OUTPUT :
1. '**training_cross_val_summary.txt**', with a summary of the average performance among all 500
                   individual predictors
2. '**training_cross_val_iterations.txt**', with a detailed performance for each of the 500
   individual predictors
3. '**training_metrics_dict.json**', with information on TPR, FPR, thresholds, recall, precision per RF
   (to create ROC and PR curves), as well as the division of training sets during cross-
   validation
4. '**feat_dict.json**', with the importance of all features between all 500 individual predictors to
  create box-plots


#### USAGE EXAMPLES :
1. See help and input parameters of the script:
`$ python VarCoPP_cross_validation.py -h`

2. Re-train the model using the same parameters and training data sets as the ones used in the paper:
`$ python VarCoPP_cross_validation.py`

2. In case you want to use a different output filename, you can specify it in the option -f.
`$ python VarCoPP_cross_validation.py -f my_filename`


**NOTE:** The running time for this script is around 16-20 hours, if run in a single core. The user should take this into
consideration before running the code.

**NOTE:** The user should keep in mind that due to the randomness of the RF procedure when training the predictor, they may
not obtain exactly the cross-validation result statistics as the ones presented in the paper.




## 6. Flexibility and Hydrophobicity annotation with calc_aa_diff.py

We also provide the in-house developed script calc_aa_diff.py to annotate variants with flexibility and hydrophobicity differences.

Make sure that you provide the complete path of the input files for this script (see below at the usage examples).

#### INPUT :

**1.** a tab-delimited variants file with one variant per line and columns with information in the following order:
translation start position, reference and alternative protein alleles and the gene name (or other identifier present in the fasta file,
see below). More columns can be available (such as zygosity, CADD etc.), but only the first 4 will be used in the script.

The protein starting position, reference allele and alternative allele should be implemented as in the examples below:

Type of change |  Protein position | Ref protein allele | Alt protein allele
-------------- |  ---------------- | ------------------ | ------------------
SNP,frameshift,lost start | 99 | G | A 
Stop gained | 99 | G | *
Stop lost      |        99       |              *         |          A
Insertion        |      99              |       G          |         AG
Deletion            |   99                |     G           |        -
Change plus deletion |  99          |           GTAA    |            A
Change plus insertion  | 99          |           G             |      ATAA
Synonymous mutation  |  99  |                   G     |              G
Non synonymous stop   | 99   |                  *      |             *
intronic variant    |   0    |                  -      |             -
splicing variant    |   0    |                  -      |             -

**2**. a .fasta file that contain the protein sequences of the genes involved, in the form:

	>GeneA (or other identifier present in the variants file)
	Protein sequence
	>GeneB (or other identifier present in the variants file)
	Protein sequence
	...

An "**example_variants.txt**" file and an "**example_sequences.fasta**" file are available at the **./example_data** directory.

#### OUTPUT : 
The output is the same input variants file with appended at the end two extra
	columns	for flexibility amino acid difference and hydrophobicity amino acid
	difference for each variant.


#### USAGE EXAMPLES :
1. See help and input parameters of the script:
`$ python calc_aa_diff.py -h`

2. Annotate the variants of the "**example_variants.txt**" file:
`$ python calc_aa_diff.py -v variants_file.txt -s fasta.txt`



## 7. Contributors

Main contributor for source code: Sofia Papadimitriou

People who also contributed to the development of the code:
* Andrea Gazzo
* Tom Lenaerts




## 8. License
This code is licensed under the GNU General Public License v3.0 (https://choosealicense.com/licenses/gpl-3.0/). 


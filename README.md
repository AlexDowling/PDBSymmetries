# Predicting Protein Symmetries

# Table of Contents
1. [Introduction](#Introduction)
2. [Data Collection](#Data-Collection)
4. [Modeling Approach](#Modeling-Approach)
5. [Conclusions and Future Directions](#Conclusions-and-Future-Directions)
6. [Description of Repository](#Description-of-Repository)

## Introduction
Proteins consist of a sequence of amino acids assembled into a linear chain. To function, it must be realized into a three dimensional structure, which typically exhibits some form of symmetry between identical or similar subunits. Identifying the symmetries of a protein structure has proved useful in understanding the function of the protein. 

We wish to identify if a protein has a certain symmetry (for this project, we restrict our attention to C2 cyclic symmetry). This information is useful to:
- Pharmaceutical scientists and biochemists are interested in determining the function of proteins.
- Biochemists often explore how subunits of these proteins can fit together.
- Applied mathematicians have an interest in how and which group structures can appear in biological processes.

Using experimental data collected on thousands of proteins, we aim to address the following:
- Which features of proteins predict their symmetry? Particular atom counts/positions, bonds, amino acids?
- Given new experimental data on a protein, does it have a cyclic C2 symmetry? How accurate are these predictions?

## Data Collection
The source of data for this project is the [Research Collaboratory for Structural Bioinformatics (RCSB) Protein Data Bank (PDB)](https://www.rcsb.org/), which provides annotated data on the three dimensional structure of over 90,000 proteins and other macromolecules. 
Each protein in the data bank has an associated PDB file which stores thousands of lines of data, including:
- names of scientists who discovered the protein,
- equipment used to perform the crystallography,
- number of repeated subunits in the protein,
- precise positions of atoms in the protein and their types
- bonds between atoms in the protein and their types
- amino acids and their sequencing

In addition, the PDB offers a [Python search API](https://search.rcsb.org/#search-api). This allows us to search for proteins which exhibit global symmetry and for those which exhibit cyclic C2 symmetry.
We restricted ourselves to chemical and geometric features to provide for our models. This included:
- the number of atoms of each type, (carbon, hydrogen, oxygen, etc.)
- the positions of these atoms relative to each other,
- the number of bonds of each type, (single, double, triple)
- the number of amino acids of each type,
- oligomeric count (number of repeated subunits).

There can be over 10,000 atoms in a particular protein, and this total can vary dramatically. There are also a variety of types for each atom and amino acid.
We extracted two data sets from this:
- point clouds (centered at 0 with farthest point scaled and rotated to (1, 0, 0)) storing the positions of all atoms in space (regardless of type),
- tabular data containing the number of atoms of each type, number of bonds of each type, number of amino acids of each type, and oligomeric count.

The code for the point cloud formatting is in src/ProteinPointCloud.py, and the code for the tabular data formatting is in src/DataFormatting.py.

## Modeling Approach
The classifiers we trained return a 0 or 1 if the protein data considered has or does not have C2 global symmetry respectively. Our baseline model predicts the mode of the training data set regardless of the input data.
Our approach to handle the point cloud data was centered on equipping the space of proteins with a suitable metric to set up a k-nearest neighbors classifier. The metric of choice was the Hausdorff metric, which is the greatest distance from a point in one set to the closest point in the other set. For a training set of 200 proteins, we set k=14.

|Model|Train accuracy|Test accuracy|
|---|:---:|:---:|
|Baseline (Mode)| 66%| 66%|
|Point Cloud KNN| 69%| 58%|

There are two crucial flaws with this model. 
1. Computational complexity: each protein's point cloud may contain thousands of elements, so measuring their distance with the Hausdorff metric is computationally expensive and immune to the standard optimization techniques.
2. High dimensionality: since the size and shape of point clouds arising from proteins can vary dramatically, the "space of proteins" is extremely high dimensional. To get a KNN model to run effectively, one would need to sample many more proteins to flesh out clusters in this space.

We used the software [XGBoost](https://xgboost.readthedocs.io/en/stable/) to generate a gradient boosted decision tree model for the tabular data. This contained 229 features (described in the data collection section) on roughly 5,000 proteins. The performance of the resultant model is shown below.

|Model|Train accuracy|Test accuracy|
|---|:---:|:---:|
|Baseline (Mode)| 75.92%| 76.15%|
|XGBoost| 99.97%| 90.44%|
|XGBoost without OC| 99.82%| 77.73%|

The gradient boosted decision tree model consistently performed 12-15 percentage points better than the baseline. Performing a feature importance analysis on the XGBoost model revealed that oligomeric count (number of repeated subunits in the protein) was the dominant feature in predictions:
|Rank|Feature|
|---|---|
|1.|Oligomeric count|
|2.|# of LYS acids|
|3.|# of ILE acids|
|4.|# of LEU acids|
|5.|# of GLY acids|
|6.|# of GLN acids|
|7.|# of TYR acids|
|8.|# of sulfur atoms|

However, the ranking of these features below oligomeric count were somewhat variate with choice of training data.

## Conclusions and Future Directions
The molecular structure of proteins is extremely complicated, with intricate chemical information occuring over thousands of atoms for each protein. Identifying C2 symmetry in proteins is therefore nontrivial. 

Some of the features we analyzed, such as the positions of atoms in the proteins, our methods were not able to incorporate into a compelling model without encountering issues with computational complexity.

However, coarser features like oligomeric count still proved effective in identifying C2 symmetry. The 90.41% accuracy boasted by the gradient boosted decision tree model proves its value to biochemists, pharmaceutical scientists, and applied mathematicians.

The PDB contains many more valuable features we did not exploit. This is not due to our estimation of their worth, but to a lack of subject knowledge. We are confident that introducing more of these chemical features would improve the models we have provided.
In addition, one could employ more advanced techniques (like those from topological data analysis or point cloud registration) to minimize computational complexity in measuring spatial information. We also believe this to be a promising direction.

## Description of Repository
```bash
|____README.md		Readme file.
|____LICENSE			License file.

|____notebooks
|    |____BoostAnalysis.ipynb					Notebook containing our GBDC classifier
|    |____PointCloudAnalysis.ipynb    Notebook containing our KNN model

|____proteins
|    |____batch-1       Folder containing 50 protein .cif files
|    |____batch-2       Folder containing 50 protein .cif files
|    |____large-batch   Folder containing 250 protein .cif files

|____src
|    |____DataFormattting.py 		  Python functions for formatting tabular data from proteins
|    |____ProteinPointCloud.py    Python functions for building regularized point clouds from proteins

|____symmetry-lists
|    |____C2_list.pkl    List of proteins possessing C2 global symmetry, retrieved from PDB query.

|____tabular-data
|    |____large_protein_df    pandas data frame for our XGBoost models


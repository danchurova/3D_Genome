# 3D_Genome
### 3D-Genome project is a Bioinformatics Institute summer internship project
We aimed to analyze the strongest interchromosomal interactions and analyze their role and functionality.

#### Main concepts

#### Data source
* Hi-C in situ data from GEO ([GSE63525](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525))
* Bisulfite-Seq data of K562 cell line ([GSM683856](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM683856))

#### Scripts
* **iterative_normalization.py** is a script that performs normalization in each chrom_pair interchromosomal contact matrices
we used all three normalizations: 
  - Knight-Ruiz (KR)
  - Vanilla Coverage (VC)
  - Square Root Vanilla Coverage (SQRTVC)

* **concatenate_data.py** is a script that look at the contact values frequences in the whole genome, calculate the threshold value for the one percent top-interactions (the most tight) and keep only those that bigger than this threshold

   - works with raw and all normalized data

* **jaccard_similarity.py** is a script that calculates jaccard similarity coefficient between two top interactions datasets (raw and normalized)

* **take_a_closer_look.py** is a script that allows to get all interactions of higher resolution that are within top-interactions of 100 kb data

#### Some results
##### K562

|metrics|raw|KRnormalized|VCnormalized|SQRTVCnormalized|
|---|---|---|---|---|
|initial number of contacts| 97719803   |97719803   |97719803   |97719803   |
|threshold   |6.0   |6.219   |12.863(?)   |6.546   |
|1%   |977198   |977198   |977198   |977198   |
|sum bigger than 1%   |1236906|977222   |977226   | 977270|


#### References
1. Rao SS, Huntley MH, Durand NC, Stamenova EK et al. A 3D map of the human genome at kilobase resolution reveals principles of chromatin looping. Cell 2014 Dec 18;159(7):1665-80. PMID: [25497547](https://www.ncbi.nlm.nih.gov/pubmed/25497547)

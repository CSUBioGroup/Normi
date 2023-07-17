# Normi
## Overview of the framework
We introduce Normi, a gene regulatory network (GRN) inference method based on non-redundant mutual information, which comprises three main components: dataset preprocessing, quantitative analysis of regulatory strength using mutual information, and filtering of redundant edges. It first preprocesses the input scRNA-seq by partitioning the whole trajectory into multiple segments, then applies sliding windows and averaging strategy on each segment to obtain smoothed representative cells. Next, compute the distance correlation (dcorr) for each gene pair to determine the optimal time lag. The mixed KSG estimator is utilized to quantify the regulatory strength among genes under optimal time lags. Finally, it adopts mRMR to remove the redundant regulatory edges in the preliminary GRN to obtain the refined network.
![Figure1](https://github.com/CSUBioGroup/Normi/assets/52996140/172a0e50-2b87-4eb5-992c-2e2318fd77c6)



## Dependencies
numpy==1.23.5
<br />pandas==1.5.3
<br />scipy==1.10.0
<br />scikit-learn==1.2.1
<br />dcor==0.6<br />

## Tutorial
We provide tutorial as shown in {demo.ipynb} for introducing the usage of Normi.

## Usage
Normi takes data as input file in csv (genes in columns and cells in rows)(e.g. [ExpressionData.csv](https://github.com/CSUBioGroup/Normi/files/12066810/ExpressionData.csv))
, and it needs an input file containing time information (e.g. [PseudoTime.csv](https://github.com/CSUBioGroup/Normi/files/12066811/PseudoTime.csv)). The output of Normi is a file in csv with edges sorted in descending order by estimated mutual information.

Command to run Normi


```
python main.py --data_file <scGNA-seq path> --time_file <time path> --tf_file <tf list path> --save_dir <output path> --window_size <window size> --slide <sliding length> --n_jobs <n_jobs>
```
<br />We also provide default hyper-parameters in main.py, Using -h option to see more.

# Normi
Normi: A framework to infer gene regulatory network from scRNA-seq

# Requirements
numpy==1.23.5
pandas==1.5.3
scipy==1.10.0

# Usage
Normi takes data as input file in csv (genes in columns and cells in rows), and it needs an input file containing time information. The output of Normi is a file in csv with edges sorted in descending order by estimated mutual information.
We also provide default hyper-parameters in main.py, Using -h option to see more.

Command to run Normi


```
python main.py --data_file <scGNA-seq path> --time_file <time path> --tf_file <tf list path> --save_dir <output path> --window_size <window size> --slide <sliding length> --n_jobs <n_jobs>
```

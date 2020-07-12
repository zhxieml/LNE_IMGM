# LNE IMGM

This repository contains MATLAB implementation of following incremental multi-graph matching methods:

- **IMGM** Tianshu Yu, Junchi Yan, Wei Liu, Baoxin Li, *Incremental Multi-graph Matching via Diversity and Randomness based Graph Clustering*, ECCV 2018.
- **LNE IMGM** Zixuan Chen, Zhihui Xie, Junchi Yan, Yinqiang Zheng, Xiaokang Yang, *Layered Neighborhood Expansion for Incremental Multiple Graph Matching*, ECCV 2020.

## Problem setting

In this codebase inline with our ECCV 2020 paper, we focus on the online setting of graph matching whereby graphs arrive one by one. This setting is nontrivial and calls for efficient mechanism.

## Dataset

Our algorithm is tested with synthetic data and real-world images ([Willow ObjectClass](https://www.di.ens.fr/willow/research/graphlearning/)). Please refer to the paper for details about generating data.

If you want to run experiment on Willow ObjectClass dataset, we provide SIFT-extracted [features](https://drive.google.com/file/d/1Wk0QAK-cey-GkvUN3qHjj9IuZ1AgHESk/view?usp=sharing) and the corresponding [ground-truth](https://drive.google.com/file/d/1i3q42Bv5eJqbwkX2v2PKmYOAE8pk0yyp/view?usp=sharing). For data configuration, please check `load_target_data.m`. 

## Experiment

Run `experiment_*.m` to reproduce results in the paper. Notice that parameters are set by the struct `target.config`.

| Experiment            | Parameters                                           | Code                        | Comments                                                     |
| --------------------- | ---------------------------------------------------- | --------------------------- | ------------------------------------------------------------ |
| Fig. 2, 5             | n_o = 0, c = 1, ϵ = 0.15, ρ = 1, (NA, NB) = (20, 50) | `experiment_online.m`       | Synthetic data                                               |
| Fig. 3                | n_o = 0, c = 1, ϵ = 0.15, ρ = 1, (NA, NB) = (20, 50) | `experiment_rawmat.m`       | Synthetic data                                               |
| Fig. 4, 5             | n_o = 4, β = 0.9, (NA, NB) = (20, 40)                | `experiment_online.m`       | WILLOW data. Use different values of `target.config.maxNumSearch` for Fig. 5 |
| Fig. 7(a)             | n_o = 0, c = 1, ϵ = 0.15, ρ = 1, (NA, NB) = (20, 52) | `experiment_ordering.m`     | Synthetic data                                               |
| Fig. 7(b)             | n_o = 4, β = 0.9, (NA, NB) = (20, 52)                | `experiment_ordering.m`     | WILLOW data (Winebottle)                                     |
| Fig. 7(c)             | n_o = 0, c = 1, ϵ = 0.15, ρ = 1, (NA, NB) = (30, 50) | `experiment_distribution.m` | Synthetic data                                               |
| Fig. 7(d)             | n_o = 4, β = 0.9, (NA, NB) = (30, 50)                | `experiment_distribution.m` | WILLOW data (Car)                                            |
| Fig. 6(a), 6(b), 6(c) | n_o = 4, β = 0.9                                     | `experiment_offline.m`      | WILLOW data                                                  |
| Fig. 6(d), 6(e)       | n_o = 0, c = 1, ϵ = 0.15, ρ = 1                      | `experiment_offline.m`      | Synthetic data                                               |


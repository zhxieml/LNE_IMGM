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

Run `experiment_*.m` to reproduce results in the paper.
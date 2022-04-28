# Paper: A self-adaptive weighted differential evolution approach for large-scale feature selection
## Available  
**Link1**: https://www.sciencedirect.com/science/article/pii/S0950705121008959#!

**Link2**: https://www.wangxubin.site/Paper/SaWDE-KBS.pdf

## Framework
![](https://github.com/wangxb96/SaWDE/blob/master/framework.png)
## Abstract
Recently, many evolutionary computation methods have been developed to solve the feature selection problem. However, the studies focused mainly on small-scale issues, resulting in stagnation issues in local optima and numerical instability when dealing with large-scale feature selection dilemmas. To address these challenges, this paper proposes a novel weighted differential evolution algorithm based on self-adaptive mechanism, named SaWDE, to solve large-scale feature selection. First, a multi-population mechanism is adopted to enhance the diversity of the population. Then, we propose a new self-adaptive mechanism that selects several strategies from a strategy pool to capture the diverse characteristics of the datasets from the historical information. Finally, a weighted model is designed to identify the important features, which enables our model to generate the most suitable feature-selection solution. We demonstrate the effectiveness of our algorithm on twelve large-scale datasets. The performance of SaWDE is superior compared to six non-EC algorithms and six other EC algorithms, on both training and test datasets and on subset size, indicating that our algorithm is a favorable tool to solve the large-scale feature selection problem. Moreover, we have experimented SaWDE with six EC algorithms on twelve higher-dimensional data, which demonstrates that SaWDE is more robust and efficient compared to those state-of-the-art methods. SaWDE source code is available on Github at https://github.com/wangxb96/SaWDE.
## Instructions
- SaWDE.m is the main function
- DataPartition.m is used to randomly partition the original data into training sets and test sets with a ratio of 7 : 3.
- CSGSTest.m is used to test the performance of each strategy.
## Cite
```
@article{wang2022self,
  title={A self-adaptive weighted differential evolution approach for large-scale feature selection},
  author={Wang, Xubin and Wang, Yunhe and Wong, Ka-Chun and Li, Xiangtao},
  journal={Knowledge-Based Systems},
  volume={235},
  pages={107633},
  year={2022},
  publisher={Elsevier}
}
```
## Contact 
wangxb19@mails.jlu.edu.cn -- Jilin University

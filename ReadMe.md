<div align="center">
<h1>A self-adaptive weighted differential evolution approach for large-scale feature selection</h1>

[**Xubin Wang**](https://github.com/wangxb96)<sup>1</sup> · **Yunhe Wang**<sup>2*</sup> · **Ka-Chun Wong**<sup>3</sup> · **Xiangtao Li**<sup>1*</sup>


<sup>1</sup>Jilin University · <sup>2</sup>Hebei University of Technology · <sup>3</sup>City University of Hong Kong

<sup>*</sup>corresponding authors

[**Paper**](https://www.wangxubin.site/Paper/SaWDE-KBS.pdf) · [**Code**](https://github.com/wangxb96/HSNOE)

</div>

## Abstract
Recently, many evolutionary computation methods have been developed to solve the feature selection problem. However, the studies focused mainly on small-scale issues, resulting in stagnation issues in local optima and numerical instability when dealing with large-scale feature selection dilemmas. To address these challenges, this paper proposes a novel weighted differential evolution algorithm based on self-adaptive mechanism, named SaWDE, to solve large-scale feature selection. First, a multi-population mechanism is adopted to enhance the diversity of the population. Then, we propose a new self-adaptive mechanism that selects several strategies from a strategy pool to capture the diverse characteristics of the datasets from the historical information. Finally, a weighted model is designed to identify the important features, which enables our model to generate the most suitable feature-selection solution. We demonstrate the effectiveness of our algorithm on twelve large-scale datasets. The performance of SaWDE is superior compared to six non-EC algorithms and six other EC algorithms, on both training and test datasets and on subset size, indicating that our algorithm is a favorable tool to solve the large-scale feature selection problem. Moreover, we have experimented SaWDE with six EC algorithms on twelve higher-dimensional data, which demonstrates that SaWDE is more robust and efficient compared to those state-of-the-art methods. 

## Framework
![](https://github.com/wangxb96/SaWDE/blob/master/framework.png)
Overview of the proposed SaWDE algorithm. The initialized population $P$ is divided equally into five sub-populations, $SubP_1$, $SubP_2$, $SubP_3$, $SubP_4$, $SubP_5$. Then, each sub-population selects an strategy EnS through the self-adaptive mechanism from the strategy pool. After that, the $P$ is evaluated and updated by the selected EnS. Meanwhile, the strategies EnSs are evaluated for additional rewards and further sub-strategy pool construction. Finally, a weighted model is proposed to assess the importance of each feature and search for the best solution by evaluating these features in a combinatorial way.


## Instructions
- SaWDE.m is the main function
  - The code "MaxFESpre = 1000000; % Function Evaluation Times" on line 63 in SaWDE.m is the number of iterations in the paper. You can modify it to suit your requirements.
  - How to load your own data? You can try the way used in SaWDE.m or the following method:
    ```
      traindata = load(['C:\Users\c\Desktop\SaWDE\train\',p_name]);
      traindata = getfield(traindata, p_name);
      data = traindata;
      feat = data(:,1:end-1); 
      label = data(:,end);
    ```
      
- DataPartition.m is used to randomly partition the original data into training sets and test sets with a ratio of 7 : 3.
- CSGSTest.m is used to test the performance of each strategy.
- CoDE1.m - CoDE10.m are the strategies used in our study. The following files are the related dependency files for our strategies:
  - getStrategyOnRoulette.m
  - gnR1R2.m
  - ldrc.m
  - randFCR.m
  - updateArchive.m
## Dependencies
- This project was developped with **MATLAB 2018b**. Early versions of MATLAB may have incompatibilities.
## Citation
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
wangxb19 at mails.jlu.edu.cn -- Jilin University

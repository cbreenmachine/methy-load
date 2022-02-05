# January 2022
`pandoc --self-contained  summary.md -o test.html`

Some changes:
- Wrap up all post-processing scripts (imputation of missing values; pca; etc.) into one bash script.
- Move to python for as much of the processing as possible
- Put python scripts in `code/python_scripts` for organization
- Take command line arguments
- Integrate with GNU parallel as much as possible

Progress
- Pools 6 and 7 are completely unpacked

Quick To-do
- Md5sum checking on first 100
- Restructure to pools instead of groups for first 100
- QC reports into their respective directories


## January 4, 2022

### Imputation May not be computationally efficient, but imputing methylated reads separately from coverage may be the most sound
- Correlation between methylated and coverage: 0.68
- Correlation between methylated and unmethylated: -0.63
- Correlation between unmethylated and coverage: 

### Two extremes
1. A sample has unique information: perhaps its a SNP or SNV; synonymous or non-synonymous? Who knows.
2. Accounting for sequencing errors
    - There tends to be bimodality. The vast majority of positions have estimates for 100% of sites; but there is some probability (small) that there is a missing estimate. Basically a binomial distribution; but then the marginal increase (the extra counts with lots of missing) are information dense--
    IDEA:
        - handle coverage and methylation estimate separately?? 
        1. Cluster each site: is it part of the binomial distribution? Or is it from the "signature methulation estimate" population?
        2. Separate the two populations.
            - If it's part of the binomial (say above 0.5 probability (this would be a user-specified parameter)) then we impute with local mean; or likelihood based; or whatever method you want (some form of rolling average or perhaps the linear regression version advocated in that paper)
            - If it's part of the signature population; set it aside; it could be part of the information-rich population of CpG traffic lights; could be a (small) island; or could be a really bad sequencing error.

# January 18 2022 (LOAD Update)

0. Scree plots
1. Do (the same) outliers show up in array data?
2. Can we get better estimates? 

## 0. Scree plots

Get about 10% of variance explained from PC1.

<img src="2021-11-pca-batch/figs-2022-01-13/scree-chr22.png" width="70%" height="70%" />

## 1. Outliers

Recall the two outliers we've seen (samples 101 and 126):
<img src="./2021-11-pca-batch/figs-2022-01-11/PC_chr22_labeled_mean_methylation.png" width="70%" height="70%" />

<br/><br/>
And now the array data. 101 is a clear outlier in array data as well.
<img src="./2021-11-pca-batch/figs-array/array_labeled_mean_methylation.png" width="70%" height="70%" />

Unfortunate;y, 126 is not in array data. Likely drop these two from analysis? Probably won't make the final decision until we can see the entire sample of 375-ish.

[Paper that reviews deconvolution methods](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5886462/pdf/ddx275.pdf)
- Nerual network based (seems excessive but don't know details)
- Reference based (used 450K reference for 850K data--since Andy's paper there has been a reference dataset published for 850K [here](https://pubmed.ncbi.nlm.nih.gov/30232931/) but its for neonatal cord blood)

## 2. Improving Estimate of Methylation

## Idea
- fit a glm at each CpG locus to create better methylation estimates.
- account for variability in <em>cell composition</em> specifically.

## Implementation notes
Binomial regression model with logit link function (in other words, logistic regression) <br/>

$log(\frac{p}{1-p}) = X\beta$

where $p$ is probability methylated at a site.

- Adding one pseudocount to make tractable/stable()
- Sampling 10,000 sites from chromosome 22 for speed

- Two models:
    1. Cell composition only (with mean_methylation)
    2. Cell composition (with mean_methylation) <i>and</i> age, sex, BMI
- Used <em>methylation estimate</em> as the response (as opposed to # methylated reads)

## Improvement
- Beta-binomial distribution with arc-sine

## Result

### PC plot from estimated / corrected methylation with only cell composition

<img src="./2022-01-better-estimates/figs/PC-cell-compositions/mean_methylation_labeled.png" width="70%" height="70%"/>

<br/>

### PC plot from estimated / corrected methylation with cell composition + age, sex, bmi
<br/>

<img src="./2022-01-better-estimates/figs/PC-cell-compositions-pheno/mean_methylation_labeled.png" width="70%" height="70%"/>

## Follow-up

Better GLM (ala Hao Wu's DSS) is Beta-Binomial model where $\mathbb{P}[\text{Methylated}] \sim \beta()$

`pandoc --wrap=preserve --self-contained 2022-01-notebook.md -o 2022-01-18.html`


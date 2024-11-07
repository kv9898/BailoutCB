# Abstract: Bail out the money printer? The impact of fiscal indemnity on central bank profitability

Since 2022, central bank losses have been prevalent in advanced economies due to previous quantitative easing and recent inflationary pressures. This paper focuses on the unique case of the United Kingdom, where the government promised in advance to cover any central bank losses arising from quantitative easing. This promise is known as the indemnity. A game-theoretical model is proposed to explain the causes and effects of such indemnity. The model’s predictions about the indemnity’s effect on central bank profitability are empirically examined. Using the novel Dynamic Multilevel Latent Factor Model (DM-LFM), the indemnity is found to have significantly boosted the Bank of England’s profits in the deflationary environment after 2008, but exacerbated its losses under the recent inflationary pressure since 2022. The theoretical model suggests the pronounced effects are due to the Bank of England’s high sensitivity to losses and the UK government’s moderate fiscal liberalism. Therefore, the British experience should not be generalised. Nevertheless, the theoretical and empirical lessons can inform policy-makers about future institutional designs concerning the fiscal-monetary interactions and the public finance-price stability trade-off.

# Instructions for reproducing the analysis
0. (optional) Clean the interest rate data by running ir/clean.R.
1. Run code/1_formal_evaluation.R for the numerical evaluation of the formal model.
2. Run code/2_emp_analysis.R for the empirical analysis of the paper.
3. Run code/3_MCMCdiag.R for the MCMC diagnosis for the Bayesian estimation in the paper.
4. Render manuscript/BailoutCB to reproduce the paper.

# Additional information
 - R (4.4.1) was used
 - Analyses were run on a Huawei Matebook 13 (2021 Intel) using Windows 11 with 16GB of memory.
 - The five 4 scripts (steps 0, 1, 2, 3) and the quarto document (step 4) should take less than 30 minutes to run. The empirical analysis (step 2) takes the longest, about 15 minutes. At the end of running all three scripts, the code, data, outputs, and stored results occupy just under 440MB of hard drive space.
 - Information on required packages is given in each R script.
 - Raw data for the central bank profits and liabilities are deliberately omitted from the repository, as they are acquired from commercial sources. They are available upon request.

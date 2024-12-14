# GWAS Analysis Pipeline in R
This repository contains an R pipeline for conducting Genome-Wide Association Studies (GWAS) using phenotype and genotype data.

## Features
- Outlier detection and data cleaning
- Genotype encoding
- Statistical regression for GWAS
- Visualization (Manhattan and Q-Q plots)
- Identification of significant SNPs



Genome-Wide Association Studies (GWAS)
The collaborator is interested in mapping genetic loci that can affect human Body Mass Index
(BMI). They want to perform a GWAS experiment and they would like you to perform the
analysis. They have collected data for a number of individuals sampled from a population and
they have provided you relative measures of BMI in the file
“Omics99_phenotypes_project1.txt” and SNP genotypes in the file
“Omics99_genotypes_project1.txt”. Note the following:
• There are two columns in the phenotype file, the first contains the name of each individual
in the sample, the second contains the BMI of each individual in the sample. Each row of
the phenotype file lists the name and BMI for an individual in the sample (n rows total.
• In the “genotypes” file, the first row contains the names of N genotypes, each of which
has
been measured for each of the n individuals in the sample. Each of the following rows contains
all of the genotypes measured for a specific individual. For each of these rows, each
consecutive PAIR OF COLUMNS represents a genotype of an individual (i.e., each genotype
name refers to a genotype defined by a pair of columns), where there are a total of N
genotypes for each individual (for the row corresponding to individual i, the 1st and 2nd
column = 1st genotype of individual i, 3rd and 4th columns = 2nd genotype of individual i,
..., rows 2N − 1 and 2N = N genotype for individual i). Each genotype of each individual
is composed of a pair of alleles (e.g., the possible genotypes at the jth genotype for an
individual i are ‘A1A1’, ‘A1A2’, ‘A2A2’). In some cases, instead of an allele ‘A1’ or ‘A2’
there is a number ‘9’ indicating that the measurement of the allele is missing.
1. (a) Plot a histogram of the phenotypes (provide your code). (b) There is one individual
with a phenotype that is clearly considerably larger than the others (an outlier).
Report the name of this individual and remove this individual from the analysis by removing
the phenotype AND the corresponding genotypes of this individual from the data set. (c)
Re- plot the histogram of the phenotypes after removal of this individual. (d) The
phenotypes should now look approximately normal. In no more than one sentence, explain
why it is important that the phenotypes be well modeled by this distribution if we are going
to use the genetic linear regression to model the relationships between genotypes and this
phenotype.
Once you have completed part ‘a’ your dataset now has n − 1 individuals each with
N genotypes. Analyze this filtered data set in Question 2.
2. (a) Write code to identify which individuals have a genotype with a missing allele. Report
the names of individuals with at least one genotype with a missing allele and report the names
of the genotypes for each of these individuals that have missing alleles. The entire list of
these numbers represent genotypes for which at least one individual has missing data. (b) For
your list in part ‘a’ remove these genotypes from the analysis, not just for the individual
where they are missing but FOR EVERY individual in the data set. That is, remove
the PAIRS OF COLUMNS corresponding to your missing genotype list.
3
s
Once you have completed part ‘d’ your dataset now has n − 1 individuals each with N−
the number of genotypes you removed. Analyze this filtered data set for Question 3-5.
3. (a) For EACH genotype, using the formulas provided in class, calculate the MLE(β ˆ ) for
the three β parameters when using the linear regression model yi = βμ + xi,aβa + xi,dβd + si,
with
si ∼ N (0, σ
2
) and plot a histogram for the estimates of each parameter = three histograms
total (provide your code! And make sure you label your plots!). (b) For EACH genotype,
calculate p-values for testing the null hypothesis H0 : βa = 0 ∩ βd = 0 versus the alternative
hypothesis HA : βa ¹ 0 ∪ βd ¹ 0 using the formulas provided in class (i.e. the predicted
value of the phenotype yˆi for an individual i, the SSM, SSE, MSM, MSE, and the F-statistic),
although you may use the function pf( ) with the option ‘lower.tail=FALSE’ to calculate the
p-value from your F-statistic. (c) Plot a histogram of ALL p-values (not -log(p-values)!
just the p-values!). If you ignore a few of the p-values that are extremely small (=highly
significant), describe how the rest of the p-values are distributed (i.e., the probability
distribution do they appear to resemble) and using no more than two sentences, explain what
this indicates about the correct statistical model for the bulk of the genotype-phenotype
relationships in the data. (d) Plot a Manhattan plot using the -log(p-values) and a Q-Q plot
(provide your code!). Explain your Q-Q plot.
4. (a) For a Type 1 error of 0.05, how many genotypes are significant from your analysis in
Question 3? (b) For a Type 1 error of 0.05 / N − 1, how many genotypes are significant
from your analysis in Question 3? (c) Using no more than two sentences, describe why the
second (the lower) of these Type 1 errors is more appropriate for identifying the location
of causal polymorphisms and justify your answer. (d) For the lower of the Type 1 errors,
report the genotypes that you returned (i.e., a list of genotypes from 1 to N− the number
of genotypes you removed in Question 2). Using no more than three sentences, explain
how many causal genotypes these two sets are likely indicating and provide a justification
of your assertion.
5. (a) For the single most significant genotype in your analysis, produce two x-y plots: Xa
vs Y and Xd vs Y (provide your code and make sure you label your plots!). (b) Choose
one of your other genotypes that has a p-value > 0.5 (i.e., choose any one of your
genotypes with a p-value > 0.5) and produce two x-y plots: Xa vs Y and Xd vs Y . (c)
Using no more than three sentences, describe how the plots in parts ‘a’ and ‘b’ differ, an
explanation as to why you believe they differ, and why the plots make sense assuming
that your belief is correct. (d) Using no more than two sentences, explain to your
collaborator why the genotype you have plotted in part ‘a’ is not likely to be the causal
polymorphism for BMI but why it indicates the position in the genome where a causal
polymorphism may be located.

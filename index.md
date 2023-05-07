# Sequoia  <img src="man/figures/sequoia_hexlogo.svg" align="right" height=225 style="float:right; height:200px; vertical-align:middle; margin:0px 20px" /> 


#### Multi-generational pedigree reconstruction from SNP data. Accounting for genotyping errors. Overlapping or discrete generations, with or without inbreeding, any proportion of genotyped parents. No lists of candidate parents needed, just birth years.


<br> 

## What it can do

#### <img src="man/figures/parents.svg" width="100" height="70" /> Parentage assignment 
Candidate parent--offspring pairs are short-listed among all genotyped individuals based on the number of SNPs at which they are opposing homozygotes. Parents are assigned based on the likelihood ratio between the pair being parent--offspring versus the most-likely alternative relationship. The pair can be oriented if their relative age is known, or if there is a complementary co-parent.
<br>


#### <img src="man/figures/siblings.svg" width="100" height="80" /> Sibship clustering
When not all parents were genotyped, clusters of half- and full-siblings are identified, and each assigned a dummy parent. Every dummy individual corresponds to a real-world, non-genotyped individual.
<br>

#### <img src="man/figures/grandparents.svg" width="100" height="70" /> Grandparent assignment
For each cluster of half-siblings, grandparents are assigned where possible, i.e. parents of the sibship's dummy-parent. These grandparents may be dummy individuals themselves, so that parent-offspring links may be established between two non-genotyped individuals. 


<div>
<img src="man/figures/Dummies_for_dummies.svg" width="300" height="300" alt="pedigree example" />
</div>


## Rationale
Pedigree reconstruction with `sequoia()` relies on the likelihood ratios between a focal relationship (e.g. parent-offspring, PO) and a myriad of alternative relationships for that pair (full siblings, aunt-niece, ..., or unrelated U). This method is inspired by the work of E.A. Thompson, such as this excerpt from her paper *'A Paradox of Genealogical Inference'* (1976):

<div>
<img src="man/figures/Thompson_1976_quote1.png" height="150" alt="Thompson 1976 quote"  style="height:150px" />
</div>


In other words, when comparing likelihood ratios LLR(PO/U) *between* candidate parents, by chance some full siblings may have a higher value than the true parent, even in absence of genotyping errors. 

One possible solution is to consider the likelihood of all assignments jointly in an MCMC(-like) approach, but the number of possible pedigree configurations to explore is enormous. 

Comparing for each pair of putative relatives many different relationships makes each assignment rather computationally intensive, but this is offset by various filtering steps (based on e.g. the age difference and Mendelian inconsistencies) and using a 'hill-climbing algorithm' rather than an MCMC. Imagine it as taking slow, careful steps up the mountain, carefully inspecting the direct surroundings before taking a new step, compared to running around the mountainside more or less at random. If there is a fairly clear path to the likelihood top (i.e. moderately high quality SNP data), `sequoia` will usually do fine. If there is not, an MCMC-based approach may be preferable. 



## Package overview
Beside the main function for pedigree reconstruction (`sequoia()`) the R package contains various other functions, amongst others to check agreement between an existing pedigree and the genotype data, or with a newly inferred pedigree. 

<a href="./reference/figures/flowchart.svg">
 <img src="man/figures/flowchart.svg" height="400" alt="flowchart" style="vertical-align:middle;margin:0px 20px; height:400px" />
</a>
<br>
<br>
[print-friendly version](reference/figures/flowchart_no_bg.pdf)

<!--
![](man/figures/flowchart.svg) 
<br>
 --->

<br>
<br>

For detailed information, please see the vignettes (rendered using [bookdown](https://bookdown.org/yihui/bookdown/)):

<a href="./articles/vignette_main/book/background.html">
 <img src="man/figures/sequoia_hexlogo_vignette_main.svg" height="150" alt="main vignette" style="vertical-align:middle;margin:0px 20px; height:150px" />
</a>
<a href="./articles/vignette_age/book/index.html">
 <img src="man/figures/sequoia_hexlogo_vignette_age.svg" height="150" alt="age vignette" style="vertical-align:middle;margin:0px 20px; height:150px"/>
</a>

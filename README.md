Subtype-aware batch correction retains biological signal of integrated  
breast cancer datasets—Supplementary Methods
================
Gil Tomás  
<gil.tomas@igmm.ed.ac.uk>

# 1—Dataset Acquisition

Raw CEL files from ten breast cancer gene expression dataset where
downloaded from [GEO](https://www.ncbi.nlm.nih.gov/geo/) and normalized
with
[fRMA](https://bioconductor.org/packages/release/bioc/html/frma.html),
[RMA](https://www.bioconductor.org/packages/release/bioc/html/oligo.html)
and
[MAS5](https://www.bioconductor.org/packages/release/bioc/html/affy.html).
Demographics of each dataset are shown in **Table
1**.

<table class="table table-striped" style="width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

Table 1—Demographics of datsets used in this study.

</caption>

<thead>

<tr>

<th style="text-align:left;">

GSE

</th>

<th style="text-align:left;">

platform

</th>

<th style="text-align:right;">

nSamples

</th>

<th style="text-align:right;">

fracER+

</th>

<th style="text-align:right;">

fracHER2+

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

GSE5327

</td>

<td style="text-align:left;">

HG-U133A

</td>

<td style="text-align:right;">

58

</td>

<td style="text-align:right;">

0.00

</td>

<td style="text-align:right;">

NA

</td>

</tr>

<tr>

<td style="text-align:left;">

GSE25065

</td>

<td style="text-align:left;">

HG-U133A

</td>

<td style="text-align:right;">

198

</td>

<td style="text-align:right;">

0.62

</td>

<td style="text-align:right;">

0.01

</td>

</tr>

<tr>

<td style="text-align:left;">

GSE2034

</td>

<td style="text-align:left;">

HG-U133A

</td>

<td style="text-align:right;">

286

</td>

<td style="text-align:right;">

0.73

</td>

<td style="text-align:right;">

NA

</td>

</tr>

<tr>

<td style="text-align:left;">

GSE17705

</td>

<td style="text-align:left;">

HG-U133A

</td>

<td style="text-align:right;">

298

</td>

<td style="text-align:right;">

1.00

</td>

<td style="text-align:right;">

NA

</td>

</tr>

<tr>

<td style="text-align:left;">

GSE16446

</td>

<td style="text-align:left;">

HG-U133Plus2

</td>

<td style="text-align:right;">

120

</td>

<td style="text-align:right;">

0.00

</td>

<td style="text-align:right;">

0.33

</td>

</tr>

<tr>

<td style="text-align:left;">

GSE17907

</td>

<td style="text-align:left;">

HG-U133Plus2

</td>

<td style="text-align:right;">

33

</td>

<td style="text-align:right;">

0.47

</td>

<td style="text-align:right;">

1.00

</td>

</tr>

<tr>

<td style="text-align:left;">

GSE21653

</td>

<td style="text-align:left;">

HG-U133Plus2

</td>

<td style="text-align:right;">

266

</td>

<td style="text-align:right;">

0.57

</td>

<td style="text-align:right;">

0.12

</td>

</tr>

<tr>

<td style="text-align:left;">

GSE5460

</td>

<td style="text-align:left;">

HG-U133Plus2

</td>

<td style="text-align:right;">

127

</td>

<td style="text-align:right;">

0.58

</td>

<td style="text-align:right;">

0.24

</td>

</tr>

<tr>

<td style="text-align:left;">

GSE2109

</td>

<td style="text-align:left;">

HG-U133Plus2

</td>

<td style="text-align:right;">

353

</td>

<td style="text-align:right;">

0.65

</td>

<td style="text-align:right;">

0.27

</td>

</tr>

<tr>

<td style="text-align:left;">

GSE23177

</td>

<td style="text-align:left;">

HG-U133Plus2

</td>

<td style="text-align:right;">

116

</td>

<td style="text-align:right;">

1.00

</td>

<td style="text-align:right;">

0.00

</td>

</tr>

</tbody>

</table>

# 2—Microarray Dataset Integration

Microarray dataset integration needs to account for technical,
non-biological, variation across multiple batches of independently
acquired datasets. The goal of conventional batch correction (BC) is to
remove batch effects while retaining biological variation conveyed by
each dataset. Several methods exist to address this task. One popular
approach is
[ComBat](https://www.bu.edu/jlab/wp-assets/ComBat/Abstract.html), which
proposes “parametric and nonparametric empirical Bayes frameworks for
adjusting data for batch effects that is robust to outliers in small
sample sizes and performs comparable to existing methods for large
samples”. However, most BC integration strategies do not account for
imbalanced subtype composition within datasets.

This manuscript introduces a novel procedure for integrating expression
profile datasets of known dissimilar composition, called **subtype-aware
batch correction** (SABC). Molecular subtypes are initially assigned to
each sample in each dataset with a publicly available single sample
predictor (SSP). Then, batch effects are resolved between datasets on a
per-subtype basis. This two-layered approach allows for the biological
specifics captured by the single sample predictor of choice to be
accounted for during the batch correction step, and thus carried over
into the integrated dataset. By ignoring distinct subtype compositions
between datasets, conventional BC incorrectly apprehends biological
variation as technical batch effect, and consequently distorts true
biological signal in the integrated dataset.

Breast cancer is widely understood to be subdivided into five intrinsic
or molecular subtypes, which can be assigned by gene expression profile
single sample predictors. Breast cancer hence provides an ideal case
stud y for the evaluation of SABC. We used conventional BC and SABC to
integrate four breast cancer datasets hybridized onto the Affymetrix
HG-U133a chip and six breast cancer datasets hybridized onto the HG-U133
Plus2 chip (**Table 1**).

To evaluate the performance of each method, we compared the
distributions of expression values for the
[205225\_at](https://genecards.weizmann.ac.il/cgi-bin/geneannot/GA_search.pl?keyword_type=probe_set_id&array=HG-U133&target=genecards&keyword=205225_at)
and the
[216836\_s\_at](https://genecards.weizmann.ac.il/cgi-bin/geneannot/GA_search.pl?keyword_type=probe_set_id&array=HG-U133&target=genecards&keyword=216836_s_at)
probesets in each chip, respectively targeting for the *ESR1* and
*ERBB2* gene transcripts. Estrogen receptor (ESR1) and Erb-B2 Receptor
Tyrosine Kinase 2 receptor (ERBB2) status are strong predictors of
breast cancer prognosis and are traditionally assessed by
immunohistochemestry (IHC). Post dataset integration, the expression
values of both genes should remain in line with the biological signal
conveyed by the independent assessment given by IHC for both proteins.
In addition, we compared the agreement between single sample predictor
class assignments for individual samples prior and post integration. The
batch correction procedure should not interfere with the molecular
subtype identity of each sample, and for that reason higher agreement
rates should be indicative of higher transcriptional fidelity of the
integrated dataset.

To guide the integration process with SABC, we used two SSPs implemented
in the
[Genefu](https://bioconductor.org/packages/release/bioc/html/genefu.html)
Bioconductor package. The first,
[sorlie2003](http://www.pnas.org/content/100/14/8418), is based on 534
diagnostic genes and is a five-subtype classifier; the second,
[desmedts2008](http://clincancerres.aacrjournals.org/content/14/16/5158?ck=nck),
is based on three genes (*ESR1*, *ERBB2* and *AURKA*), and is a
three-subtype classifier. To assess classifier agreement prior and post
integration, and in order to circumvent the redundancy caused by using
the same single sample predictor to integrate and to validate the
integration process, we used the Genefu implementation of the
[PAM50](http://ascopubs.org/doi/abs/10.1200/JCO.2008.18.1370) single
sample predictor, based on 50 genes and yielding a five-subtype
classifier.

    ##    user  system elapsed 
    ##  22.646   8.356  31.081

## 2.1—Distortion of Molecular *ESR1* and *ERBB2* Measurements

We compared the distributions of expression values prior and post
dataset integration for the 205225\_at probeset in 840 breast tumours,
hybridised onto the HG-U133a chip, from four datasets with distinct
fractions of ER+ samples (**Table 1** and **Fig.1**). Regardless of the
normalisation method, BC integration significantly distorts *ESR1*
expression measurements in samples from datasets with extreme fractions
of ER+ samples (GSE5327 and GSE17705), to the point where the two
distributions no longer can tell the difference between ER– and ER+
samples (**Fig.2**, second column of panels). In both cases, SABC
integration, whether driven by a 3-subtype SSP (SABC3) or a 5-subtype
SSP (SABC5), succeeds in retaining the biological signal conveyed by the
IHC status in datasets with extreme compositions (**Fig.2**, third and
fourth columns of panels). Further comparison of ESR1 transcript
abundance in these two datasets prior and post integration reveals that
expression values depart significantly from original measurements when
datasets are integrated with BC, yet are preserved from extreme
distortion by SABC integration
(**Fig.3**).

<img src="/Volumes/igmm/sims-lab/Gil/composition-draft/composition/readme_files/figure-gfm/er-a-boxplots-1.png" title="Figure 1---205225_at probeset measurements hybridised onto the HG-U133a chip broken by dataset prior and post integration.  Integration was done using standard batch correction (BC, with ComBat) and subtype-aware batch correction (SABC3, driven by a three-subtype SSP; and SABC5, driven by a five-subtype SSP).  The distributions are further split by ER status, independently assessed by IHC on fresh frozen specimens.  Raw data was normalised with FRMA, RMA and MAS5." alt="Figure 1---205225_at probeset measurements hybridised onto the HG-U133a chip broken by dataset prior and post integration.  Integration was done using standard batch correction (BC, with ComBat) and subtype-aware batch correction (SABC3, driven by a three-subtype SSP; and SABC5, driven by a five-subtype SSP).  The distributions are further split by ER status, independently assessed by IHC on fresh frozen specimens.  Raw data was normalised with FRMA, RMA and MAS5." style="display: block; margin: auto auto auto 0;" />

<img src="/Volumes/igmm/sims-lab/Gil/composition-draft/composition/readme_files/figure-gfm/er-a-extreme-datasets-by-er-status-1.png" title="Figure 2---Distributions of probeset 205225_at measurements from datasets GSE5327 (n=58, all ER--) and GSE17705 (n=298, all ER+), taken from **Fig.1**, are compared side by side. See **Fig.1** for details." alt="Figure 2---Distributions of probeset 205225_at measurements from datasets GSE5327 (n=58, all ER--) and GSE17705 (n=298, all ER+), taken from **Fig.1**, are compared side by side. See **Fig.1** for details." style="display: block; margin: auto auto auto 0;" />

<img src="/Volumes/igmm/sims-lab/Gil/composition-draft/composition/readme_files/figure-gfm/er-a-extreme-datasets-by-integration-1.png" title="Figure 3---Distributions of probeset 205225_at measurements from datasets GSE5327 (n=58, all ER--) and GSE17705 (n=298, all ER+), taken from **Fig.1**, are each compared prior and post dataset integration.  See **Fig.1** for details." alt="Figure 3---Distributions of probeset 205225_at measurements from datasets GSE5327 (n=58, all ER--) and GSE17705 (n=298, all ER+), taken from **Fig.1**, are each compared prior and post dataset integration.  See **Fig.1** for details." style="display: block; margin: auto auto auto 0;" />

Although less pronounced, a similar trend is observed when comparing the
expression values prior and post- dataset integration for the
216836\_s\_at probeset, in 1015 breast tumours hybridised onto the
HG-U133Plus2 chip, from six datasets with distinct proportions of HER2+
samples (**Table 1** and **Fig.4**). Regardless of the normalisation
protocol, The ERBB2 probeset distributions clearly reflect the IHC HER2
receptor status in the two datasets in our analysis with extreme HER2
compositions (GSE23177, all HER2–; and GSE17907, all HER2+). When the
five datasets are integrated with BC, this biological signal is erased,
yet preserved (albeit to a lesser extent than in the original datasets),
when integrated with SABC3 and SABC5
(**Fig.4**).

<img src="/Volumes/igmm/sims-lab/Gil/composition-draft/composition/readme_files/figure-gfm/her2-p2-extreme-datasets-by-her-status-1.png" title="Figure 4---216836_s_at probeset measurements hybridised onto the HG-U133Plus2 chip are shown for datasets GSE23177 (n=116, all HER2--) and GSE17907 (n=37, all HER2+), in the leftmost column.  These two datasets were integrated with GSE16446, GSE21653, GSE5460 and GSE2109 (**Table 1**), with BC, SABC3 and SABC5.  ERBB2 expression values for the samples in the two datasets with extreme HER2+ compositions are then shown post-integration with each of these methods, broken by normalisation procedure (cf. **Fig.1** for more details)." alt="Figure 4---216836_s_at probeset measurements hybridised onto the HG-U133Plus2 chip are shown for datasets GSE23177 (n=116, all HER2--) and GSE17907 (n=37, all HER2+), in the leftmost column.  These two datasets were integrated with GSE16446, GSE21653, GSE5460 and GSE2109 (**Table 1**), with BC, SABC3 and SABC5.  ERBB2 expression values for the samples in the two datasets with extreme HER2+ compositions are then shown post-integration with each of these methods, broken by normalisation procedure (cf. **Fig.1** for more details)." style="display: block; margin: auto auto auto 0;" />

The different degree to which subtype-aware batch correction
successfully retains biological signal from these two molecular
correlates of breast cancer biology, in datasets with extreme
compositions, could be explained by how well the classifiers used to
drive integration capture the underlying biology of the ER and HER2
receptors. Because the clinical HER2 breast cancer phenotype is the
result of a gene amplification, it is possible that gene expression
classifiers are less apt to model binary gene expression distributions
(HER2) rather than continuous ones (ER). In addition, the top level
split highlighted by most molecular characterisations of breast tumours
is lead by ER status, consigning most SSP derived from molecular data to
be particularly sensitive to biological signal conveyed by this marker.

## 2.2—Single Sample Predictor Agreement

Single sample predictor assignments for the Genefu implementation of the
PAM50 breast cancer classifier were computed for each of the 1015
samples in GSE2109, GSE21653, GSE5460, GSE16446, GSE23177 and GSE17907,
after FRMA normalisation (**Table 1**). We then integrated the six
datasets with BC, SABC3 and SABC5, and computed each sample’s PAM50
subtype post-integration. Comparisons of subtype assignments prior- and
post-integration can be seen on
**Fig.5**.

<img src="/Volumes/igmm/sims-lab/Gil/composition-draft/composition/readme_files/figure-gfm/pam50-concordance-plot-1.png" title="Figure 5---PAM50 subtype assignments for 1015 samples from six datasets (**Table 1**) hybridised on the HG-U133Plus2 chip and normalised with FRMA.  Subtype assignements were computed prior to dataset integration (first row) and post integration with BC, SABC3 and SABC5 (subsequent rows)." alt="Figure 5---PAM50 subtype assignments for 1015 samples from six datasets (**Table 1**) hybridised on the HG-U133Plus2 chip and normalised with FRMA.  Subtype assignements were computed prior to dataset integration (first row) and post integration with BC, SABC3 and SABC5 (subsequent rows)." style="display: block; margin: auto auto auto 0;" />

<table class="table table-striped" style="width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

Table 2—Interrater agreement of PAM50 subtype assessment prior- and
post-dataset integration (see text for details).

</caption>

<thead>

<tr>

<th style="text-align:left;">

GSE

</th>

<th style="text-align:right;">

nSamples

</th>

<th style="text-align:right;">

fracER+

</th>

<th style="text-align:right;">

fracHER2+

</th>

<th style="text-align:right;">

Cohen’s Kappa: BC

</th>

<th style="text-align:right;">

Cohen’s Kappa: SABC3

</th>

<th style="text-align:right;">

Cohen’s Kappa: SABC5

</th>

<th style="text-align:right;">

Cohen’s Kappa: HRMN5

</th>

<th style="text-align:right;">

Cohen’s Kappa: BCCOV5

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

all

</td>

<td style="text-align:right;">

1015

</td>

<td style="text-align:right;">

NA

</td>

<td style="text-align:right;">

NA

</td>

<td style="text-align:right;">

0.72

</td>

<td style="text-align:right;">

0.78

</td>

<td style="text-align:right;">

0.80

</td>

<td style="text-align:right;">

0.72

</td>

<td style="text-align:right;">

0.83

</td>

</tr>

<tr>

<td style="text-align:left;">

GSE2109

</td>

<td style="text-align:right;">

353

</td>

<td style="text-align:right;">

0.65

</td>

<td style="text-align:right;">

0.27

</td>

<td style="text-align:right;">

0.95

</td>

<td style="text-align:right;">

0.85

</td>

<td style="text-align:right;">

0.93

</td>

<td style="text-align:right;">

0.92

</td>

<td style="text-align:right;">

0.96

</td>

</tr>

<tr>

<td style="text-align:left;">

GSE21653

</td>

<td style="text-align:right;">

266

</td>

<td style="text-align:right;">

0.57

</td>

<td style="text-align:right;">

0.12

</td>

<td style="text-align:right;">

0.80

</td>

<td style="text-align:right;">

0.90

</td>

<td style="text-align:right;">

0.85

</td>

<td style="text-align:right;">

0.80

</td>

<td style="text-align:right;">

0.85

</td>

</tr>

<tr>

<td style="text-align:left;">

GSE5460

</td>

<td style="text-align:right;">

127

</td>

<td style="text-align:right;">

0.58

</td>

<td style="text-align:right;">

0.24

</td>

<td style="text-align:right;">

0.72

</td>

<td style="text-align:right;">

0.75

</td>

<td style="text-align:right;">

0.76

</td>

<td style="text-align:right;">

0.74

</td>

<td style="text-align:right;">

0.72

</td>

</tr>

<tr>

<td style="text-align:left;">

GSE16446

</td>

<td style="text-align:right;">

120

</td>

<td style="text-align:right;">

0.00

</td>

<td style="text-align:right;">

0.33

</td>

<td style="text-align:right;">

0.21

</td>

<td style="text-align:right;">

0.63

</td>

<td style="text-align:right;">

0.50

</td>

<td style="text-align:right;">

0.20

</td>

<td style="text-align:right;">

0.70

</td>

</tr>

<tr>

<td style="text-align:left;">

GSE23177

</td>

<td style="text-align:right;">

116

</td>

<td style="text-align:right;">

1.00

</td>

<td style="text-align:right;">

0.00

</td>

<td style="text-align:right;">

0.30

</td>

<td style="text-align:right;">

0.25

</td>

<td style="text-align:right;">

0.40

</td>

<td style="text-align:right;">

0.32

</td>

<td style="text-align:right;">

0.43

</td>

</tr>

<tr>

<td style="text-align:left;">

GSE17907

</td>

<td style="text-align:right;">

33

</td>

<td style="text-align:right;">

0.47

</td>

<td style="text-align:right;">

1.00

</td>

<td style="text-align:right;">

0.58

</td>

<td style="text-align:right;">

0.79

</td>

<td style="text-align:right;">

0.61

</td>

<td style="text-align:right;">

0.65

</td>

<td style="text-align:right;">

0.69

</td>

</tr>

</tbody>

</table>

With the exception of the largest dataset, GSE2109, PAM50 interrater
agreement was always higher with SABC integration than with conventional
integration. Incidentally, the datasets that showed lesser subtype
agreement post BC integration are the ones with most extreme ER and HER2
compositions (GSE16446, GSE23177 and GSE17907). For these datasets, SABC
integration was able to increase subtype assignment agreement (with the
exception of SABC3 for GSE23177).

Because sample size can further confound the effect of dataset
composition towards retention of biological signal post-integration, we
compared the PAM50 interrate agreement prior and post-integration for
datasets GSE2109, GSE21653, GSE5460, GSE16446 and GSE23177, sampling in
each dataset, 116 samples with the same IHC ER+ proportions as the
original datasets. Comparisons of subtype assignements prior and
post-integration in this experiment can be seen in **Fig.6** and **Table
3**.

## 2.3—Comparison with other Methods

SABC is not the first method to account for external biological
variables or clinical covariates during dataset integration. ComBat
offers the possibility of using clinical variables as covariates during
dataset integration.
[Harman](https://www.ncbi.nlm.nih.gov/pubmed/27585881 "Risk-conscious correction of batch effects: maximising information extraction from high-throughput genomic datasets"),
a more recent integration procedure using a PCA and a constrained
optimisation technique, also calls for factor coding for an
“experimental grouping variable” to drive integration. We ran ComBat
and Harman to integrate the 1015 breast tumours hybridised onto the the
HG-U133Plus2 chip, using the the output of the 5-class sorlie2003 SSP as
covariate and experimental grouping variable
respectively.

<img src="/Volumes/igmm/sims-lab/Gil/composition-draft/composition/readme_files/figure-gfm/er-a-boxplots-other-1.png" title="Figure X---205225_at probeset measurements hybridised onto the HG-U133a chip broken by dataset prior and post integration.  Integration was done using subtype-aware batch correction (SABC5, driven by a five-subtype SSP), ComBat (BCCOV5, using a five-subtype SSP as covariate) and Harman (HRMN5, using a five-subtype SSP as a grouping variable).  The distributions are further split by ER status, independently assessed by IHC on fresh frozen specimens.  Raw data was normalised with FRMA, RMA and MAS5." alt="Figure X---205225_at probeset measurements hybridised onto the HG-U133a chip broken by dataset prior and post integration.  Integration was done using subtype-aware batch correction (SABC5, driven by a five-subtype SSP), ComBat (BCCOV5, using a five-subtype SSP as covariate) and Harman (HRMN5, using a five-subtype SSP as a grouping variable).  The distributions are further split by ER status, independently assessed by IHC on fresh frozen specimens.  Raw data was normalised with FRMA, RMA and MAS5." style="display: block; margin: auto auto auto 0;" />

<img src="/Volumes/igmm/sims-lab/Gil/composition-draft/composition/readme_files/figure-gfm/er-a-extreme-datasets-by-er-status-other-1.png" title="Figure Y---Distributions of probeset 205225_at measurements from datasets GSE5327 (n=58, all ER--) and GSE17705 (n=298, all ER+), taken from **Fig.X**, are compared side by side. See **Fig.X** for details." alt="Figure Y---Distributions of probeset 205225_at measurements from datasets GSE5327 (n=58, all ER--) and GSE17705 (n=298, all ER+), taken from **Fig.X**, are compared side by side. See **Fig.X** for details." style="display: block; margin: auto auto auto 0;" />

<img src="/Volumes/igmm/sims-lab/Gil/composition-draft/composition/readme_files/figure-gfm/pam50-concordance-plot-other-1.png" title="Figure Z---PAM50 subtype assignments for 1015 samples from six datasets (**Table 1**) hybridised on the HG-U133Plus2 chip and normalised with FRMA.  Subtype assignements were computed prior to dataset integration (first row) and post integration with BC, SABC3 and SABC5 (subsequent rows)." alt="Figure Z---PAM50 subtype assignments for 1015 samples from six datasets (**Table 1**) hybridised on the HG-U133Plus2 chip and normalised with FRMA.  Subtype assignements were computed prior to dataset integration (first row) and post integration with BC, SABC3 and SABC5 (subsequent rows)." style="display: block; margin: auto auto auto 0;" />

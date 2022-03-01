
# Background

From University of Illinois Sequencing Center (via email exchange with Reid and Alvaro from December 6, 2021):

> There are two differnent spike ins. All libraries up to date, including the early tests, have been sequenced in teh same way. I

> Internal non-methylated spike-in. This is the one you told us that you were getting something around 0.1% of the reads. This is in each and every one of the libraries. See explanation below.

> The runs are spiked with ~ 5% of PhiX. However, the PhiX library is indexed and the barcodes are not included in the demultiplexing, so there should be no PhiX reads in the sequences you get.

From NEB support:

> We do recommend that each sample receive control DNA spike-ins. The controls are important because they assess the conversion efficiency in each sample. Even though you might be tempted to add the control only to 1 sample because the rest of the samples might have been processed from the same master mixes and reagents, I will caution that it will be difficult to assess for sample-to-sample variability based on sample quality or variability introduced by the user. In terms of what levels of conversion to expect from each control, please take a look at the FAQ from the product page:
> https://www.neb.com/faqs/2019/06/04/what-levels-of-conversion-are-typical-with-the-control-dnas-supplied-in-the-em-seq-kit
 
> For data processing we do have some resources available, please see below:
> This paper/preprint (refers to our data analysis repository).
> https://www.biorxiv.org/content/10.1101/2019.12.20.884692v2

In the Github repository the most relevant file is linked below:
https://github.com/nebiolabs/EM-seq/blob/master/em-seq.nf

For libraries sequenced to a depth of 2–4 M paired end reads, approximately 5,000 x 76 base paired end reads of unmethylated lambda and 500 x 76 base paired end reads of CpG methylated pUC19 are needed to give enough reads for accurate conversion estimates. If these same libraries are sequenced to a higher depth of 200–400 M reads per library, then the number of reads associated with the controls would be in vast excess, 500,000 for unmethylated lambda and 50,000 for pUC19.

The two controls are added to the libraries at different ng concentrations in an attempt to approximately give similar coverage depths. A minimum of 10X coverage is ok, but 16-20 x is better and does not require many more reads. 


# Problem Setting

The `gemBS` pipeline 

# Data

We will want to re-run mapping with five samples while specifying the control sequence. Pool01-group01 has a potential outlier (101 stands out in PCA plots)

# gemBS configuration
http://statgen.cnag.cat/gemBS/UserGuide/_build/html/pipelineCall.html

> For this to work (a) there have to be conversion control sequences in the sequencing library, (b) these have to have been included in the reference sequence and (c) the names of the control sequences have to had been supplied to gemBS when the mapping was performed. If one or both of the values can not be estimated then gemBS will use defaults (0.01 and 0.05) for the missing value(s).

Action item: since the PhiX are not included in the de-multiplexing we have the (under?)-conversion but not over-conversion sequence. Thus, we cannot used the Bayes' network to estimate these numbers. Stick to defaults: 0.05, 0.01
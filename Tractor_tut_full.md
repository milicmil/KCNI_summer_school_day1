---
title: Local Ancestry Deconvolution, Extract and Plotting
filename: Tractor_tut_full.md
---

&nbsp;  
&nbsp;  

Welcome to the first tutorial of the in person summer school 2022. Today we will use Tractor tool in order to include admixed individuals in association studies by leveraging local ancestry information. We will then run a GWAS on this mixed dataset and visualise the results.

# Preliminary steps. Logging in and navigating to the work environment.

We assume that you have set up your Sci-net account and that you have applied and been approved to use the **NIAGARA** cluster.


#### Logging into Niagara sci-net cluster via Jupyterhub

You can log into the Niagara cluster via Jupyrterhub.

Open a chrome browser and go to https://jupyter.scinet.utoronto.ca/. Use your Compute canada credential to log in. Click on the **JupyterLab** button in the top right and open a **Terminal** in the Other options.

&nbsp;  
&nbsp;

#### Navigating to the Sci-net work environment

You should now be located in the **$LOGIN** node. This is where we will not be conducting the tutorial. We will conduct the tutorial in the **$SCRATCH** node. To go to the **$SCRATCH** node type in ```cd $SCRATCH```

the absolute file path of the **$SCRATCH** space is  ```/gpfs/fs0/scratch/sponsor/SPONSOR_NAME/YOUR_SCINET USERNAME/```

#### Obtaining tutorial dataset, Option 1, download from github

The example cohort dataset we are going to use here consists of chromosome 22 for 61 African American individuals from the [Thousand Genome Project](https://www.internationalgenome.org/). These individuals are two-way admixed with components from continental Europe (EUR) and continental Africa (AFR). We simulated phenotypes for these individuals for use in the GWAS.

**DOWNLOAD URL**
The example dataset is found on the following link [example dataset](https://github.com/Atkinson-Lab/Tractor-tutorial/blob/main/tutorial-data.zip) that you may analyze to follow along with this tutorial.

save and unzip the data in a folder in your scratch space called **tutorial_data**.

#### Obtaining tutorial dataset, Option 2, using a Scinet local copy of the data

If you are having trouble downloading or unzipping the data, we have a local copy of the dataset on Sci-net.

While in the main scratch folder for your username, the command below will copy the data to your scratch space.

```
cp -r /gpfs/fs1/home/d/dfelsky/milicmil/tutorial_data/ /gpfs/fs0/scratch/SPONSOR_FIRST_INITIAL/SPONSOR_NAME/YOUR_SCINET_USERNAME/
```

We also provide a haplotype reference panel that can be downloaded from  [Shapeit](https://mathgen.stats.ox.ac.uk/impute/data_download_1000G_phase1_integrated_SHAPEIT2_16-06-14.html), which will be used for phasing, along with phased reference VCF files for local ancestry inference. Here is a complete list of the files:

&nbsp;  
&nbsp;


```
tutorial_data
│
└───ADMIX_COHORT
│       |ASW.unphased.vcf.gz
│   
└───AFR_EUR_REF
|       |YRI_GBR.phased.vcf.gz
|       |YRI_GBR.tsv
|       |chr22.genetic.map.modified.txt
|
└───HAP_REF    
|       |chr22.hap.gz
|       |chr22.legend.gz
|       |ALL.sample
|       |chr22.genetic.map.txt
└───PHENO   
        |Phe.txt
        |Phe2.txt
```

&nbsp;  
&nbsp;  

# Software required for the tutorial
Now that the data is in your scratch space, You will need to load the following modules and their dependencies.

[Shapeit2: v2.r837](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html)  

[Rfmix version 2](https://github.com/slowkoni/rfmix/blob/master/MANUAL.md)  

[Tractor (requires the following python packages: numpy, statsmodel, pandas)](https://github.com/Atkinson-Lab/Tractor) https://github.com/Atkinson-Lab/Tractor


#### Loading Python module with numpy and pandas installed

The commands below tell you what modules on Sci-net you should load.
Load python 3.8.5 (It has numpy and pandas already installed)
```
module load python/3.8.5
```
#### installing statsmodel

The command below will install **statsmodels** pythjon package (user tag is applied automatically).

```
pip install statsmodels
```
#### Loading Shapeit Rfmix and copyting Tractor script folder

We have made a custom environment where Shapeit and Rfmix  are installed.
For it to work we need to change the pointer for the two modules.

```
module use /gpfs/fs1/home/d/dfelsky/dfelsky/modules/
```
Once the command is finished, run the following two lines of code in order to load Shapeit and Rfmix.
```
module load shapeit
```
```
module load rfmix
```

To get Tractor, copy it to your folder from the following location.

```
cp -r /gpfs/fs1/home/d/dfelsky/milicmil/Tractor/ /gpfs/fs0/scratch/*/SPONSOR_NAME/YOUR_SCINET_USERNAME/
```

&nbsp;  
&nbsp;

#### Starting debug nodes

Some of the code blocks in the tutorial require more resources than what the default login node would supply. While we can submit our code as a batch job, depending on the queue at the time the batch job may have to wait until it is executed.

In this tutorial we will use the debug node.

This is a specific node that gives you more computing resources without entering the batch job queue. **NOTE: The debug node lasts 60 minutes** after which you will have to repeat the `debugjob` command in order to reinitailise it.

loads the ddt module that allows for debug node initialization.
```
module load ddt
```

This command starts the debug node which gives additional computing resources to complete the tutorial commands.
```
debugjob
```

&nbsp;  
&nbsp;


# Part 1. Statistical Phasing

#### Phasing and Local Ancestry inference
Before running the Tractor GWAS method, data will need to be phased and have their local ancestry inferred. Today we will use the provided example code for running the full pipeline including these preliminary steps, assuming QC'ed but unphased cohort data.



#### Aligning common variants between query vcf and reference panel

The goal of this step is to find the common variants of the haplotype reference panel and our admixed population.

**All commands should be run in the tutorial_data folder in your SCRATCH space.**


**Aligning common variants between query vcf and reference panel (FROM TRACTOR TUTORIAL)**
```       
shapeit -check \
        --input-vcf ADMIX_COHORT/ASW.unphased.vcf.gz \
        --input-map HAP_REF/chr22.genetic.map.txt \
        --input-ref HAP_REF/chr22.hap.gz HAP_REF/chr22.legend.gz HAP_REF/ALL.sample \
        --output-log alignments
```

The program will print some information and also throw some error message. You may pay attention to some of these messages:

```       
* 61 samples included
* 217153 SNPs included


# if reference panel sites are different from query vcf file, shapeit2 might print:
* #Missing sites in reference panel = 123456
* #Misaligned sites between panels = 123456
* #Multiple alignments between panels = 123456
```

&nbsp;  
&nbsp;  

#### Phasing (this is the most computational expensive step).

We will perform the actual phasing in this step. Notice we should pass argument `--exclude-snp alignments.snp.strand.exclude` to the program if you encounter errors with your own dataset, so that conflict variants are excluded. In this command, we direct Shapeit to our cohort data (`--input-vcf`) which is in a directory with the input data, the reference panel and recombination map files, and tell it where to place the output (`-O`).

After running this command, you should find that two file (`ASW.phased.haps` & `ASW.phased.sample`) have been created in the `ADMIX_COHORT` output directory.

**Phasing command(this is the most computational expensive step)**
```       
shapeit  --input-vcf ADMIX_COHORT/ASW.unphased.vcf.gz \
      --input-map HAP_REF/chr22.genetic.map.txt \
      --input-ref HAP_REF/chr22.hap.gz HAP_REF/chr22.legend.gz HAP_REF/ALL.sample \
      -O ADMIX_COHORT/ASW.phased

      # add this if shapeit throw error message:   --exclude-snp alignments.snp.strand.exclude
```    
&nbsp;  
&nbsp;  

In case we are pressed for time, one can copy the output from the location below. The `ASW.phased` file needs to be copied to **ADMIX_COHORT/** folder.

```
cp /gpfs/fs1/home/d/dfelsky/milicmil/phase/ASW.phased.haps /gpfs/fs0/scratch/*/SPONSOR_NAME/YOUR_SCINET_USERNAME/tutorial_data/ADMIX_COHORT/

cp /gpfs/fs1/home/d/dfelsky/milicmil/phase/ASW.phased.sample /gpfs/fs0/scratch/*/SPONSOR_NAME/YOUR_SCINET_USERNAME/tutorial_data/ADMIX_COHORT/
```
&nbsp;  
&nbsp;       


#### Converting  shapeit output to vcf format.

Shapeit provides a convenient function to convert from its `haps`/`sample` file format to `vcf` format. Notice in the new vcf file we just created, slashes have changed to the pipe or vertical bar to seperate genotype.

**Converting shapeit output to vcf format**
```       
shapeit -convert \
        --input-haps ADMIX_COHORT/ASW.phased\
        --output-vcf ADMIX_COHORT/ASW.phased.vcf
```   

#### Local Ancestry Inference

Here we use Rfmix to perform local ancestry inference, given its high accuracy in multi-way admixed populations. RFmix takes the following argument and will produce 4 files `ASW.deconvoluted.fb.tsv`, `ASW.deconvoluted.msp.tsv`, `ASW.deconvoluted.rfmix.Q`, `ASW.deconvoluted.sis.tsv`.



**Local Ancestry Inference**
```
rfmix -f ADMIX_COHORT/ASW.phased.vcf\
        -r AFR_EUR_REF/YRI_GBR.phased.vcf.gz \
        -m AFR_EUR_REF/YRI_GBR.tsv \
        -g AFR_EUR_REF/chr22.genetic.map.modified.txt \
        -o ADMIX_COHORT/ASW.deconvoluted \
        --chromosome=22
```

---
Now you have your local ancestry calls and are ready for Tractor! In the next step we extract risk allele information from our data.

&nbsp;  
&nbsp;

# Part 2. Extract Tracts
**https://github.com/Atkinson-Lab/Tractor-tutorial/blob/main/Extract.md**


We provide a script that can simultaneously extract risk allele information and local ancestry information. Simply type the following command in terminal:
**Extract Tracts**
```
python Tractor/ExtractTracts.py \
      --msp ADMIX_COHORT/ASW.deconvoluted \
      --vcf ADMIX_COHORT/ASW.phased \
      --num-ancs 2
```

6 new files will be generated in the ADMIX_COHORT folder:
```
ASW.phased.anc0.dosage.txt
ASW.phased.anc0.hapcount.txt
ASW.phased.anc0.vcf
ASW.phased.anc1.dosage.txt
ASW.phased.anc1.hapcount.txt
ASW.phased.anc1.vcf
```

&nbsp;  
&nbsp;

# Part 3. Local Tractor Run
**https://github.com/Atkinson-Lab/Tractor-tutorial/blob/main/Local.md**


To perform linear regression on a (simulated) continuous phenotype in our admixed cohort, simply type the following command in your terminal:

**Using Tractor to perform Linear regression**
```
python RunTractor.py --hapdose ADMIX_COHORT/ASW.phased --phe PHENO/Phe.txt --method linear --out SumStats.tsv
```
Output is **SumStats.tsv** in the **tutorial_data** folder.

&nbsp;  
&nbsp;

# Part 4. Visualise our GWAS results using R in a Jupyter Notebook**

In Jupyterlab, you can open https://jupyter.scinet.utoronto.ca/ a new R console by pressing the + button in the top left. Open a notebook of **R(4.0.3)**.

Note that **qqman** package is not present in the main R library and you will need to provide a folder in your scratch space to install the package outside of the main R library(users cannot install packages in the main Sci-net R library).

in your main scratch space make a directory for the R packages ```mkdir r_packages```

```
install.packages("qqman", repos="http://cran.r-project.org", lib="/gpfs/fs0/scratch/SPONSOR_FIRST_INITIAL/SPONSOR_NAME/YOUR_SCINET_USERNAME/r_packages")
```

To run the plot, use the code below in R.

```
setwd('/gpfs/fs0/scratch/SPONSOR_FIRST_INITIAL/SPONSOR_NAME/YOUR_SCINET_USERNAME/tutorial_data')
library("qqman", lib="/gpfs/fs0/scratch/SPONSOR_FIRST_INITIAL/SPONSOR_NAME/YOUR_SCINET_USERNAME/r_packages")
sumstats = read.csv("SumStats.tsv", sep = "\t")

par(mfrow=c(1,2))
manhattan(sumstats[!is.na(sumstats$ANC0P),], chr="CHROM", bp="POS", snp="ID", p="ANC0P",
          xlim = c(min(sumstats$POS), max(sumstats$POS)), ylim = c(0,15), main = "AFR")
manhattan(sumstats[!is.na(sumstats$ANC1P),], chr="CHROM", bp="POS", snp="ID", p="ANC1P",
          xlim = c(min(sumstats$POS), max(sumstats$POS)), ylim = c(0,15), main = "EUR")
```

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

#### Obtaining tutorial dataset, Option 2, using a google drive link

in order to download the data using google drive, we will use the python package **gdown**.
First you will need to load python using this command `module load python/3.8.5` and then install gdown using `pip install gdown`.

When gdown is installed, make a python file using nano `nano download_data.py` and copy the following code below into it.

replace **url** with this link in QUOTATIONS: https://drive.google.com/drive/folders/1z9NAOs1E37xHRjAfUbXUvfU4L7YMuXbf?usp=sharing

```
import gdown
gdown.download_folder(url, quiet=True)
```
press `ctrl+x` to save and exit the file. Run the file `python download_data.py`.

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
python 3.8.5 (It has numpy and pandas already installed) (This should already be done during the tutorial_data copy step)
```
module load python/3.8.5
```
#### installing statsmodel

The command below will install **statsmodels** pythjon package (user tag is applied automatically).

```
pip install statsmodels
```
#### Loading Shapeit Rfmix and copying Tractor script folder

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

If you type in `module list` you should see Python, shapeit and rfmix loaded in your environment.


To get Tractor, use the same **gdown** script you used to download the **tutorial_data** except this time us this url. `https://drive.google.com/drive/folders/1L9FdQZcDbY7e2J8aFsgTObDoYnMASI9-?usp=sharing`

**NOTE: Run the script from the tutorial_data folder. The Tractor folder should be a subdirectory of the tutorial_data folder.**

&nbsp;  
&nbsp;

#### Starting debug nodes

Some of the code blocks in the tutorial require more resources than what the default login node would supply. While we can submit our code as a batch job, depending on the queue at the time the batch job may have to wait until it is executed.

In this tutorial we will use the debug node.

This is a specific node that gives you more computing resources without entering the batch job queue. **NOTE: The debug node lasts 60 minutes** after which you will have to repeat the `debugjob` command in order to reinitailise it.


This command starts the debug node which gives additional computing resources to complete the tutorial commands.
```
debugjob
```

You may have to wait a little bit before the resources are allocated to your debug node. You can see the sample output below.

```
debugjob: Requesting 1 node(s) with 40 task(s) for 60 minutes and 0 seconds
SALLOC: Pending job allocation 7634116
SALLOC: job 7634116 queued and waiting for resources
SALLOC: job 7634116 has been allocated resources
SALLOC: Granted job allocation 7634116
SALLOC: Waiting for resource configuration
SALLOC: Nodes nia0061 are ready for job

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
The job should run about 150 seconds in debug node.

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
**NOTE: Chances are, there wont be enough computational time in debug node (60 min) to complete this step**. In that case you should either copy results or run the batch script below.

&nbsp;  
&nbsp;  

In case we are pressed for time, one can copy the output using the gdown script below with these URLs:

download script for "ASW.phased.haps"

https://drive.google.com/file/d/1fPa50ZhY9PY_gRTBrg1CVAE4MOPAHCLO/view?usp=sharing
https://drive.google.com/file/d/1fPa50ZhY9PY_gRTBrg1CVAE4MOPAHCLO/view?usp=sharing
```
import gdown
url = "https://drive.google.com/file/d/1fPa50ZhY9PY_gRTBrg1CVAE4MOPAHCLO/view?usp=sharing"
output = "ASW.phased.haps"
gdown.download(url, output, quiet=True, fuzzy=True)

```

download script for "ASW.phased.sample"
```
import gdown
url = "https://drive.google.com/file/d/14Ts40tRi1cIcB8K04AGb79ZP6zpdX6dK/view?usp=sharing"
output = "ASW.phased.sample"
gdown.download(url, output, quiet=True, fuzzy=True)

```


When downloaded, move the `ASW.phased.haps` and `ASW.phased.sample` files to the **ADMIX_COHORT** folder in **tutorial_data**

**Using a local copy of the oputput**

In case we are pressed for time, one can copy the output from the location below. The `ASW.phased` file needs to be copied to **ADMIX_COHORT/** folder.

```
cp /gpfs/fs1/home/d/dfelsky/milicmil/phase/ASW.phased.haps /gpfs/fs0/scratch/SPONSOR_INITIAL/SPONSOR_NAME/YOUR_SCINET_USERNAME/tutorial_data/ADMIX_COHORT/

cp /gpfs/fs1/home/d/dfelsky/milicmil/phase/ASW.phased.sample /gpfs/fs0/scratch/SPONSOR_INITIAL/SPONSOR_NAME/YOUR_SCINET_USERNAME/tutorial_data/ADMIX_COHORT/
```


&nbsp;  
&nbsp;       

One can create also create a batch job out of this command. In order to do that, you will need to make an empty bash file which will store your command. To make the file, type the following code below:

```
nano phase_job.sh
```
Now that you are in nano text editor select the information in the code block below and paste it using *ctrl + insert*. you can use the arrow keys to navigate to the different lines. When the **USERNAME** has been modified, press **ctrl + x** in order to exit. You will get a prompt to confirm the save. Press **Enter** to confirm the save.

In the code below, you will need to specify the location of your output (all the print statements and whatnot that the batch job generates) as well as to change the directory in the script to your `tutorial_data` folder. **Note that you need to load the modules again in a batch job in Sci-net.**

```
#!/bin/bash
#SBATCH --nodes=2
#SBATCH --cpus-per-task=40
#SBATCH --time=1:30:00
#SBATCH --job-name=phsing_tractor_USERNAME
#SBATCH --output=/DIRECTORY PATH TO OUTPUT/omp_output_%j.txt
#SBATCH --mail-type=FAIL

module use /gpfs/fs1/home/d/dfelsky/dfelsky/modules/

module load shapeit

cd /gpfs/fs0/scratch/SPONSOR_INITIAL/SPONSOR_NAME/YOUR_SCINET_USERNAME/tutorial_data

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

shapeit  --input-vcf ADMIX_COHORT/ASW.unphased.vcf.gz \
      --input-map HAP_REF/chr22.genetic.map.txt \
      --input-ref HAP_REF/chr22.hap.gz HAP_REF/chr22.legend.gz HAP_REF/ALL.sample \
      -O ADMIX_COHORT/ASW.phased

```
This batch job should take around 70 minutes on Sci-net.

To start your batch job type in `sbatch phase_job.sh`.
To monitor the status of your job type in ```squeue -u USERNAME``` to see the status of your batch job.

Depending on the number of requests, there may be a few hour wait before your job is run.

#### Converting  shapeit output to vcf format.

Shapeit provides a convenient function to convert from its `haps`/`sample` file format to `vcf` format. Notice in the new vcf file we just created, slashes have changed to the pipe or vertical bar to seperate genotype.

**Converting shapeit output to vcf format**
```       
shapeit -convert \
        --input-haps ADMIX_COHORT/ASW.phased\
        --output-vcf ADMIX_COHORT/ASW.phased.vcf
```   
This command takes 7 seconds to complete

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
The above command takes around 5 minutes using debug node.
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
The above command takes around 5 minutes using debug node.


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
python Tractor/RunTractor.py --hapdose ADMIX_COHORT/ASW.phased --phe PHENO/Phe.txt --method linear --out SumStats.tsv
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

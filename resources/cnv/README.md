# Calculating CNVs using Illumina 450k data

In order to get the CNV status we can use the relative intensities of each probe against a common set of controls. It's a slow running algorithm and so it might be best to split up the task, and luckily each sample can be calculated independently without loss of accuracy. These scripts allow you to split the data in two ways.


## Setup

You need the `.idat` files all present in a directory. For example, look at this example directory:

    GSM1338100_6057825094_R01C01_Grn.idat
    GSM1338100_6057825094_R01C01_Red.idat
    GSM1338101_6057825094_R01C02_Grn.idat
    GSM1338101_6057825094_R01C02_Red.idat
    GSM1338102_6057825094_R02C01_Grn.idat
    GSM1338102_6057825094_R02C01_Red.idat
    GSM1338103_6057825094_R02C02_Grn.idat
    GSM1338103_6057825094_R02C02_Red.idat
    GSM1338104_6057825094_R03C01_Grn.idat
    GSM1338104_6057825094_R03C01_Red.idat

These 10 files represent 5 samples - a red and green channel for each sample. The script is going to try to come up with a unique sample ID from the file names, in this case it will assume that `GSM133810*` is the sample name, and everything after the underscore will be assumed to be batch information. 

You also need the following CRAN libraries:

- `plyr`
- `splitstackshape`

The following Bioconductor libraries:

- `CopyNumber450k`
- `CopyNumber450kData`

And the `meffil` package. e.g.:

    install.packages("plyr")
    install.packages("splitstackshape")
    install.packages("devtools")
    source("http://bioconductor.org/biocLite.R")
    biocLite("CopyNumber450k")
    biocLite("CopyNumber450kData")
    library(devtools)
    install_github("perishky/meffil")

That should be all that's needed.

## 1. Using multiple cores on a single node

Navigate to the `cnv` folder:

    cd godmc/resources/cnv

To calculate the CNVs using 16 cores, for `.idat` files located at `/path/to/idatfiles`, and to save them to the filename `cnv_filename.RData`, run the following command:

    R --no-save --args /path/to/idatfiles/ ../../raw_data/cnv/cnv_filename.RData 16 < estimate_cnv.R

It takes about 5 minutes per sample, so depending on your sample size this may take a while. Alternatively you can distribute using a job scheduler...


## 2. Using multiple nodes based on a job scheduler

An example of how to do this in parallel on our PBS system is provided in the `godmc/resources/cnv/estimate_cnv.sh` script. You will need to edit this script to be specific to your submission system and file paths etc. The instructions on how to do this are written as comments in the script itself.

Once you have edited the file you can test if it is going to work by running it interactively:
    
    cd godmc/resources/cnv
    ./estimate_cnv.sh 1

This will begin running the first batch. Replace `1` with another number to run that batch instead. If it seems to be working ok then submit the script. *e.g.* with PBS I would run:

    qsub estimate_cnv.sh

This will result in a number of files being generated in the `godmc/raw_data/cnv/` folder. For example, if we had an `output_name` called `cnv_filename` and there were 5 batches, then we would have the following files: 

    cnv_filename_1.RData
    cnv_filename_2.RData
    cnv_filename_3.RData
    cnv_filename_4.RData
    cnv_filename_5.RData

Check that all the files that you expect are indeed there.

The next step is to aggregate these files into a single data file. To do this run the `godmc/resources/cnv/aggregate_cnv.R` script:

    R --no-save --args ../../raw_data/cnv/cnv_filename < aggregate_cnv.R

Here we are specifying the prefix to the batch filenames.
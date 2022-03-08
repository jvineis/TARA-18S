# TARA-18S
### Most of the work outlined here is from my work with Bess Ward and others in her lab at Princeton. Its just a reference for us to recreate our analysis, but could be useful to others. Its  unlikely that I will maintain the scripts in this directory, but I might be able to help if you are really stuck. Keep in mind that 18S is probs really not very good for much and many organisms have tons of copies of this gene within a single cell. Anyway. Its fun to analyze and sometimes we learn things. 
### I use conda to download almost all software needed for the analysis below unless otherwise noted.
#### 1.  Download the data from SRA.All you need is an SRA number and sratools to get download the fastq files for analysis. The SRA numbers should be in a list like the one below. 

    ERR562490
    ERR562503 
    ERR562426
    ERR562657
    ERR562473
    ERR562622

##### then you can run a script like this to download the data.. this will create a file for both the read1 and read2 fastq files

    #!/bin/bash
    for i in `cat x_sra-18S-names-to-download.txt`; do fastq-dump $i --split-files; done

####  2. Create .ini files for each of the samples which will be used to run the merging script.

    ls *_1.fastq > 1
    ls *_2.fastq > 2
    ls *_1.fastq | cut -f 1 -d "_" > 3
    paste 3 1 2 > x_file-for-iu-gen-configs.txt

####  you need to add in a header for the x_file-for-iu-gen-configs.txt, so your file looks like this

    sample	r1	r2
    ERR562370	ERR562370_1.fastq	ERR562370_2.fastq
    ERR562382	ERR562382_1.fastq	ERR562382_2.fastq
    ERR562390	ERR562390_1.fastq	ERR562390_2.fastq
    ERR562426	ERR562426_1.fastq	ERR562426_2.fastq
    ERR562473	ERR562473_1.fastq	ERR562473_2.fastq
    
#### 3. Generate an ini file for each of the samples in your x_file-for-gen-configs.txt

    iu-gen-configs x_file-for-gen-configs.txt
    
####  4. Then you can merge your reads.  NOTE: these sequences are not directional (forward and reverse primers can be found in both read1 and read2.. so you can't trim the adaptor using the integrated adapter filtering flag in the merging script. This will need to be done separately. 

    #!/bin/bash
    #
    #SBATCH --nodes=1
    #SBATCH --tasks-per-node=20
    #SBATCH --mem=100Gb
    #SBATCH --time=05:00:00
    #SBATCH --array=1-25

    SAMPLE=$(sed -n "$SLURM_ARRAY_TASK_ID"p samples.txt)

    iu-merge-pairs ${SAMPLE}.ini --enforce-Q30-check -o ${SAMPLE} --ignore-deflines
    
####  5. This will merge your reads and create a ton of output files. Now you are going to want to filter out the primer sequences. Which you can do using the script "


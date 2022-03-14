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
    
####  5. This will merge your reads and create a ton of output files. Now you are going to want to filter out reads that don't contain the primer sequences. Which you can do using the script below

    #!/bin/bash#
    #SBATCH --nodes=1
    #SBATCH --tasks-per-node=1
    #SBATCH --mem=100Gb
    #SBATCH --time=00:20:00
    #SBATCH --array=1-25
    
    SAMPLE=$(sed -n "$SLURM_ARRAY_TASK_ID"p samples.txt)
    
    python ~/scripts/filter-for-primer.py --i ${SAMPLE}_MERGED --o ${SAMPLE}_MERGED-primer-filtered.fa

#### 6. Now that you have the quality and ensured that all of your amplicons contain the proper adaptor, Its time to run swarm.  This consists of a few parts; 1. concatenate all of the filtered fasta files, 2. use vsearch to dereplicate all sequences 3. run swarm 4. Create a table of counts of each SWARM in each sample. 

    #!/bin/bash
    #SBATCH --nodes=1
    #SBATCH --tasks-per-node=20 
    #SBATCH --mem=200Gb
    #SBATCH --time=24:00:00

    cat *MERGED-primer-filtered.fa > pooled-samples.fa
    vsearch --derep_fulllength pooled-samples.fa --sizeout --output pooled-samples-derep.fa
    swarm -d 1 -f -t 40 -z pooled-samples-derep.fa -s pooled-samples-derep-stats.txt -w pooled-samples-node-representatives.fa -o pooled-samples-node-table.txt
    python ~/scripts/mu-swarms-to-ASVs-table-for-tarra.py -s pooled-samples-node-table.txt -o swarm-min50-count.txt -l samples.txt -n pooled-samples-node-representatives.fa -min 50

#### 7. Run the taxonomy on the representative nodes using vsearch

    #!/bin/bash
    #
    #SBATCH --nodes=1
    #SBATCH --tasks-per-node=1
    #SBATCH --mem=100Gb
    #SBATCH --time=00:20:00

    vsearch --usearch_global reduced-node-fasta-min50.fa --db /scratch/gpfs/WARD/DBs/Database_W2_v9_pr2.fasta --blast6out NODE-HITS-min50.txt --id 0.6
    
#### 8. We need to transpose the swarm-min50-count.txt table and create a tree file based on the relative abundance of each swarm in the table - euclidian distances based on bray-curtis dissimilarity. The R script "x_rscript-to-build-tree-from-node-table.R" can do both of these things and is executed using the bash script below. You will need to edit the R script script so that your file names are contained in the "dat" and "write.tree" and "write.table" lines. 

    #!/bin/bash
    #SBATCH --nodes=1
    #SBATCH --tasks-per-node=1
    #SBATCH --mem=100Gb
    #SBATCH --time=05:00:00

    Rscript x_rscript-to-build-tree-from-node-table.R
    
#### 9. The output of R script file doesn't have the proper header. You need to open the file and add "ASVs" to the first row, first column. The first row should look something like this. The "ERR.." column headers are the sample names.

    ASVs	ERR562370	ERR562382	ERR562390	ERR562426   ....

#### 10. Merge the taxonomy, and count matrix to create a beautiful anvio table for data exploration using the script "convert-node-hits-to-tax-node-table.py"

    python convert-node-hits-to-tax-node-table.py -n NODE-HITS-min50.txt -o swarm-taxonomy-and-counts.txt -r W2_v9_pr2-tax.txt -a swarm-min50-count-for-anvio.txt
    
##### you should run this on a server node and that sbatch script should look something like this.

    #!/bin/bash
    #SBATCH --nodes=1
    #SBATCH --tasks-per-node=1
    #SBATCH --mem=100Gb
    #SBATCH --time=00:20:00
    python convert-node-hits-to-tax-node-table.py -n NODE-HITS-min50.txt -o swarm-taxonomy-and-counts.txt -r W2_v9_pr2-tax.txt -a swarm-min50-count-for-anvio.txt
    
#### 11. Now you should be able to load the files into anvio and visualize the abundance and taxonomy of your swarms (ASVs). Its helpful to run it from the server. In which case you will need to ssh in a special way. like thus. 

    ssh -L 8083:localhost:8083 jv2474@della.princeton.edu

##### then cd to your directory where you have been doing all of your good 18S work, activate anvio and then run the command to get the interactive display up and running. The files that you specify with the -d and -t flags come from steps 10 and 8 respectively. 

    anvi-interactive -d swarm-taxonomy-and-counts.txt -p swarm-taxonomy-and-counts.db -t swarm-min50-count.tre --manual-mode --server-only -P 8083

##### now in a web browser, type in the following and the display will magically appear! Enjoy!

    http://0.0.0.0:8083




    

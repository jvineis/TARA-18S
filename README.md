# TARA-18S
### Most of the work outlined here is from my work with Bess Ward and others in her lab at Princeton. Its just a reference for us to recreate our analysis, but could be useful to others. Its  unlikely that I will maintain the scripts in this directory, but I might be able to help if you are really stuck. Keep in mind that 18S is probs really not very good for much and many organisms have tons of copies of this gene within a single cell. Anyway. Its fun to analyze and sometimes we learn things. 
### I use conda to download almost all software needed for the analysis below unless otherwise noted.
### SRA numbers for the TARA oceans data can be found here for most samples https://www.ncbi.nlm.nih.gov/Traces/study/?query_key=8&WebEnv=MCID_623338d644660f5bac78969b&o=acc_s%3Aa and you can download metadata for samples of interest using sratools for any or all of the samples if you like. The metadata is included in this git for you as well. 

#### 1.  Download the data from SRA. All you need is an SRA number and sratools to get download the fastq files for analysis. The SRA numbers should be in a list like the one below which can be found in the metadata (link above) for all the samples that you want to analyze. 

    ERR562490
    ERR562503 
    ERR562426
    ERR562657
    ERR562473
    ERR562622

##### then you can run a script like this to download the data.. this will create a file for both the read1 and read2 fastq files. If you are working on the princeton cluster, you need to ssh jv2474@tigressdata.princeton.edu

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
    
    
####  5. Step 4 will merge your reads and create a ton of output files. Now you are going to want to filter out reads that don't contain the primer sequences from the high quality read files that end with MERGED. You can do using the script below. The script below also contains the step to dereplicate the reads for each of your samples. Make sure that you have an active bioconda environment (e.g. "conda activate bioconda") prior to running python script through the sbatch. Also make sure that you have vsearch active prior to running the derepliate step for each sample. You can use the hash (#) character in front of the lines that you don't want to run.

    #!/bin/bash#
    #SBATCH --nodes=1
    #SBATCH --tasks-per-node=1
    #SBATCH --mem=100Gb
    #SBATCH --time=00:20:00
    #SBATCH --array=1-25
    
    SAMPLE=$(sed -n "$SLURM_ARRAY_TASK_ID"p samples.txt)
    
    python ~/scripts/filter-for-primer.py --i ${SAMPLE}_MERGED --o ${SAMPLE}_MERGED-primer-filtered.fa
    vsearch --quiet --derep_fulllength ${SAMPLE}_MERGED-primer-filtered.fa --sizeout --fasta_width 0 --relabel_sha1 --output ${SAMPLE}-primer-derep.fa

##### Lets do a little checking to make sure that our output makes sense. Lets first look at the number of sequences that passed merging/quality and primer filtering. Here is a way to get these numbers 
##### 1. Save the output of these two commands below as 1 and 2
     
    for i in *primer-filtered.fa ; do echo $i | cut -f 1 -d "_"; done
    for i in *primer-filtered.fa ; do grep ">" $i | wc -l; done
    
##### 2. paste the two together 

    paste 1 2 > x_Quality-and-primer-filtered-reads-per-sample.txt

##### This file should look something like this - just a sanity check to make sure we didn't loose all of our reads. 

    ERR562370	868285
    ERR562382	2028652
    ERR562390	1027926
    ERR562426	801888
    ERR562473	699206


#### 6. Now that you have the quality and ensured that all of your amplicons contain the proper adaptor, Its time to run swarm.  This consists of a few parts that are outlined in the bash script below. 

    #!/bin/bash
    #
    #SBATCH --nodes=1
    #SBATCH --tasks-per-node=1
    #SBATCH --mem=150Gb
    #SBATCH --time=04:00:00

    ### These steps replicate the work here
    ### https://github.com/frederic-mahe/swarm/wiki/Fred's-metabarcoding-pipeline#global-dereplication-clustering-and-chimera-detection
    ### The cat and vsearch steps should be run with the conda vsearch envionment
    ### The sarm should be run with this environment "conda activate /home/jv2474/.conda/envs/swarm-v3.1"
    ### The python script should be run with the bioconda environment.

    ## 1. Concatenate the merged and filtered sequences
    #cat *-primer-derep.fa > pooled-samples.fa
    ## 2. Dereplicaete the concatenated sequences
    #vsearch --derep_fulllength pooled-samples.fa --fasta_width 0 --sizeout --sizein --output pooled-samples-derep.fa
    ## 3. Cluster the sequences
    #swarm -d 1 -f -t 40 -z pooled-samples-derep.fa -i pooled-samples-derep-struct.txt -s pooled-samples-derep-stats.txt -w pooled-samples-node-representatives.fa -o pooled-samples-node-table.txt
    ## 4. Sort the clustered node representatives
    #vsearch --fasta_width 0 --sortbysize pooled-samples-node-representatives.fa --output pooled-samples-node-representatives-sorted.fa
   
#### 7. Run the taxonomy on the representative nodes.

    #!/bin/bash
    #
    #SBATCH --nodes=1
    #SBATCH --tasks-per-node=1
    #SBATCH --mem=100Gb
    #SBATCH --time=00:20:00

    vsearch --usearch_global pooled-samples-node-representatives-sorted.fa --db /scratch/gpfs/WARD/DBs/Database_W2_v9_pr2.fasta --blast6out NODE-HITS.txt --id 0.6
    
#### 8. Convert the swarm output to a contingency table. Then create two tables from the resulting output. 1. contains the metadata for each swarm including taxonomy, length of the representative sequence etc. 2. a contingency table of samples and swarms. The script also create a tree file for both the sample organization and the amplicon organization. you could also reconstruct a phylogenetic tree if you wanted to  
    #!/bin/bash
    #
    #SBATCH --nodes=1
    #SBATCH --tasks-per-node=1
    #SBATCH --mem=150Gb
    #SBATCH --time=04:00:00 
    
    ## 5. Combine the swarms into an ASV table.
    #python ~/scripts/mu-swarms-to-ASVs-table-for-tarra.py -repfa pooled-samples-node-representatives-sorted.fa -stats pooled-samples-derep-stats.txt -swarms pooled-samples-node-table.txt -l samples-primer-derep-names.txt > x_SWARM-contingency-table.txt
    ## 6. Filter out the low abundance SWARMS, and create the file for anvio visualization.
    #python ~/scripts/convert-node-hits-to-tax-node-table.py -n NODE-HITS.txt -o x_SWARMS-and-tax-for-anvio -r W2_v9_pr2-tax.txt -s x_SWARM-contingency-table.txt -min 50

    
#### 8. We need to transpose the swarm-min50-count.txt table and create a tree file based on the relative abundance of each swarm in the table - euclidian distances based on bray-curtis dissimilarity. The R script "x_rscript-to-build-tree-from-node-table.R" can do both of these things and is executed using the bash script below. You will need to edit the R script script so that your file names are contained in the "dat" and "write.tree" and "write.table" lines. The script will produce a x_SWARMS-and-tax-for-anvio-relative-abundance-samples.tre file and a x_SWARMS-and-tax-for-anvio-relative-abundance.tre which are for your sample and ASV organization respectively. 

    #!/bin/bash
    #SBATCH --nodes=1
    #SBATCH --tasks-per-node=1
    #SBATCH --mem=100Gb
    #SBATCH --time=05:00:00

    Rscript x_rscript-to-build-tree-from-node-table.R

#### 9. Lets add some more detail to the display. The R script "x_rscript-to-build-tree-from-node-table.R" is run through "x_build-tree-from-node-table.shx" and will create a newick style tree that you can make readable by ANVIO and then display the samples and ASVs in a biologicaly meaningful order. The file needs to look something like below. One way to do this woud be to open the tree file and paste in all the text except for the newick tree section. easy peasy. If you are not sure of the name of your tree file.. just run ls \*.tre and you will find the file that you need to open and edit in this way. 

    item_name	data_type	data_value
    sample_order	basic	(((ERR562382:0.3305044361,ERR562490:0.3305044361):0.1622693878,((ERR562495:0.3165081272,ERR562721:0.31
    65081272):0.1529479025,((ERR562556:0.3701376285,ERR562622:0.3701376285):0.08540201303,(((ERR562483:0.2220207845,ERR562574:0.22
    20207845):0.08438854788,(ERR562620:0.2762364185,(ERR562539:0.1766348813,ERR562730:0.1766348813):0.09960153716):0.03017291388):
    0.1233936379,(ERR562722:0.3852795611,((ERR562576:0.1943682858,ERR562644:0.1943682858):0.1354090604,(ERR562370:0.2469266734,ERR
    562598:0.2469266734):0.08285067274):0.05550221498):0.04452340914):0.02573667126):0.01391638816):0.02331779422):0.005857568796,
    (ERR562390:0.4984154387,((ERR562517:0.3260482863,(ERR562616:0.2913913532,(ERR562553:0.247530561,ERR562672:0.247530561):0.04386
    07922):0.03465693311):0.1469173896,(ERR562657:0.3818427369,(ERR562503:0.2317135468,(ERR562426:0.2187443165,ERR562473:0.2187443
    165):0.01296923032):0.1501291901):0.09112293894):0.02544976278):0.0002159540589);

##### To add this information to an anvio database, you first need to create one, simply by loading the data that you have already made into an anvio interactive session. Here is how to do that.  

##### First start a fresh ssh using a login like the one you see bleow.
    
    ssh -L 8083:localhost:8083 jv2474@della.princeton.edu

##### then cd to your directory where you have been doing all of your good 18S work, activate anvio and then run the command to get the interactive display up and running. The files that you specify with the -d and -t flags come from steps 10 and 8 respectively. 

    anvi-interactive -d x_SWARMS-and-tax-for-anvio-relative-abundance.txt -t x_SWARMS-and-tax-for-anvio-relative-abundance.tre -p x_SWARMS-and-tax-for-anvio-relative-abundance.db --manual-mode --server-only -P 8083

##### now in a web browser, type in the following and the display will magically appear (after you click on Draw).

    http://0.0.0.0:8083

##### you can then stop the anvi-interactive display in your terminal using cntrl c

##### Now you can add the sample tree information to your display like this. (make sure that you have an active anvio conda environment.)

    anvi-import-misc-data x_SWARMS-and-tax-for-anvio-relative-abundance-samples.tre -p x_SWARMS-and-tax-for-anvio-relative-abundance.db --target-data-table layer_orders

##### Add any metadata to your anvio database. 

##### Lets say that I have a file x_SWARMS-and-tax-for-anvio-metadata.txt that was produced by the convert-node-hits-to-tax-node-table.py way up in step 7. I can now add this to the layers anvio database that will allow me to visualize the taxonomy and swarm node details of each swarm(ASV). After you import this information, you can reload your anvio display as you did above.

    anvi-import-misc-data x_SWARMS-and-tax-for-anvio-metadata.txt -p x_SWARMS-and-tax-for-anvio-relative-abundance.db --target-data-table items

##### The metadata file, x_SWARMS-and-tax-for-anvio-metadata.txt, should look something like this. 

    OTU	total	cloud	amplicon	length	abundance	spread	taxid	taxonomy	tax1	tax2	tax3	tax4	tax5	tax6	tax7	tax8	tax9	tax10 
    s_1	1367871	6663	b268468dd19dac760691ce28812afcd10270f7b7	167	442316	25	HO778380.2.1586_U	Eukaryota|Archaeplastida|Streptophyta|Embryophyceae|Embryophyceae_X|Embryophyceae_XX|Phaseolus|Phaseolus+acutifolius|na|na	Eukaryota Archaeplastida	Streptophyta	Embryophyceae	Embryophyceae_X	Embryophyceae_XX	Phaseolus	Phaseolus+acutifolius	na	na
    s_2	1060749	7677	df0716f210b687c98b2e052e82f5c73e3a1c146b	165	713680	25	HO778380.2.1586_U	Eukaryota|Archaeplastida|Streptophyta|Embryophyceae|Embryophyceae_X|Embryophyceae_XX|Phaseolus|Phaseolus+acutifolius|na|na	Eukaryota
	Archaeplastida	Streptophyta	Embryophyceae	Embryophyceae_X	Embryophyceae_XX	Phaseolus	Phaseolus+acutifolius	na	na



    

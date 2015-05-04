# 14_crispr
Code and processed data reported in Shi et al. (2015), Nat. Biotechnol. [in press]

#####################
RUNNING THE PIPELINE
#####################

To run the prevalence analysis on an example Cas9 data set, do the following.

1. Change to the "example" directory. This containes the following files:

- barcodes.txt: Lists the different barcodes used. Contains the following columns:
     column 1: barcode sequence
     column 2: barcode name

- experiments.txt: Lists which barcodes and which regions correspond to which experiments and which timepoints. Contains the following columns:
     column 1: experiment name
     column 2: region name
     column 3: barcodes for control sample, with ‘x’ denoting “any”.
     column 4,5,…: barcodes for timepoint samples. The corresponding header for the column must be “dXX” where XX is the day number, e.g. d3, d7, d10, etc.

- regions.txt: Lists the regions interrogated in the experiment.
     column 1: region name
     column 2: forward primer sequence
     column 3: reverse primer sequence
     column 4: sequence of the wild-type region beginning with the 5’ end of the forward primer and ending with the 5’ end of the reverse primer. Exonic nucleotides must be capitalized, with the left-most capitalized nucleotide indicating the start of the reading frame.

- r1.fastq and r2.fastq: the forward and reverse read fastq files. Currently the pipeline requires that all forward reads be in a single file, and similarly for all reverse reads. Starting with the *R1*.fastq.gz and *R2*.fastq.gz files provided by the illumina sequences, this can be accomplished by executing the following commands:

$ zcat *R1*.fastq.gz > r1.fastq
$ zcat *R2*.fastq.gz > r2.fastq

- infection_efficiency.txt: Lists the measured infection efficiency observed for each sample.
     column 1: experiment name
     column 2: infection efficiency observed for control sample
     column 3,4,.... infection efficiency observed for timepoint samples

2. Edit config.py appropriately. Specifically, set the correct location of the python you want to use in the first line, and change the following appropriately:

pipeline_dir 
python_to_use 
input_dir 

3. To make sure the analysis code can be run, change to the analysis directory and execute the following command ($ indicates the prompt):

$ chmod u+x config.py

4. To run the pipeline, edit config.py appropriately. Then do

$ ./config.py

This will create a subdirectory called "output" and fill it with the appropriate analysis stuff. Note: if you actually go through this analysis on the example data set, no plots will be created because there won't be enough reads to track. To actually get the plots that appear in the paper, run analysis on data in published_results/dataset_1/ and published_results/dataset_2/. 
 
#####################
PIPELINE DETAILS
#####################

config.py calls pipeline.py, which uses the rufus module and is divided into the following stages:

1. Split r1.fastq and r2.fastq into files containing 100K reads each. This allows the read processing to be efficiently farmed out to nodes

2. Create directories for pipeline output

3. Process reads from each pair of split fastq files into observed sequences

4. Combine observed sequences from each timepjoint into a single file and count the number of occurrences of each

5. For each experiment, combine the observed sequences obtained for every timepoint into a single file 

6. Tally sequence occurrences [WHY IS THIS NECESSARY?]

7. Summarize the sequences observed in each experiment. This is the step in which each sequence is classified as an in-frame (IFN), out-of-frame (OFN), etc. mutation. The algorithm that does this is in /data/kinney/jkinney/projects/14_crispr/pipeline/routine_summarize_seqs.py. 

8. Collect the experiment summaries into one place (typically output/summaries/). 

It is important to understand how sequences are classified in stage 7. The sequence classes annotated are:

     WTS: wild-type sequence     
     WTE: wild-type exon    
     IFM: in-frame sense mutation
     IFM: in-frame nonsense mutation
     OFM: out-of-frame mutation

The decision tree used to assign classes is as follows:

     Q: Is the observed exon length a multiple of 3?
     N -> OFM (out-of-frame mutation)
     Y -> Q: Is the observed exon identical to the wild-type exon?
               Y -> Q: Is the observed sequence identical to the wild-type sequence?
                         Y -> WTS (wild-type sequence)
                         N -> WTE (wild-type exon)
               N -> Q: Does the observed exon contain a stop codon?
                         Y -> IFN (in-frame nonsense mutation)
                         N -> IFM (in-frame sense mutation)

In particular, note that a sequence that has no exon is said to have exon length 0, and therefore all mutations observed within this sequence will be recorded as WTE. Also, the “frameshift” mutations described in the paper refer to both OFM and IFN mutations. 


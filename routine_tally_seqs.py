import scipy as sp
import glob
import os, sys, time

# Get arguments
experiment_dir = sys.argv[1]
experiments_file = sys.argv[2]

# Get list of timepoints
f = open(experiments_file)
atoms = f.readline().split()
timepoint_names = atoms[2:]
f.close()

# Load unique seqs into memory
experiment_unique_seqs_file = '%s/unique_seqs'%experiment_dir
seq_to_counts_dict = {}
f = open(experiment_unique_seqs_file)
for line in f.readlines():
    atoms = line.split()
    counts_list = [int(atoms[1])] + [0]*len(timepoint_names)
    seq = atoms[2]
    seq_to_counts_dict[seq] = counts_list
f.close()

# Determine coutns of each timepoint
for n, tp in enumerate(timepoint_names):
    sample_unique_seqs_file = '%s/%s/unique_seqs'%(experiment_dir,tp)
    f = open(sample_unique_seqs_file)
    for line in f.readlines():
        atoms = line.split()
        counts = int(atoms[1])
        seq = atoms[2]
        seq_to_counts_dict[seq][n+1] = counts
    f.close()

# Get list of sequences ordered by their timepoints
total_counts = [x[0] for x in seq_to_counts_dict.values()]
seqs = seq_to_counts_dict.keys()
indices = sp.argsort(total_counts)[::-1]

# Write results
tallied_seqs_file = '%s/tallied_seqs'%experiment_dir
f = open(tallied_seqs_file,'w')
line = 'total\t' + '\t'.join(timepoint_names) + '\tseq\n'
f.write(line)
for i in indices:
    seq = seqs[i]
    counts = seq_to_counts_dict[seq]
    line = '\t'.join([str(x) for x in counts]) + '\t' + seq + '\n'
    f.write(line)
f.close()



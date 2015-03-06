#!/usr/bin/python
import scipy as sp
import numpy as np
import string
import timeit
import os,sys


# Get input files
r1_file = sys.argv[1]
r2_file = sys.argv[2]
barcodes_file = sys.argv[3]
regions_file = sys.argv[4]
experiments_file = sys.argv[5]
output_dir = sys.argv[6]

extension = '.'.join(r1_file.split('.')[-2:])
extension2 = '.'.join(r2_file.split('.')[-2:])
assert(extension == extension2)

total_time_start = timeit.default_timer()

time_dict = {}
time_dict['align'] = 0
time_dict['trim'] = 0
time_dict['consensus'] = 0
time_dict['read_fastq'] = 0
time_dict['match'] = 0
time_dict['align_compute_scores'] = 0

# Function to compute the reverse complement of a sequence
complement = string.maketrans('ATCGN', 'TAGCN')
def rc(seq):
    return seq.upper().translate(complement)[::-1]

# Return the next read from a fastq file
time_get_next_read_from_fastq = 0.0
def get_next_read_from_fastq(f):
    start_time = timeit.default_timer()
    f.readline()
    read = f.readline().strip()
    f.readline()
    f.readline()
    time_dict['read_fastq'] += timeit.default_timer() - start_time
    return read

# Finds the barcode corresponding to each sequence
def match_barcode(seq, barcodes_dict,search_area=20):
    start_time = timeit.default_timer()
    tag_length = 0
    region = False
    for barcode in barcodes_dict.keys():
        k = seq.find(barcode,0,search_area)
        if k >= 0:
            region = barcodes_dict[barcode]
            tag_length = len(barcode)+k
            #continue
    time_dict['match'] += timeit.default_timer() - start_time            
    return (region, tag_length)

def findchar(s, ch):
    return [i for i, letter in enumerate(s) if letter == ch]

# Performs a gapless alginment of seq1 and seq2. Returns sequences padded with
# dashes as appropriate
def gapless_alignment(seq1, seq2):
    start_time = timeit.default_timer()
    
    L1 = len(seq1)
    L2 = len(seq2)
    dash_val = ord('-')
    
    # Convert sequences to an array of integers
    nseq1 = sp.array([ord(c) for c in seq1],dtype='int')
    nseq2 = sp.array([ord(c) for c in seq2],dtype='int')

    alignment_array_length = 2*L1+L2-2
    #positions_to_test = L1+L2-1

    # Position nseq2 in the middle of a2
    a2 = sp.ones(alignment_array_length,dtype='int')*ord('-')
    a2[L1:L1+L2] = nseq2

    # Find best alginment position

    # First, try finding alignment using heuristics
    a2_seq = ''.join([chr(c) for c in a2])
    k = -1
    i = 0
    while k < 0 and i < L1-end_length:
        k = a2_seq.find(seq1[i:end_length+i])-i
        i += 1
    
    # If heuristics found a match, use that alignment    
    if k >= 0:
        kbest = k
        
    # Otherwise, do costly alignment
    else:   
        #scores = [sum(nseq1 == a2[k:L1+k]) for k in range(positions_to_test)]
        #scores = compute_alignment_scores(nseq1,a2) 
        #kbest = sp.argmax(scores)
        kbest = 0

    # Position nseq1 in the optimal place of a1
    a1 = sp.ones(alignment_array_length,dtype='int')*ord('-')
    a1[kbest:kbest+L1] = nseq1

    # Trim excess '-' from ends
    indices = (a1 != dash_val) | (a2 != dash_val)
    a1_trimmed = a1[indices]
    a2_trimmed = a2[indices]

    # Convert back to string
    seq1_aligned = ''.join([chr(c) for c in a1_trimmed])
    seq2_aligned = ''.join([chr(c) for c in a2_trimmed])
    
    time_dict['align'] += timeit.default_timer() - start_time  

    # Return aligned sequences
    return seq1_aligned, seq2_aligned

## Compute alignment scores quickly
def compute_alignment_scores(nseq1,array2):
    start_time_2 = timeit.default_timer()
    L1 = len(nseq1)
    alignment_array_length = len(array2)
    positions_to_test = alignment_array_length - L1
    
    scores = sp.zeros(positions_to_test)
    for k in range(positions_to_test):
        scores[k] = sum(nseq1 == array2[k:L1+k])
    
    time_dict['align_compute_scores'] += timeit.default_timer() - start_time_2    
    return scores

# Gets consensus sequence of two aligned sequences, using the higher quality
# one for the overlap region
def get_consensus(aligned_seq1,aligned_seq2):
    start_time = timeit.default_timer()
    
    # Make sure sequences are the same length
    L = len(aligned_seq1)
    assert(L==len(aligned_seq2))
    
    # Convert sequences to an array of integers
    nseq1 = sp.array([ord(c) for c in aligned_seq1],dtype='int')
    nseq2 = sp.array([ord(c) for c in aligned_seq2],dtype='int')

    N_val = ord('N')
    dash_val = ord('-')
    overlap_indices = (nseq1 != dash_val)*(nseq2 != dash_val)
    
    num_Ns_1 = sum(nseq1[overlap_indices] == N_val)
    num_Ns_2 = sum(nseq2[overlap_indices] == N_val)    

    # Compute the three types of overlap region
    indices_only1 = (nseq2==dash_val) & (nseq1!=dash_val)
    indices_only2 = (nseq1==dash_val) & (nseq2!=dash_val)
    indices_overlap = (nseq1!=dash_val) & (nseq2!=dash_val)
    
    # Fill in values for consensus sequence
    nconsensus = sp.ones(L,dtype='int')*ord('-')
    nconsensus[indices_only1] = nseq1[indices_only1]
    nconsensus[indices_only2] = nseq2[indices_only2]
    nconsensus[indices_overlap] = nseq1[overlap_indices] if num_Ns_1 <= num_Ns_2 else nseq2[overlap_indices]
    
    # Determine whether to keep sequence
    num_overlap_matches = sum(nseq1[overlap_indices] == nseq2[overlap_indices])
    num_overlap_positions = sum(overlap_indices)
    frac_match = 1.0*num_overlap_matches/(num_overlap_positions + 1E-6)
    
    # Set consensus sequence
    if frac_match >= 0.5 or num_overlap_positions < 10:
        consensus = ''.join([chr(c) for c in nconsensus])
    else:
        consensus = ''
    
    time_dict['consensus'] += timeit.default_timer() - start_time  
    
    # Return consensus sequence
    return consensus
    
# Trimms 3' junk, assuming the second read has been reverse complemented
def trim_3p_junk(aligned_seq1, aligned_seq2):
    start_time = timeit.default_timer()
     
    # Make sure sequences are the same length
    L = len(aligned_seq1)
    assert(L==len(aligned_seq2)) 
    
    dash_val = ord('-')
    
    # Convert sequences to an array of integers
    nseq1 = sp.array([ord(c) for c in aligned_seq1],dtype='int')
    nseq2 = sp.array([ord(c) for c in aligned_seq2],dtype='int')
    
    junk1 = min(sp.where(nseq1!=dash_val)[0])
    junk2 = max(sp.where(nseq2!=dash_val)[0])
    
    seq1_trimmed = aligned_seq1[junk1:junk2+1]
    seq2_trimmed = aligned_seq2[junk1:junk2+1]
    
    time_dict['trim'] += timeit.default_timer() - start_time  
    return seq1_trimmed, seq2_trimmed

# Set other analysis parameters
end_length = 15

# Load barcodes into dictionary
f = open(barcodes_file)
barcodes_dict = {}
lines = f.readlines()
for line in lines[1:]:
    atoms = line.strip().split()
    if len(atoms) != 2:
        continue
    sample = atoms[0]
    barcode = atoms[1]
    barcodes_dict[barcode] = sample
f.close()

# Load regions into dictionary
f = open(regions_file)
f.readline() # Skip headder line
region_to_seq_dict = {}
read_end_to_region_dict = {}
for line in f.readlines():
    atoms = line.strip().split()
    if len(atoms) != 4:
        continue
    name = atoms[0]
    seq = atoms[3].upper()
    # User primers to get sequences of ends
    end_5p = atoms[1][:end_length]
    end_3p = atoms[2][:end_length]
    
    region_to_seq_dict[name] = seq
    read_end_to_region_dict[end_5p] = name+'_5p'
    read_end_to_region_dict[end_3p] = name+'_3p'

# Load experiment information
f = open(experiments_file,'r')
line = f.readline()
atoms = line.strip().split()
timepoints = atoms[2:]
output_files_by_key = {}

observed_seq_files = {}

for line in f.readlines():

    if len(line.strip())==0:
        continue;
    
    atoms = line.strip().split()
    experiment = atoms[0]
    region = atoms[1]
    barcode_names = atoms[2:]
    
    for timepoint_num, timepoint_name in enumerate(timepoints):
        barcode_name = barcode_names[timepoint_num]
        key = '%s_%s'%(region,barcode_name)
        output_file_name = '%s/%s/%s/observed_seqs.%s'%(output_dir,experiment,timepoint_name,extension)
        output_files_by_key[key] = open(output_file_name,'w')
        
        #print 'Creating %s...'%output_file_name

valid_keys = output_files_by_key.keys()

#
# Process sequences one-by-one: This is where the main processing happens
#

r1_f = open(r1_file)
r2_f = open(r2_file)
stop = False
num_successes = 0
num_reads = 0
num_fails_sample = 0
num_fails_region = 0
num_fails_alignment = 0
while not stop:
    # Increment number of reads
    if num_reads%1000 == 0 and num_reads > 0:
        total_time = timeit.default_timer() - total_time_start  
        successes_pct = 100*num_successes/num_reads
        fails_sample_pct = 100*num_fails_sample/num_reads
        fails_region_pct = 100*num_fails_region/num_reads
        fails_alignment_pct =  100*num_fails_alignment/num_reads
        print 'Total time: %d sec, Total reads: %d, Successes: %d%%, Sample fails: %d%%, Region fails: %d%%, Alignment fails: %d%%'%(total_time, num_reads, successes_pct, fails_sample_pct, fails_region_pct, fails_alignment_pct)
        #print time_dict
    
    num_reads += 1
    
    # Get reads
    read1 = get_next_read_from_fastq(r1_f)
    read2 = get_next_read_from_fastq(r2_f)
    
    if len(read1) == 0 or len(read2) == 0:
        stop = True
        continue
    
    # Determine sequence sample by matching barcode, and clip barcode
    sample1, tag_length1 = match_barcode(read1, barcodes_dict, search_area=15)
    sample2, tag_length2 = match_barcode(read2, barcodes_dict, search_area=15)    
    if not (sample1 and sample2):
        num_fails_sample += 1
        continue
    sample = sample1 if sample1 else sample2  
    read1_clipped = read1[tag_length1:]
    read2_clipped = read2[tag_length2:]

    ###
    ### I need a more robust way of determining region. 
    ###

    # Determine region by matching front of read1_clipped
    region1, tag_length1 = match_barcode(read1_clipped, read_end_to_region_dict, search_area=15)
    region2, tag_length2 = match_barcode(read2_clipped, read_end_to_region_dict, search_area=15)  
    #if not region:
    if (not region1) or (not region2) or (region1[-3:] == region2[-3:]):
        num_fails_region += 1
        continue
    region = region1 if region1 else region2 
    
    # Clip excess nucleotides from each end
    read1_clipped = read1_clipped[tag_length1-end_length:]
    read2_clipped = read2_clipped[tag_length2-end_length:]
    
    # Test gapless_alignment
    aligned1, aligned2 = gapless_alignment(read1_clipped, rc(read2_clipped))

    # Test trim_junk
    trimmed1, trimmed2 = trim_3p_junk(aligned1, aligned2)

    # Gest get_consensus
    consensus = get_consensus(trimmed1, trimmed2)
    if len(consensus) == 0:
        num_fails_alignment += 1
        continue
    
    # RC consensus if read1 matches to 3' end of region
    if '_3p' in region1:
        consensus = rc(consensus)    

    # Verify that first end_length and last end_length positions match region seq
    region_name = region[:-3]
    region_seq = region_to_seq_dict[region_name]
    if (not region_seq[:end_length]==consensus[:end_length]) or (not region_seq[-end_length:]==consensus[-end_length:]):
        num_fails_alignment += 1
        continue
   
    # Store reads associated with region
    key = region_name + '_x.' + sample.split('.')[1]
    if key in valid_keys:
        output_files_by_key[key].write(consensus+'\n')
        num_successes += 1
    else:
        num_fails_region += 1
        
for f in output_files_by_key.values():
    f.close()

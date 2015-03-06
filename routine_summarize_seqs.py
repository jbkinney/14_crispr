import scipy as sp
import os, sys
import glob

#print 'In summarize_seqs.py'

regions_file = sys.argv[1]
experiments_file = sys.argv[2]
experiment_dir = sys.argv[3]

reverse_dict = {'A':'T','C':'G','G':'C','T':'A','N':'N','a':'t','c':'g','g':'c','t':'a','n':'n'}

rc = lambda x: ''.join([reverse_dict[b] for b in x[::-1]])

mis_dist = 2 # Used to filter out SNPs

# Create region to sequence dictionary
f = open(regions_file)
region_to_seq_dict = {}
for line in f.readlines():
    atoms = line.strip().split()
    if len(atoms) == 0:
        continue;
    name = atoms[0]
    seq = atoms[3]
    region_to_seq_dict[name] = seq

# Create experiment ot region dictionary
f = open(experiments_file)
experiment_to_region_dict = {}
for line in f.readlines():
    atoms = line.strip().split()
    if len(atoms) == 0:
        continue
    experiment = atoms[0]
    region = atoms[1]
    experiment_to_region_dict[experiment] = region
    
# Get wt sequence for this sample
atoms = experiment_dir.split('/')
experiment = atoms[-1]
region = experiment_to_region_dict[experiment]
wt_seq_annotated = region_to_seq_dict[region]
wt_seq = wt_seq_annotated.upper()
wt_L = len(wt_seq)

# Get wt seq for comparison to each observe seq
wt_nseq = sp.array([ord(c) for c in wt_seq],dtype='int')

#print 'region: ' + region
#print 'experiment: ' + experiment
#print 'wt_seq_annotated: ' + wt_seq_annotated

# Define input and output files
in_file = experiment_dir+'/tallied_seqs'
out_file = experiment_dir+'/summarized_seqs'
#print 'Processing %s -> %s....'%(in_file,out_file)

# Open input file
f_in = open(in_file)

# Get header for tallied_seqs file
line = f_in.readline()
atoms = line.split()
count_names = '\t'.join(atoms[:-1])

# Open output file
f_out = open(out_file,'w')

# Write two headerlines, the first containing the wt sequence, the second listing field names
f_out.write('%s\t%s\n'%(region,wt_seq))
f_out.write('%s\tclass\tdL\tL_del\tL_ins\tsnps\tfbp\trbp\tindel_seq\n'%count_names)

# Compute exon start and exon stop positions
wt_exon_stop = len(wt_seq_annotated.rstrip('acgtn'))
wt_exon_start = len(wt_seq_annotated) - len(wt_seq_annotated.lstrip('acgtn'))
wt_exon = wt_seq[wt_exon_start:wt_exon_stop]
wt_exon_L = len(wt_exon)

codons = [wt_exon[i:i+3].upper() for i in range(0,wt_exon_L,3)]
stop_codons = [c in ['TAG','TAA','TGA'] for c in codons]
if any(stop_codons):
    print 'WT exon contains a stop codon! Abort!'
    for n, c in enumerate(codons):
        if c in ['TAG','TAA','TGA']:
            codons[n] = c.lower()
    print 'wt_exon = %s'%('.'.join(codons))
    #raise

#print 'preceeding exon: %s'%wt_seq_annotated[:wt_exon_start]
#print 'exon sequence: %s'%wt_seq_annotated[wt_exon_start:wt_exon_stop]
#print 'following exon: %s'%wt_seq_annotated[wt_exon_stop:]

# Process each sequence in unique_seq_file
for line in f_in.readlines():
    atoms = line.strip().split()
    if len(atoms) < 3:
        continue
    
    # Get observed sequence
    obs_seq = atoms[-1]
    obs_L = len(obs_seq)
    
    # Get multiplicities
    mult = '\t'.join(atoms[:-1])
    
    # Get indel length
    indel_length = len(obs_seq) - len(wt_seq)
    
    # Make arrays  and set variables for computing mismatch positions
    obs_nseq = sp.array([ord(c) for c in obs_seq],dtype='int')
    N_val = ord('N')
    L = min(len(wt_seq), len(obs_seq))
    
    # Compute compute fwd and rev mismatch positions
    fwd_mismatches = (obs_nseq[:L]!=wt_nseq[:L]) & (obs_nseq[:L]!=N_val)
    rev_mismatches = ((obs_nseq[-L:]!=wt_nseq[-L:]) & (obs_nseq[-L:]!=N_val))[::-1]
    breakpoint_found = False
    
    # Identify front break point fbp
    tmp = sp.where(fwd_mismatches)[0]
    i=0
    fbp = -1
    while (fbp < 0) and (i < len(tmp)-1):
        if tmp[i+1] <= tmp[i]+mis_dist: # To avoid SNPs, demand second mismatch within mis_dist of first
            fbp = tmp[i]
            breakpoint_found = True
        i += 1
    
    # Identify rear break point
    tmp = sp.where(rev_mismatches)[0]
    i=0
    rbp = -1
    while (rbp < 0) and (i < len(tmp)-1):
        if tmp[i+1] <= tmp[i]+mis_dist: # To avoid SNPs, demand second mismatch within mis_dist of first
            rbp = tmp[i]
        i += 1
    
    # If breakpoint is found
    if breakpoint_found:
        
        # if fbp occurs after rbp on wt_seq, set rbp to fbp position 
        if fbp > wt_L-rbp:
            rbp = wt_L-fbp
        
        # Compute indel sequence
        indel_seq = '%s'%obs_seq[fbp:-rbp]

        # Compute number of SNPs
        snps = sum(fwd_mismatches[:fbp]) + sum(rev_mismatches[:rbp])
        
    # Otherwise, if no breakpoint is found
    else:
        if not (fbp==-1) and (rbp==-1):
            print fbp, rbp
            raise
        
        # Set empty indel_seq
        indel_seq = ''
        
        # Compute number of SNPs
        snps = sum(fwd_mismatches)
    
    # Step 1: Is either exon boundary disrupted?
    if breakpoint_found and \
        (fbp <= wt_exon_start < wt_L - rbp or \
        fbp <= wt_exon_stop < wt_L - rbp):
            mut_type = 'EBD'
            
    # Otherwise, examine the observed exon
    else:
        
        # If there is no breakpoint...
        if not breakpoint_found:
            start = wt_exon_start
            stop = wt_exon_stop
            
        # Otherwise, if breakpoint ends before exon starts
        elif wt_L-rbp <= wt_exon_start:
            start = wt_exon_start + (obs_L - wt_L)
            stop = wt_exon_stop + (obs_L - wt_L)

        # Otherwise, if breakpoint starts after exon ends
        elif fbp >= wt_exon_stop:
            start = wt_exon_start
            stop = wt_exon_stop
            
        # Otherwise, is breakpoints are confined to exon
        elif wt_exon_start <= fbp and wt_L-rbp < wt_exon_stop:
            start = wt_exon_start
            stop = wt_exon_stop + (obs_L - wt_L)
            
        # Get sequence of observed exon
        obs_exon = obs_seq[start:stop]
        obs_exon_L = len(obs_exon)
    
        # If the observed exon is not of coding length
        if len(obs_exon)%3 != 0:    
            mut_type = 'OFM' # Out-of-frame mutation
            
        # Otherwise, if the observed exon is of a coding length,
        else:
            
            # Is the observed exon the same as wt?
            if obs_exon == wt_exon:
                
                # Is the observed sequence the same as wt?
                if obs_seq == wt_seq:
                    mut_type = 'WTS' # wild type sequence
                    
                #  If not, call wild type exon
                else: 
                    mut_type = 'WTE' # wild type exon
                    
            # Otherwise, does the observed exon contain a stop codon?
            else:
                codons = [obs_exon[i:i+3].upper() for i in range(0,obs_exon_L,3)]
                stop_codons = [c in ['TAG','TAA','TGA'] for c in codons]
                if any(stop_codons):
                    mut_type = 'IFN' # in frame nonsense mutation
                else:
                    mut_type = 'IFM' # in frame sense mutation

    # Compute deletion length                
    del_length = -(indel_length - len(indel_seq))
    
    # Compute insertion length
    ins_length = len(indel_seq)
    
    # Compute which reverse break point to display (-1 if none found)
    rbp_disp = wt_L-rbp if rbp >= 0 else -1
    
    # Write result to file
    output_line = '%s \t%s \t%d \t%d \t%d \t%d \t%d \t%d \t[%s]\n'%\
        (mult, mut_type, indel_length, del_length, ins_length, snps, fbp, rbp_disp, indel_seq)
    f_out.write(output_line)
    
# After processing each line, close file
f_out.close()
    
# After processing each file, signal completion 
#print 'Done!'





   

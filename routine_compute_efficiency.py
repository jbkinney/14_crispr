import scipy as sp
import matplotlib.pyplot as plt
from matplotlib.mlab import find
import sys
import glob
plt.close('all')

# Get commandline input. Requires 4 arguments
assert len(sys.argv) == 4
summary_dir = sys.argv[1]
infection_file = sys.argv[2]
efficiency_file = sys.argv[3]

# Number of reads required for each kept sequence
track_min_num_reads = 100
detect_min_num_reads = 10

# Read in set of experiments to analyze
f = open(infection_file)
infection_dict = {}
lines = f.readlines()

# Read in days
days = lines[0].strip().split()[1:]
num_days = len(days)

# Read in rates
for line in lines[1:]:
    atoms = line.strip().split()
    if len(atoms)==4:
        for n in range(3):
            key = atoms[0]+'_'+days[n]
            infection_dict[key] = float(atoms[n+1])/100.

# File to write output to
g = open(efficiency_file,'w')
string='experiment'.ljust(10) + '\tday\tI\tM\tm\tE\t3n+0\t3n+1\t3n+2\tIFM\tOFM\tshort|IFM'
g.write(string+'\n')

# Perform analysis for each summary file
summary_files = glob.glob(summary_dir+'/*.txt')
for summary_file in summary_files:
    
    # Get experiment name
    experiment = summary_file.split('/')[-1].split('_summary.txt')[0]

    # Read summary file into memory
    f = open(summary_file)
    summary_file_lines = f.readlines()

    # Get region from header line in file
    region = summary_file_lines[0].strip().split()[0]

    # Get days in summary file
    columns = summary_file_lines[1].strip().split()
    class_col_index = columns.index('class')
    dLs_col_index = columns.index('dL')
    these_days = columns[2:class_col_index]
    assert cmp(days,these_days) == 0
    assert class_col_index == num_days + 2
    assert dLs_col_index == num_days + 3
    
    # Read in counts and mutation types from file
    mut_types = []
    counts_list = []
    dLs = []
    for line in summary_file_lines[2:]:
        atoms = line.split()
        counts_list.append([float(x) for x in atoms[1:class_col_index]])
        mut_types.append(atoms[class_col_index])
        dLs.append(int(atoms[dLs_col_index]))
    num_seqs = len(mut_types)
    
    # Convert counts to numerical matrix
    counts = sp.array(counts_list)
    
    # Specify columns
    ctl_col = 0
    ctl_counts = counts[:,ctl_col]
    
    # Process sample for each day
    for day_num in range(num_days):

        # Get fration of cells infected
        if 'Rosa' in experiment:
            key = 'Rosa_' + days[day_num]
        else:
            key = experiment+'_'+days[day_num]
        
        # Check to see if have infection data. If not, just continue
        if not key in infection_dict.keys():
            continue;

        # Get infection rate
        I = infection_dict[key]
    
        t_col = day_num+1
        t_counts = counts[:,t_col]    
            
        # Compute epsilon, the mutation rate due to PCR
        epsilon = 1.0 - ctl_counts[0]/sp.sum(ctl_counts) 
    
        # Compute M, the mutation rate observed for t0
        M = 1.0 - t_counts[0]/sp.sum(t_counts) 
    
        # Note: element 0, the wt element, is skipped. 
        N_IFM = sp.sum([t_counts[i] for i in range(num_seqs) if mut_types[i]=='IFM']) 
        N_OFM = sp.sum([t_counts[i] for i in range(num_seqs) if mut_types[i]=='OFM']) 
        N_short_IFM = sp.sum([t_counts[i] for i in range(num_seqs) if (mut_types[i]=='IFM' and abs(dLs[i])<=9)])

	N_3np0 = sp.sum([t_counts[i] for i in range(1,num_seqs) if dLs[i]%3==0])
	N_3np1 = sp.sum([t_counts[i] for i in range(1,num_seqs) if dLs[i]%3==1])
	N_3np2 = sp.sum([t_counts[i] for i in range(1,num_seqs) if dLs[i]%3==2])

        # Compute crispr_efficiency
        E = (M-epsilon)/(I*(1.-epsilon))
    
        N = N_3np0 + N_3np1 + N_3np2
        f_3np0 = N_3np0/N
        f_3np1 = N_3np1/N
        f_3np2 = N_3np2/N
    
        N = N_IFM + N_OFM
        eps = 1E-10
        f_IFM = N_IFM/(N+eps)
        f_OFM = N_OFM/(N+eps)
        f_short_given_IFM = N_short_IFM/(N_IFM+eps)
    
        # Compute the fraction of each type of mutation
        experiment_padded = experiment.ljust(10)
        string =  '%s\t%s\t%d%%\t%d%%\t%d%%\t%d%%\t%d%%\t%d%%\t%d%%\t%d%%\t%d%%\t%d%%'%\
            (experiment_padded, days[day_num], 100.*I, 100.*M, 100.*epsilon,100*E, 100*f_3np0, 100*f_3np1, 100*f_3np2, \
            100*f_IFM, 100*f_OFM, 100*f_short_given_IFM)
        g.write(string+'\n')
g.close()
    

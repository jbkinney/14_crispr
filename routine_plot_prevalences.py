#!/opt/hpc/pkg/python-2.7/bin/python -W ignore
import scipy as sp
import matplotlib.pyplot as plt
from matplotlib.mlab import find
import glob
import sys
#from scipy.stats import mannwhitneyu
plt.close('all')

# arguments: 
# 1. summary file
# 2. output directory
# 3. display type

# EXAMPLE INPUT
# > python plot_prevalences.py summaries/Smarca4_e2.1_summary.txt plots/ area

# Number of reads required for each kept sequence
track_min_num_reads = 100
detect_min_num_reads = 10

# Plotting parameters
red = [.9,.2,.2]
green = [.2,.7,.2]
blue = [0,0,1]
black = [0,0,0]
alpha = 0.15
fontsize=12
fig_size_inches = [3.5, 3.5]
#plt.rc('font', family='serif', size=fontsize)

# Get summary file from input
assert len(sys.argv) == 4
summary_file_name = sys.argv[1]
output_dir = sys.argv[2]
disp_type = sys.argv[3] 

assert (disp_type in ['area','lines']), "Invalid display type"

#summary_file_name = 'summaries/Smarca4_e2.1_summary.txt'

#disp_type = 'lines'
#disp_type = 'area'

# Open file for reading
f = open(summary_file_name)
<<<<<<< Updated upstream
experiment = '_'.join(summary_file_name.split('/')[-1].split('_')[:-1])
=======
experiment = summary_file_name.split('/')[-1].split('_summary.txt')[0]
>>>>>>> Stashed changes

# Read in wt sequence and region name
line = f.readline()
atoms = line.strip().split()
region = atoms[0]
wt_seq = atoms[1]

# Check to see if region is noncoding
is_nc = '_nc' in region

# Read in timepoints
line = f.readline()
atoms = line.strip().split()
class_col = atoms.index('class')
timepoint_names = atoms[2:class_col]
num_timepoints = len(timepoint_names)

# Read in timepoint values
timepoints = sp.zeros(num_timepoints)
for i, tn in enumerate(timepoint_names):
    assert tn[0]=='d'
    timepoints[i] = float(tn[1:])
    
# Now read in all the counts
mut_types = []
counts_list = []
for line in f.readlines():
    atoms = line.split()
    counts_list.append([float(x) for x in atoms[1:class_col]])
    mut_types.append(atoms[class_col])
num_seqs = len(mut_types)
    
# Convert counts to numerical matrix
counts = sp.array(counts_list) + 0.5
    
ctl_col = 0
t0_col = 1
sample_cols = range(1,num_timepoints+1)

ctl_counts = counts[:,ctl_col]
t0_counts = counts[:,t0_col]
sample_counts = counts[:,sample_cols] 

# Get index of wt sequence
seq_is = sp.arange(num_seqs)
wt_is = [mut_types[i] in ['WTS'] for i in seq_is]
    
# Compute enrichment wrt t0
enrichment = sample_counts/t0_counts[:,sp.newaxis]
wt_index = find(wt_is)[0]
    
wt_enrichment = enrichment[wt_index,:] 
norm_enrichment = enrichment/wt_enrichment[sp.newaxis,:]
    
ctl_count_ceiling = 5
t0_count_floor = 100
    
ctl_ok_is = ctl_counts < ctl_count_ceiling
t0_ok_is = t0_counts >= t0_count_floor

# if region is noncoding, only set nc_is
nc_is = sp.array([mut_types[i] in ['WTE'] for i in seq_is])
mis_is = sp.array([mut_types[i] in ['IFM'] for i in seq_is])
non_is = sp.array([mut_types[i] in ['OFM','IFN'] for i in seq_is])

black_is = nc_is * t0_ok_is
red_is = non_is * ctl_ok_is * t0_ok_is
green_is = mis_is * ctl_ok_is * t0_ok_is
    
# Create figure
plt.figure(figsize=fig_size_inches)
ax = plt.subplot(1,1,1)
    
# black
if is_nc:
    black_enrichment = norm_enrichment[black_is,:]
    num_black = black_enrichment.shape[0]
    low_black_enrichment = sp.percentile(black_enrichment,q=25,axis=0)
    med_black_enrichment = sp.percentile(black_enrichment,q=50,axis=0)
    high_black_enrichment = sp.percentile(black_enrichment,q=75,axis=0)

else:        
    # green
    green_enrichment = norm_enrichment[green_is,:]
    num_green = green_enrichment.shape[0]
    assert num_green > 0, "There are no inframe mutations to track"
    low_green_enrichment = sp.percentile(green_enrichment,q=25,axis=0)
    med_green_enrichment = sp.percentile(green_enrichment,q=50,axis=0)
    high_green_enrichment = sp.percentile(green_enrichment,q=75,axis=0)

    # red
    red_enrichment = norm_enrichment[red_is,:]
    num_red = red_enrichment.shape[0]
    assert num_red > 0, "There are no frameshift mutations to track"
    low_red_enrichment = sp.percentile(red_enrichment,q=25,axis=0)
    med_red_enrichment = sp.percentile(red_enrichment,q=50,axis=0)
    high_red_enrichment = sp.percentile(red_enrichment,q=75,axis=0)
    
# If want lines, plot lines
if disp_type == 'lines':
    if is_nc:
        plt.semilogy(timepoints, black_enrichment.T, '-', color=black, linewidth=0.5, alpha=alpha)        
    else:
        plt.semilogy(timepoints, green_enrichment.T, '-', color=green, linewidth=0.5, alpha=alpha)
        plt.semilogy(timepoints, red_enrichment.T, '-', color=red, linewidth=0.5, alpha=alpha)

# If want area, plot area
elif disp_type == 'area':
    if is_nc:
        ax.fill_between(timepoints, low_black_enrichment, high_black_enrichment, facecolor=black, edgecolor='none',alpha=0.3)
    else:
        ax.fill_between(timepoints, low_green_enrichment, high_green_enrichment, facecolor=green, edgecolor='none',alpha=0.3)
        ax.fill_between(timepoints, low_red_enrichment, high_red_enrichment, facecolor=red, edgecolor='none',alpha=0.3)
   
t_min = min(timepoints)
t_max = max(timepoints)                          
x_placement = t_min 
      
# Plot solid lines
if is_nc:
    plt.semilogy(timepoints, med_black_enrichment.T, '-o', color=black, linewidth=1.5, markeredgecolor='none', markerfacecolor=black, markersize=4)
    plt.text(x_placement,0.3,'%d'%num_black, color=black, fontsize=fontsize)

else:
    plt.semilogy(timepoints, med_green_enrichment.T, '-o', color=green, linewidth=1.5, markeredgecolor='none', markerfacecolor=green, markersize=4)
    plt.semilogy(timepoints, med_red_enrichment.T, '-o', color=red, linewidth=1.5, markeredgecolor='none', markerfacecolor=red, markersize=4)
    
    # Do test               
    # for i in range(1,num_timepoints):
    #    u,p = mannwhitneyu(red_enrichment[:,i], green_enrichment[:,i])    
    #    yplacement = 1.5
    #    if p < 0.001:
    #        plt.text(timepoints[i],yplacement,'***',horizontalalignment='center',fontsize=fontsize)
    #    elif p < 0.01: 
    #        plt.text(timepoints[i],yplacement,'**',horizontalalignment='center',fontsize=fontsize)
    #    elif p < 0.05:
    #        plt.text(timepoints[i],yplacement,'*',horizontalalignment='center',fontsize=fontsize)
    
    plt.text(x_placement,0.3,'%d'%num_green, color=green, fontsize=fontsize)
    plt.text(x_placement,0.1,'%d'%num_red, color=red, fontsize=fontsize) 
        
pct_editing = (1. - counts[0,1]/sum(counts[:,1]))*100
pct_editing_control = (1. - counts[0,0]/sum(counts[:,0]))*100

plt.title('%s'%experiment)
plt.ylim([2E-3,1E1])
plt.xlim([t_min-0.5,t_max+0.5])
plt.xticks(timepoints,timepoint_names)
        
ax.xaxis.set_tick_params(size=2)
ax.yaxis.set_tick_params(size=2)

#plt.tight_layout()
#plt.show()

plt.savefig('%s/prevalence_%s_%s.pdf'%(output_dir,disp_type,experiment))

import os, sys, time, commands
import glob
from ruffus import *

# The first thing the pipeline does is run the specified configuration file
#assert len(sys.argv) == 2
#config_file = sys.argv[1]
#execfile(config_file)

# Set names for various files
read1_split_file_glob = reads_dir + '/r1.fastq.*'
read2_split_file_glob = reads_dir + '/r2.fastq.*'
observed_seq_files_split = results_dir+'/*/*/observed_seqs.*'
observed_seq_files_by_sample = results_dir+'/*/*/observed_seqs'
unique_seq_files_by_sample = results_dir+'/*/*/unique_seqs'
observed_seq_files_by_experiment = results_dir+'/*/observed_seqs'
unique_seq_files_by_experiment = results_dir+'/*/unique_seqs'
tallied_seq_files = results_dir+'/*/tallied_seqs'
summarized_seq_files = results_dir+'/*/summarized_seqs'

summarize_seqs_script = python_to_use  + ' ' + pipeline_dir + '/routine_summarize_seqs.py'
parse_seqs_script = python_to_use  + ' ' + pipeline_dir + '/routine_parse_seqs.py'
tally_seqs_script = python_to_use  + ' ' + pipeline_dir + '/routine_tally_seqs.py'
plot_prevalences_script = python_to_use  + ' ' + pipeline_dir + '/routine_plot_prevalences.py'
compute_efficiencies_script = python_to_use  + ' ' + pipeline_dir + '/routine_compute_efficiency.py'

# main_pipeline.py must be run from a "run_pipeline.py" script. run_pipeline.py defines a small number of global variables that are essential for the operation of main_pipeline.py. 

# Useful way to give feedback
def give_feedback(feedback):
    func_name =sys._getframe(1).f_code.co_name
    print '\nIn '+func_name+': '+feedback,
    sys.stdout.flush()

# Submits a list of scripts to be run as separate jobs
# Waits for all jobs to complete before continuing
def submit_and_complete_jobs(scripts, use_multiple_nodes):
    
    # Clear scripts dir
    os.system('rm -r %s'%scripts_dir)
    os.system('mkdir %s'%scripts_dir)

    if use_multiple_nodes:

        # Write and execute scripts
        for n, script in enumerate(scripts):
            script_name = scripts_dir + '/script_num_%d.sh'%n 
            f = open(script_name, 'w')
            f.write(script)
            f.close()
            os.system('chmod +x '+script_name)
            os.system('qsub -cwd -e %s -o %s %s > .junk' % (script_name+'.e', script_name+'.o', script_name))
        
        # Monitor jobs 
        give_feedback('Waiting for jobs to complete...')
        jobs_remaining = len(scripts)
        wait_time = 60
        start_time = time.time()
        while jobs_remaining > 0:
            give_feedback(str(jobs_remaining)+' jobs remaining; waiting '+str(wait_time)+' seconds...')
            time.sleep(wait_time)
            jobs_remaining = int(commands.getoutput('qstat | grep script_num | wc -l'))
            
        give_feedback('All jobs finished after %.1f seconds'%(finish_time - start_time))
            
    else:
        for n, script in enumerate(scripts):
            print 'Executing script %d of %d...'%(n+1,len(scripts))
            os.system(script)

    # Announce job completion
    finish_time = time.time()
    give_feedback('Done.\n')

###################################################################################

# Split reads into files containing 100K reads each
@files([read1_file_glob, read2_file_glob],[read1_split_file_glob, read2_split_file_glob])
def stage_one(ins, outs):
    os.system('rm -r %s'%output_dir)
    os.system('mkdir %s'%output_dir)
    os.system('rm -r %s'%reads_dir)
    os.system('mkdir %s'%reads_dir)
    
    # Get list of r1 and r2 files
    read1_files = glob.glob(read1_file_glob)
    read2_files = glob.glob(read2_file_glob)
    
    # Make sure there are the same number of such files
    assert(len(read1_files) == len(read2_files))
    num_paired_read_files = len(read1_files)
    assert num_paired_read_files > 0
    
    scripts = []
    for pair_num in range(num_paired_read_files):
        read1_file_name = read1_files[pair_num]
        read2_file_name = read2_files[pair_num]
        
        read1_split_file_prefix = read1_split_file_glob[:-2]+'.%03d'%pair_num
        read2_split_file_prefix = read2_split_file_glob[:-2]+'.%03d'%pair_num
        
        give_feedback('Splitting fastq files...\n')   
        script = '''
        source ~/.bash_profile  
        split -l 400000 %s %s.
        split -l 400000 %s %s.
        ''' % (read1_file_name, read1_split_file_prefix, read2_file_name, read2_split_file_prefix) 
        scripts.append(script)
    
    # Submit jobs    
    submit_and_complete_jobs(scripts, use_multiple_nodes)

# Create directories for experiments 
@follows(stage_one)
@files(experiments_file, '%s/*/*'%results_dir)
def stage_two(ins, outs):

    # Load experiment information
    f = open(experiments_file,'r')
    line = f.readline()
    atoms = line.strip().split()
    timepoints = atoms[2:]
    os.system('rm -r %s'%results_dir)
    os.system('mkdir %s'%results_dir)
    for line in f.readlines():
        if len(line.strip())==0:
            continue;
        atoms = line.strip().split()
        experiment = atoms[0]
        os.system('mkdir %s/%s'%(results_dir,experiment))
        for tp in timepoints:
            os.system('mkdir %s/%s/%s'%(results_dir,experiment,tp))  

    give_feedback('Done!\n')

# Process reads in each fastq file into observed_seqs
@follows(stage_two)
@files([read1_split_file_glob, read2_split_file_glob],observed_seq_files_split)
def stage_three(ins,outs):
    give_feedback('Reconstructing observed_seqs from reads...\n')

    # Get list of extensions
    extensions = ['.'.join(x.split('.')[-2:]) for x in glob.glob(read1_split_file_glob)]
    
    # For each extension, farm out a parse_seqs.py job
    scripts = []
    for extension in extensions:
        r1_file = '%s.%s'%(read1_split_file_glob[:-2],extension)
        r2_file = '%s.%s'%(read2_split_file_glob[:-2],extension)
        command = '%s %s %s %s %s %s %s'%(parse_seqs_script,r1_file, r2_file, barcodes_file, regions_file, experiments_file, results_dir)
        script = '''
        source ~/.bash_profile
        %s
        '''%command
        scripts.append(script)
    submit_and_complete_jobs(scripts, use_multiple_nodes)

# Combine observed seqs at each timepoint into a single file, then count the number
# of occurances of each one
@follows(stage_three)
@files(observed_seq_files_split, [observed_seq_files_by_sample, unique_seq_files_by_sample])
def stage_four(ins,outs):
    give_feedback('Computing results/[experiments]/[timepoints]/[observed_seqs, unique_seqs]...\n')
    sample_dirs = glob.glob('%s/*/*'%results_dir)
    scripts = []
    for sample_dir in sample_dirs:
        script = '''
        source ~/.bash_profile
        cat %s/observed_seqs.* > %s/observed_seqs
        cat %s/observed_seqs | sort | uniq -c | sort -nr | nl -b t > %s/unique_seqs
        '''%(sample_dir, sample_dir, sample_dir, sample_dir)
        scripts.append(script)
    submit_and_complete_jobs(scripts, use_multiple_nodes)

# Combine observed seqs at all timepoints into a single file, then count the number
# of occurances of each one
@follows(stage_four)
@files(observed_seq_files_by_sample, [observed_seq_files_by_experiment, unique_seq_files_by_experiment])
def stage_five(ins,outs):
    give_feedback('Computing results/[experiments]/[observed_seqs, unique_seqs]...\n')
    experiment_dirs = glob.glob('%s/*'%results_dir)
    scripts = []
    for experiment_dir in experiment_dirs:
        script = '''
        source ~/.bash_profile
        cat %s/*/observed_seqs > %s/observed_seqs 
        cat %s/observed_seqs | sort | uniq -c | sort -nr | nl -b t > %s/unique_seqs
        '''%(experiment_dir, experiment_dir, experiment_dir, experiment_dir)
        scripts.append(script)
    submit_and_complete_jobs(scripts, use_multiple_nodes)

# Tally seqs
@follows(stage_five)
@files([unique_seq_files_by_experiment, unique_seq_files_by_sample], tallied_seq_files)
def stage_six(ins,outs):
    give_feedback('Tallying seq occurances at each timepoint...\n')
    experiment_dirs = glob.glob('%s/*'%results_dir)
    scripts = []
    for experiment_dir in experiment_dirs:
        script = '''
        source ~/.bash_profile
        %s %s %s
        '''%(tally_seqs_script, experiment_dir,experiments_file)
        scripts.append(script)
    submit_and_complete_jobs(scripts, use_multiple_nodes)

# Summarize seqs
@follows(stage_six)
@files(tallied_seq_files, summarized_seq_files)
def stage_seven(ins,outs):
    give_feedback('Summarizing sequences observed in each experiment...\n')
    experiment_dirs = glob.glob('%s/*'%results_dir)
    scripts = []
    for experiment_dir in experiment_dirs:
        script = '''
        source ~/.bash_profile
        %s %s %s %s
        '''%(summarize_seqs_script, regions_file, experiments_file, experiment_dir)
        scripts.append(script)    
    submit_and_complete_jobs(scripts, use_multiple_nodes)

# Collect summarized_seqs
@follows(stage_seven)
@files(summarized_seq_files, '%s/*.txt'%summary_dir)
def stage_eight(ins,outs):
    give_feedback('Collecting summarized sequences...\n')
    os.system('rm -r %s'%summary_dir)
    os.system('mkdir %s'%summary_dir)
    experiment_dirs = glob.glob('%s/*'%results_dir)
    for experiment_dir in experiment_dirs:
        experiment = experiment_dir.split('/')[-1]
        command = 'cp %s/summarized_seqs %s/%s_summary.txt'%(experiment_dir,summary_dir,experiment)
        os.system(command)
 
# Compute CRISPR efficiency
@follows(stage_eight)
@files('%s/*.txt'%summary_dir, efficiency_file)
def stage_nine(ins,outs):
    give_feedback('Estimating CRISPR efficiency\n')

    script = '''
    source ~/.bash_profile
    %s %s %s %s
    '''%(compute_efficiencies_script, summary_dir, infection_file, efficiency_file)
    #submit_and_complete_jobs([script])
    os.system(script)
    
# Make prevalence plots
@follows(stage_nine)
@files('%s/*.txt'%summary_dir, '%s/*.pdf'%plots_dir)
def stage_ten(ins,outs):
    give_feedback('Plotting allele prevalence over time...\n')
    # Clean out plots dir
    os.system('rm -r %s'%plots_dir)
    os.system('mkdir %s'%plots_dir)
    
    # Plot prevalences for every summary file
    summary_file_names = glob.glob(summary_dir + '/*.txt')
    scripts = []
    for summary_file_name in summary_file_names:
        print 'Processing %s...'%summary_file_name

        # Make plots showing individual lines and area
        script_lines = '''
        source ~/.bash_profile 
        %s %s %s lines
        '''%(plot_prevalences_script, summary_file_name, plots_dir)
        scripts.append(script_lines)

        script_area = '''
        source ~/.bash_profile 
        %s %s %s area
        '''%(plot_prevalences_script, summary_file_name, plots_dir)
        scripts.append(script_area)

    # For some reason, the nodes take forever to do this. 
    #submit_and_complete_jobs(scripts)
    # Instead, just run each script here
    for script in scripts:
        os.system(script)

    give_feedback('Done.\n')


### Run pipeline
print '##################################################'
print 'Begin read mapping pipeline'
pipeline_start_time = time.time()
pipeline_run([stage_ten])
interval = time.time() - pipeline_start_time
print '\nPipeline Done! Runtime: %.1f min' % (interval/60.0)
print '##################################################'


    
    

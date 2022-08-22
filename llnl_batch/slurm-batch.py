#! /usr/bin/env python2.7
from __future__ import print_function
import os
import shutil
import argparse
import subprocess

# -------------------------- ARGUMENT PARSING ------------------------ #

prog_description = 'Creates files needed for slum jobs and optionally submits jobs'
prog_epilogue    = 'And that\'s all, folks!'

# Initialize argument parser
parser = argparse.ArgumentParser(
    description=prog_description,
    usage='%(prog)s [options]',
    epilog=prog_epilogue,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)


# Arguments 
#
parser.add_argument('-d', '--dir',type=str,
                    default='/p/lscratchh/ftorales/AI-codesign/',
                    help='top level scratch directory')
parser.add_argument('-j', '--jobset',type=str,default='codesign-example',
                    help='an identifier for the jobset, jobs submitted from a subdir of this name')
parser.add_argument('-n', '--njob', type=int,default='10',
                    help='number of jobs to submit')
parser.add_argument('-s', '--submit',action='store_true',
                    help='submit jobs with sbatch after creating')
parser.add_argument('-m', '--merge',action='store_true',
                    help='submit separate merge job with dependencies on this one')
parser.add_argument('-v', '--verbose',action='store_true',
                    help='turn on diagnostic print statements')
parser.add_argument('-ne', '--nevents',type=str,default='1000',
                    help='number of events to run per subjob')
parser.add_argument('-p', '--particle',type=str,default='pion+',
                    help='particle species: pion+, pion-, pion0, muon+, muon-, etc')
parser.add_argument('-pmin', '--pmin',type=str,default='0.0',
                    help='minimum generated momentum')
parser.add_argument('-pmax', '--pmax',type=str,default='100.0',
                    help='maximum generated momentum')


# parser.add_argument('-c', '--config',type=str,default='configuration-1',
#                     help='base for config files')

# Parse arguments
args = parser.parse_args()

if not os.path.isdir(args.dir):
    raise ValueError('directory {} not found.'.format(args.dir))

#parameters for this example
bank='mlodd'
# subjob_script='run-subjob.sh'
subjob_script='run_gun.sh'
merge_script='run-merge.sh'
# time_est_subjob='00:10:00' #ten hours
# time_est_merge='00:02:00' #two hours
time_est_subjob='00:00:10' #ten minutes
time_est_merge='00:00:02' #two minutes


# Further example of how to make sure enviroment is set
# if os.getenv('HIPLOCAL'):
#     hip_local = os.getenv('HIPLOCAL')
# else:
#     raise ValueError('You must have $HIPLOCAL defined to run this script.')


#create the submission directory args.dir/args.jobset
submit_dir = os.path.join(args.dir,args.jobset)
if (args.verbose):
    print('Submission directory: ', submit_dir)

if not os.path.isdir(submit_dir):
    os.mkdir(submit_dir) 
    os.chmod(submit_dir,0o2770)

#also create subdirectory for subjob log files
logdir = os.path.join(submit_dir,'log')
if not os.path.isdir(logdir):
    os.mkdir(logdir) 
    os.chmod(logdir,0o2770) 

#...and output
outputdir = os.path.join(submit_dir,'output')
if not os.path.isdir(outputdir):
    os.mkdir(outputdir) 
    os.chmod(outputdir,0o2770) 

#copy scripts to submission directory
#could alternatively put them in your path
extra_files=[subjob_script,merge_script]
for script in extra_files:
    script_src = os.path.join("./",script)
    script_dest = os.path.join(submit_dir,script)
    shutil.copyfile(script_src,script_dest)
    shutil.copystat(script_src,script_dest)


#Further example of how to copy a (possibly job-specific) configuration file
#in this case the config file is located in the same directory you execute this script from 
#and has a suffix of .config
# configfile = args.config + '.config'
# configfile_src = os.path.join("./",configfile)
# configfile_dest = os.path.join(submit_dir,configfile)
# shutil.copyfile(configfile_src,configfile_dest)
# shutil.copystat(configfile_src,configfile_dest)


#Create a job script for slurm
basefile = os.path.join(submit_dir,'job')
batchfile = basefile + '.sh'
f = open(batchfile,'w')
f.write('#!/bin/bash\n')
f.write('#SBATCH -n 1\n')
f.write('#SBATCH -t %s\n' % time_est_subjob)
f.write('#SBATCH --job-name=%s\n' % args.jobset)
f.write('#SBATCH -p pbatch\n')
f.write('#SBATCH -A %s\n' % bank)
f.write('#SBATCH --array=0-%d\n' % (args.njob-1))
#the %A's are slurm syntax to insert the numeric job ID and subjob index
f.write('#SBATCH -o %s%s.%s.out\n' % (os.path.join(submit_dir,'log/'),'%A','%a')) 
f.write('#SBATCH -e %s%s.%s.err\n' % (os.path.join(submit_dir,'log/'),'%A','%a'))
f.write('%s/%s -n %s -p %s -j %s --pmin %s --pmax  %s' %(submit_dir,subjob_script,
    args.nevents,args.particle,args.jobset,args.pmin,args.pmax))
# ^This line parses all the arguments, and interfaces to the subjob script

f.close()

#Further example, replace the last f.write(...) with the following to additionally pass config to subjob_script
#f.write('%s/%s %s %s' % (submit_dir,subjob_script,submit_dir,args.config))




#submit job
if (args.submit):
    submit_command = 'sbatch ' + batchfile
    print (submit_command)
    os.system(submit_command)
else:
    print('%s/%s -n %s -p %s -j %s --pmin %s --pmax  %s -t test' % (submit_dir,subjob_script,args.nevents,args.particle,args.jobset,args.pmin,args.pmax))
    os.system('%s/%s -n %s -p %s -j %s --pmin %s --pmax  %s -t test' % (submit_dir,subjob_script,args.nevents,args.particle,args.jobset,args.pmin,args.pmax))

#Now optionally create job script for merge job and possibly submit
if ( args.merge ):
    batchfile = basefile + '.merge.sh'
    f = open(batchfile,'w')
    f.write('#!/bin/bash\n')
    f.write('#SBATCH -n 1\n')
    f.write('#SBATCH -t %s\n' % time_est_merge)
    f.write('#SBATCH --job-name=%s\n' % args.jobset)
    f.write('#SBATCH -p pbatch\n')
    f.write('#SBATCH -A %s\n' % bank)
    f.write('#SBATCH -o %smerge.out\n' % os.path.join(submit_dir,'log/'))
    f.write('#SBATCH -e %smerge.err\n' % os.path.join(submit_dir,'log/'))
    f.write("#SBATCH --dependency=singleton\n")
    f.write('%s/%s %s' % (submit_dir,merge_script,submit_dir))
    f.close()
    if (args.submit):
        submit_command = 'sbatch ' + batchfile
        print (submit_command)
        os.system(submit_command)

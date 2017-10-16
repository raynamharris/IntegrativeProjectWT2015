# Cutadapt: Filter and Trim Low Quality Reads

Lets remove adapters (trim) and low quality reads (filter) with the program cutadapt.

## Navigate to your raw data directory on stampede, 

~~~ {.bash}
ssh <username>@stampede.tacc.utexas.edu
cd $SCRATCH/JA16444/00_rawdata
~~~


## Using cutadapt on TACC

The cutadapt program is not maintained by TACC. We need to use the program stored in a colleagues working directory. To access the program, set the path. 

~~~ {.bash}
PATH=/work/01184/daras/bin/cutadapt-1.3/bin:$PATH
~~~

## Write a cutadapt commands file 

This loop will create a command to process the paired-end reads. 

~~~ {.bash}
rm 02_filtrimmedreads.cmds # precaution just incase you've done this bevore
for R1 in *R1_001.fastq.gz
do
    R2=$(basename $R1 R1_001.fastq.gz)R2_001.fastq.gz
    R1filtrim=$(basename $R1 fastq.gz)filtrim.fastq.gz
    R2filtrim=$(basename $R2 fastq.gz)filtrim.fastq.gz
    echo $R1 $R2 $R1filtrim $R2filtrim
    echo "cutadapt -q 15,10 -a GATCGGAAGAGCACACGTCTGAACTCCA -A ATCGTCGGACTGTAGAACTCTGAACGTG -m 22 -o $R1filtrim -p $R2filtrim $R1 $R2" >> 02_filtrimmedreads.cmds
done
~~~

### Option 1: Submit a job on Stampede.
Create a launcher script and launch the job

~~~ {.bash}
launcher_creator.py -t 4:00:00 -n 02_filtrimmedreads -j 02_filtrimmedreads.cmds -l 02_filtrimmedreads.slurm -A NeuroEthoEvoDevo -q 'normal'  
sbatch 02_filtrimmedreads.slurm
~~~

### Option 2: Use an interactive compute node
Request compute time, makde cmd file executable, run commands.

~~~ {.bash}
idev -m 120
chmod a+x 02_filtrimmedreads.cmds
bash 02_filtrimmedreads.cmds
~~~

## Clean up

Now, let's make our processed reads read only so we don't accidentally modify them. 

~~~ {.bash}
chmod a-w *filtrim.fastq.gz 
~~~

Now, move the processed reads to a new file.

~~~ {.bash}
mkdir ../02_filtrimmedreads
mv *filtrim.fastq.gz ../02_filtrimmedreads
~~~


## References
- Cutadapt: http://cutadapt.readthedocs.io/en/stable/guide.html
- BioITeam Launcher Creator: https://wikis.utexas.edu/display/bioiteam/launcher_creator.py
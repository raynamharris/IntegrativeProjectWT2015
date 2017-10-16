# Fastqc: Quality Control Assessment of Processed Reads

Now, let's see if filtering and trimming improved our read quality.

## Navigate to your working directory on stampede, 

~~~ {.bash}
cd $SCRATCH/JA16444/02_filtrimmedreads
~~~

## Write a fastqc commands file 

Create a commands file using a for loop. This will create a file with one command per line for performing the fastqc function on every read in the directory. We use the `>>` function to append the new line to the existing file. In case we made an error and are rerunning the loop, its always good to start with the `rm` command, just incase there is a bad version of this file around.

~~~ {.bash}
rm 03_fastqc.cmds 
for file in *.fastq.gz
do
     echo $file
     echo "fastqc $file" >> 03_fastqc.cmds
done
~~~

Check to see that the commands file looks like it should

~~~ {.bash}
cat 03_fastqc.cmds
~~~

### Option 1: Submit a job on Stampede.
Create a launcher script and launch the fastqc job

~~~ {.bash}
launcher_creator.py -t 0:30:00 -n 03_fastqc -j 03_fastqc.cmds -l 03_fastqc.slurm -A NeuroEthoEvoDevo -m 'module load fastqc/0.11.5'
sbatch 03_fastqc.slurm
~~~

### Option 2: Use an interactive compute node
Request compute time, load module, make cmd file executable, run commands.

~~~ {.bash}
idev -m 120
module load fastqc/0.11.5
chmod a+x 03_fastqc.cmds
bash 03_fastqc.cmds
~~~

## Clean up

Then, I moved all the output files to a separate folder where we will store the fastqc results.

~~~ {.bash}
mkdir ../03_fastqc
mv *.html ../03_fastqc
mv *.zip ../03_fastqc
~~~

## save locallaly

One must save the data locally in order to view the html files. 

In a new terminal window:

~~~ {.bash}
cd <pathtoplaceonyourpersonalpc>
cd $RNAseqProject/$RNAseqJob
scp <username>@stampede.tacc.utexas.edu:$SCRATCH/$RNAseqProject/$RNAseqJob/03_fastqc/*html .
~~~

## References
FastQC: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
BioITeam Launcher Creator: https://wikis.utexas.edu/display/bioiteam/launcher_creator.py
FastQC Overview: https://wikis.utexas.edu/display/bioiteam/FASTQ+Quality+Assurance+Tools

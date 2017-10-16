# Fastqc: Quality Control Assessment of Raw Reads

We must check the quality of our reads to see if we need to do any processing to improve the quality. We will use the function called `fastqc` (aka quality control of fastq files).

## Login to stampede and navigate to the rawdata

~~~ {.bash}
ssh <username>@stampede.tacc.utexas.edu
cd $SCRATCH/JA16444/00_rawdata
~~~

## Write a fastqc commands file 

Create a commands file using a for loop. This will create a file with one command per line for performing the fastqc function on every read in the directory. We use the `>>` function to append the new line to the existing file. In case we made an error and are rerunning the loop, its always good to start with the `rm` command, just incase there is a bad version of this file around.

~~~ {.bash}
for file in *.fastq.gz
do
     echo $file
     echo "fastqc $file" >> 01_fastqc.cmds
done
~~~

Check to see that the commands file looks like it should

~~~ {.bash}
cat 01_fastqc.cmds
~~~

### Option 1: Submit a job on Stampede.
Create a launcher script and launch the fastqc job

~~~ {.bash}
launcher_creator.py -t 0:30:00 -n 01_fastqc -j 01_fastqc.cmds -l 01_fastqc.slurm -q normal -m 'module load fastqc/0.11.5'
sbatch 01_fastqc.slurm
~~~

### Option 2: Use an interactive compute node
Request compute time, makde cmd file executable, load modules, run commands.

~~~ {.bash}
idev -m 120
module load fastqc/0.11.5
chmod a+x 01_fastq.cmds
bash 01_fastq.cmds
~~~


## Exploring the raw fastq data

First, let's move the output files to a separate folder where we will store and view the fastqc results.

~~~ {.bash}
mkdir ../01_fastqc
mv *_fastqc ../01_fastqc
cd ../01_fastqc
~~~

To view the text files from the fastq output, first unzip the zip file

~~~ {.bash}
# unzip all the fastqc files 
for file in *.zip
do
unzip $file
done
~~~

Now, use a one-liner to extract the read lenght and read counts for each file

~~~ {.bash}
## get read lenght and count for each file
for file in *R1*fastqc
do
echo $file
cat $file/fastqc_data.txt | grep -w -A 1 "Length" | grep -v "Length"
done
~~~

You can also save this output to a file that can be downloaded and imported into R or python for plotting and stats

~~~ {.bash}
## save read lenght and count to a new file
for file in *R1*fastqc
do
cat $file/fastqc_data.txt | grep -w -A 1 "Length" | grep -v "Length" >> readcounts.txt
done
~~~

## References
- FastQC: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
- BioITeam Launcher Creator: https://wikis.utexas.edu/display/bioiteam/launcher_creator.py
- FastQC Overview: https://wikis.utexas.edu/display/bioiteam/FASTQ+Quality+Assurance+Tools

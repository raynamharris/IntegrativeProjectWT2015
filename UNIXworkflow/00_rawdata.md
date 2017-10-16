# Get, Store, and/or Retrieve Raw Data

Includes the process and commands needed to download the data to scratch on Stampede with `00_gsaf_download`. In theory, this only has to be done once, but it could be done again if data are lost or compromised. 


## 00_gsaf_download 

The first step in the bioinformatics pipeline is to download the data from the sequencing facility, in this case GSAF. The GSAF provides the script for downloading the data from there amazon web server. 

1. Setup project and data directories on Stampede.
2. Copy script to download data: Save a script to download the data (called "00_gsaf_download.sh") in each subdirector  
3. Download the data using TACC: For this three-part step I created a commands file (called "00_gsaf_download.cmd") with the command to execute the download script,  created a launcher script to execute the commands file to execute the launcher script (I know, sounds like a lot of steps, but this is so I use TACC's compute power not my own), and launched the job on TACC  
4. Repeat steps 2 and 3 for all RNAseq jobs 

### Setup project and RNAseq job directories 

Login to TACC using ssh with your password and authentication credentials. Replace "<username>" with your TACC user name. 

~~~ {.bash}
ssh <username>@stampede.tacc.utexas.edu
~~~

On scratch, create the project directory using the job number (e.g. JA16444) and subsubdirectory called 00_rawdata. The argument `-p` will create the parent and subdirectories if they do not already exist.

~~~ {.bash}
mkdir -p $SCRATCH/JA16444/00_rawdata
cd $SCRATCH/JA16444/00_rawdata
~~~

### Copy script to download data 
Copy the text for the gsaf download script found here: https://wikis.utexas.edu/display/GSAF/How+to+download+your+data. Use nano to save it as `00_gsaf_download.sh`. Then, make the file executable.

~~~ {.bash}
chmod a+x 00_gsaf_download.sh
~~~


### Download the data using TACC
Create a luncher scipt to execute the data on TACC. To so, use nano t create a new file `00_gsaf_download.cmds` with the following line of text `00_gsaf_download.sh "<amazonwebaddress>"

Now, I use `launcher_creator.py` to create a launcher script that will tell how to launch this job on TACC. The arguments are defined clearly on this website: https://wikis.utexas.edu/display/bioiteam/launcher_creator.py. Then I will use `sbatch 00_gsaf_download.slurm` to launch the job.

~~~ {.bash}
launcher_creator.py -t 12:00:00 -n 00_gsaf_download -j 00_gsaf_download.cmds -l 00_gsaf_download.slurm -q normal
sbatch 00_gsaf_download.slurm
~~~

### Make all fastq.gz files READ ONLY!

To protect our data from accidental deletion or overriding, let's make all the files read only!

~~~ {.bash}
chmod a-w *fastq.gz
ls -l
~~~

### Clean up a little

Create a new directory called `wget_log` and move all the output files associated with downloading the data.

~~~ {.bash}
mkdir wget_log
mv *.wget.log wget_log
ls
~~~

### Repeat for all jobs. 
The great thing about using TACC for this is that you can go about doing other things while the files are download. Reset the environment variables above and repeat the process for all other RNAseq jobs.

### Save summary file for metadata for later

Samples often get named different things along the way. This is an attempt to make a rossetta stone for future use. 

~~~ {.bash}
for R1 in *R1_001.fastq.gz
do   
	R2=$(basename $R1 R1_001.fastq.gz)R2_001.fastq.gz
	echo "$R1" | awk -F '_' '{print "15-" $1 "," "15" $1 "," $2 "," $1 "_" $2 "_" $3 ","  $R1 }' >> 00_rossettastone.csv
done
~~~

NCBI

For NCBI, I need to send them a zipped folder of all my fastqc and mdsums files. 

For workflow purposes, I keep launcher scripts and output files in the directory where they are run and created. I know its against best practices, but the hierarchical organization schemes that I have tried have failed me. 

So, first I'll move these extranous files to a temporary folder. Then, I'll zip the raw files. Then I move the extracous files back into the un-zipped director of raw data.

~~~{.bash
mv 00_rawdata/*.slurm temp/
mv 00_rawdata/*.e8* temp/
mv 00_rawdata/*.o8* temp/
mv 00_rawdata/*.cmds temp/
mv 00_rawdata/*.sh temp/
mv 00_rawdata/multiqc* temp/
~~~

Then zip

~~~{.bash
tar -cvpf 00_rawdata.tar 00_rawdata
~~~

Then turn the extra files to their home. 

~~~{.bash
mv temp/* 00_rawdata
rm temp
~~~

Then I'll do the same for kallisto. FIrst, let's rename the files so that they have the sample name in the file name. 

~~~{.bash
#cd 02_kallistoquant
#for abund in */abundance.tsv; do
#    newAbund=$(echo $abund | perl -pe 's|(.*)/.*|\1/abundance_\1.tsv|')
#    mv $abund $newAbund
#done
for abund in */abundance.h5; do
    newAbund=$(echo $abund | perl -pe 's|(.*)/.*|\1/abundance_\1.h5|')
    mv $abund $newAbund
done
for abund in */run_info.json; do
    newAbund=$(echo $abund | perl -pe 's|(.*)/.*|\1/runinfo_\1.json|')
    mv $abund $newAbund
done
~~~

Now, lets move the h5 and jason files.

~~~{.bash
mkdir ../hd5runinfo
mv */*h5 ../hd5runinfo
mv */*json ../hd5runinfo
~~~

~~~{.bash
cd ..
tar -cvpf 02_kallistoquant.tar 02_kallistoquant
~~~

Now use ftp to transfer the files (00_rawdata.tar, 02_kallistoquant.tar, and FMR1MetaData.xls) to NCBI. First, I copy the files to the Hofmann lab server so I can use Filzilla to transfer the files. The instructions are at: https://www.ncbi.nlm.nih.gov/geo/info/seq.html#data

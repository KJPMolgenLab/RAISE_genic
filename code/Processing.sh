#!/bin/bash
#SBATCH --job-name=TrimmFQ
#SBATCH --partition=fuchs
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1200
#SBATCH --time=00:30:00
#SBATCH --mail-type=ALL 


################################
#		               #
#	Preprocessing & QC     #
#			       #
################################

#QC was performed in FastQC 
#Filtering was performed in trimmomatic

homedir=$(pwd)	   
datadir=$1

cd $datadir

if [ 1 == 1 ]
then
    
   #Read filter
   n=0

   for filename in $(realpath *.fq.gz)
   do

       # first, make the base by removing fastq.gz
       base=$(basename $filename .fq.gz)
       #echo $base
       java -jar ~/Software/Trimmomatic-0.39/trimmomatic-0.39.jar \
	    SE ${base}.fq.gz \
	    ${base}.qc.fq.gz \
            ILLUMINACLIP:/home/fuchs/agmisc/chiocchetti/Software/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 \
            LEADING:3 TRAILING:3 \
            SLIDINGWINDOW:4:15 \
            MINLEN:36 &
       n=$(($n+1))
       if [ $(($n % 20)) == 0 ]
       then
	   wait
       fi
       
       
   done

   wait
fi

R --slave --no-save --args $datadir < ~/Projects/RAISE_genic/code/01_Mapping_Alignment.R 

wait

cd $homedir

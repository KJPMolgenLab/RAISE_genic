#!/bin/bash -x
#SBATCH --job-name=FeatNet
#SBATCH --partition=fuchs
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --mem-per-cpu=6000
#SBATCH --cpus-per-task=1
#SBATCH --time=20:00:00
#SBATCH --mail-type=ALL 


cd /home/fuchs/agchiocchetti/chiocchetti/Projects/RAISE_genic/
inputfile="/scratch/fuchs/agchiocchetti/public/data/ASC_Data_2022/ccdg_asc_ndd_daly_talkowski_goethe_asd_exome/LGD_ccdg_asc_ndd_daly_talkowski_goethe_asd_exome.vcf.gz.raw" 
inputbim= "/scratch/fuchs/agchiocchetti/public/data/ASC_Data_2022/ccdg_asc_ndd_daly_talkowski_goethe_asd_exome/LGD_ccdg_asc_ndd_daly_talkowski_goethe_asd_exome.vcf.gz.bim"


linestot=$(< $inputfile wc -l)
linespart=$(expr $linestot / 10)

tail -n +2 $inputfile | split -l $linespart - $inputfile.part_

    
for file in $inputfile.part_*
do
    filename=$(basename -- "$file") &&
    extension="${filename##*.}" &&
    echo $extension &&
    head -n 1 $inputfile > "$file.tmp_file" &&
    cat $file >> $file.tmp_file &&
    mv -f $file.tmp_file "$file" &&
    Rscript --slave --no-save --no-restore code/06_2_Networkfeatures_only.R $file $inputbim $extension &
done

wait



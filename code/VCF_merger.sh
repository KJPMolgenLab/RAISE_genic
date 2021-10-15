#!/bin/bash
#SBATCH --job-name=merger
#SBATCH --partition=fuchs
#SBATCH --nodes=2
#SBATCH --ntasks=20
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1200   
#SBATCH --time=04:00:00
#SBATCH --no-requeue
#SBATCH --mail-type=ALL


## estimated run time ~1.5h


vcffolder="/scratch/fuchs/agmisc/chiocchetti/RAISEGenicData/VCF_Files/"
refdb="/scratch/fuchs/agmic/chiocchetti/annovar/humandb/"
homedir=$(pwd)

cd $vcffolder

echo "merge Brussel files"

if [ ! -f "Brussels/Brussels_merged.vcf.gz" ]
then    
    cd "$vcffolder/Brussels"
    echo "compress files in Brussels"
    # get uncompressed files 
    compressfiles=$(find . -type f | grep .vcf$)
    for file in $compressfiles
    do
	#compress each 
	echo $file
	bgzip -c $file > $file.gz &
    done
    wait
    
    echo "index files in Brussels"
    for file in $compressfiles
    do
	#create index 
	echo "indexing $file.gz"
	bcftools index $file.gz &
    done
    wait
    
    echo "merge files in Brussels"
    files=$(find . -type f | grep .vcf.gz$)
    bcftools merge  $files --threads 20 -0 -Oz -o Brussels_merged.vcf.gz
    bcftools index Brussels_merged.vcf.gz
    wait
    
    echo "done"
else
    echo "Brussels file merged exits nothing done"
fi


echo "index Frankfurt files"

if [ ! -f "DEUPUM_RAISE-GENIC_210526.vcf.gz.csi" ]
then
    cd "$vcffolder/Frankfurt"
    files=$(find . -type f | grep .vcf.gz$)
    bcftools index $files --threads 20
    echo "done"
else
    echo "Frankfurt index files exist,  skipped"
fi


echo "index Kuopio files"

if [ ! -f "AnVIL_CCDG_Broad_NP_Epilepsy_FINKPH_EPIL_CO_MORBIDI_MDS_WES_y2-only-raisegenic.vcf.gz.csi" ]
then    
    cd "$vcffolder/Kuopio"
    files=$(find . -type f | grep .vcf.gz$)
    bcftools index $files threads 20
    echo "done"
else
    echo "Kuopio index files exist, skipped"
fi


cd $vcffolder

echo "convert SNPfiles to plinkables and perform quality control "

FrankfurtSNP=Frankfurt/Epi25_DEUPUM_GSA-MD_Year1_Year2_v1.vcf.gz

# define function for internal quality control

function qc(){
    if [ ! -f $1 ]
    then
	echo "$1 soes not exist"
	exit
    fi
    
    #filter 95% genotyping sample  95%SNP genotyping, 2% minor allele frequency and HWE cutoffs
    echo "$1 passing through quality control"
    homedir=$(pwd)
    outdir=$(dirname $1)
    filebase=$(basename $1)
    
    nohup plink --vcf $1 --geno 0.05 --mind 0.05 --hwe 10e-8 --maf 0.02\
	  --make-bed --out $outdir/QC1 > QC1_output.log 

    cd $outdir

    echo "$1 check gender"
    nohup plink --bfile QC1 --check-sex --out QC2_Sexcheck > QC2_Plinkoutput.log

    # count errors 
    grep "PROBLEM" QC2_Sexcheck.sexcheck | awk '$3 == 0' > QC2_NAGender.txt
    grep "PROBLEM" QC2_Sexcheck.sexcheck | awk '$3 != 0' > QC2_BadGender.txt

    read NAgend filename <<< $(wc -l QC2_NAGender.txt)
    read Badgend filename <<< $(wc -l QC2_BadGender.txt)

    echo "$1 update missing gender"
    if [ $NAgend -gt 0 ]
    then
	awk {'print($1 "\t"  $2 "\t" $4)'} QC2_NAGender.txt > temp;
	plink --bfile QC1 --update-sex temp --make-bed --out QC2 > QC3_Plinkoutput_updatesex.log;
	rm temp;
	echo ".....gender updated"
    else
	rename QC1 QC2 QC1.*;  
    fi;

    # exclude wrong gender annotated 
    if [ $Badgend -gt 0 ]
    then
	awk {'print($1 "\t"  $2)'} QC2_BadGender.txt > QC2_BadGender_removed.txt;
	echo ".....Samples with genderincongruencies being removed";
	plink --bfile QC2 --remove QC2_BadGender_removed.txt --make-bed --out QC3 > QC3_Plinkoutput_removebadgender.log;
    else
	echo ".....no gender discrepancies detected"
	rename QC2 QC3 QC2.*
    fi

    echo "$1 check inbreeding and contamination"
    nohup plink --bfile QC3 --het --out Breeding_output > QC4_Plinkoutput_breeding.log
    awk '$6>0.2'  Breeding_output.het | awk '$6!="F"'  > QC4_inbredIDs.txt
    awk '$6<-0.15'  Breeding_output.het > QC4_contaminatedIDs.txt

    read cont x  <<< $(wc -l QC4_contaminatedIDs.txt)
    read inbrd x <<< $(wc -l QC4_inbredIDs.txt)

    echo "$1 exclude contaminated samples "
    if [ $cont -gt 0 ]
    then
    awk {'print $1 "\t" $2'} QC4_contaminatedIDs.txt > QC4_Contaminated_samples_removed.txt
    plink --bfile QC3 --remove Contaminated_samples_removed --make-bed --out QC4 > QC4_Plinkoutput_rm_cont.txt
    else
	echo ".....no samples needed to be removed due to contamination"
	rename QC3 QC4 QC3.*;
    fi

    echo "$1 remove inbreeding individuals"
    if [ $inbrd -gt 0 ]
    then
	awk {'print $1 " " $2'} QC4_inbredIDs.txt > QC4_Inbred_samples_removed.txt
	plink --bfile QC4 --remove QC4_Inbred_samples_removed.txt --make-bed --out QC5 > QC5_Plinkoutput_rmInbred.log
    else
	echo "$1 no samples needed to be removed due to inbreeding"
	rename QC4 QC5 QC4.*
    fi

    echo"$1 cyclic pedigirees are not tested as no family data"
    rename QC5 QC6 QC5.*

    echo "$1 relatedness check"
    nohup plink --bfile QC6  --genome --out QC6_Relcheck > QC6_Plinkoutput_relcheck.log
    awk '$9 >0.8 && $10 >0.8'  QC6_Relcheck.genome > QC6_Relcheck.duplicates

    read dups x <<< $(wc -l QC6_Relcheck.duplicates)

    if [ $dups -gt 1 ]
    then
       echo "$1 remove duplicated"
       awk {'print $1 " " $2'} QC6_Relcheck.duplicates > QC6_duplicatesDeleted.txt
       nohup plink --bfile QC6 --remove QC6_duplicatesDeleted.txt --make-bed --out QC7 > QC6_Plinkoutput_rmDupls.log
    else
	echo "no duplicated samples detected"
	rename QC6 QC7 QC6.*
    fi
    
    echo "$1 checked for cross relatedness nothing done at the moment but written out"
    nohup plink --bfile QC7 --genome --out QC7_Relcheck > QC7_Plinkoutput_Relcheckround2.log
    awk '($1!=$3 && $10>0.1)' QC7_Relcheck.genome > QC7_CrossfamRelations.txt

    rename QC7 $filebase.qc QC7.*
    
    cd $homedir
}

qc $FrankfurtSNP


# qc $BrusselsSNP
# qc $KuopioSNP
# plink --bfile $FrankfurtSNPs.qc --bmerge $BrusselsSNPs.qc $KuopioSNPs.qc --make-bed --out RAISE_Genic_merged.qc


echo "merge WES data"
Brussel=Brussels/Brussels_merged.vcf.gz
Frankfurt=Frankfurt/DEUPUM_RAISE-GENIC_210526.vcf.gz
Kuopio=Kuopio/AnVIL_CCDG_Broad_NP_Epilepsy_FINKPH_EPIL_CO_MORBIDI_MDS_WES_y2-only-raisegenic.vcf.gz


echo "calculate stats"
bcftools stats --threads 20 $Brussels > $Brussels.vchk
bcftools stats --threads 20 $Frankfurt > $Kuopio.vchk
bcftools stats --threads 20 $Kuopio > $Kuopio.vchk

echo "todo liftover to hg38"




echo "merging all"
# todo add Brussels
bcftools merge  $Frankfurt $Kuopio --threads 20 -0 -Oz -o RAISE_Genic_WES_Frankfurt_Kuopio_merged.vcf.gz
bcftools index --threads 20 RAISE_Genic_WES_Frankfurt_Kuopio_merged.vcf.gz
bcftools stats --threads 20 RAISE_Genic_WES_Frankfurt_Kuopio_merged.vcf.gz > RAISE_Genic_WES_Frankfurt_Kuopio_merged.vcf.gz.vchk

bcftools view --threads 20 -i  'MIN(FMT/DP>10) && MIN(FMT/GQ>15)' RAISE_Genic_WES_Frankfurt_Kuopio_merged.vcf.gz -Oz -o RAISE_Genic_WES_Frankfurt_Kuopio_merged_filtered.vcf.gz

wait 
echo "merged"

cd $homedir

echo "done"

# todo convert name to macht annotation script

#!/bin/bash
#SBATCH --job-name=merger
#SBATCH --partition=fuchs
#SBATCH --nodes=2
#SBATCH --ntasks=20
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1200   
#SBATCH --time=8:00:00
#SBATCH --no-requeue
#SBATCH --mail-type=ALL

modus="test"
reffasta="/scratch/fuchs/agmisc/chiocchetti/ReferenceGenomes/hg38.fa"
refdb="/scratch/fuchs/agmisc/chiocchetti/annovar/humandb/"

filename=$1

echo "the length of filename is ${#filename}"

if [[ ${#filename} == 0  ]] | [[ ! -f $filename ]];
then
    echo "file does not exist or is not specified";
    exit;
else
    echo "$filename is used"
    filedir=$(dirname $filename)
    filename=$(basename $filename)
fi

homedir=$(pwd)
cd $filedir

if [ "$modus" == "test" ]
then
    echo "modus is test: subset is used"
    if [ -f testset.vcf.gz ]
    then
	echo "using existing testset file; 
if you want to use anotherone delete or rename testset.vcf.gz in $filedir"
	filename="testset.vcf.gz"
    else
	echo "testset is created"
	bcftools query -l $filename | grep 21 > samples.tmp #selects all files with a 21 in it is kind of arbitrary
	bcftools view $filename -S samples.tmp -Oz -o testset.vcf.gz
	filename="testset.vcf.gz"
	rm samples.tmp
    fi
fi

bcftools index -f $filename
bcftools view --regions chr1 $filename --threads 20 -o tmp.vcf


echo "run normalization step1"
bcftools norm -m-both --threads 20 -o 01_merged.vcf.gz tmp.vcf

echo "run normalization step2"
bcftools norm tmp_01_merged.vcf -f $reffasta --threads 20  -Oz -o  tmp_02_merged.vcf.gz
bcftools index tmp_02_merged.vcf.gz
## todo hg 38
echo "run annotation"
## annotate_variation.pl -vcfinput -out Annotated -build hg19 02_merged.vcf.gz humandb/

table_annovar.pl tmp_02_merged.vcf.gz  $refdb \
		 -buildver hg38 \
		 -out Annotated \
		 -remove \
		 -protocol refGene,avsnp147,dbnsfp30a\
		 -operation g,f,f\
		 -nastring . \
		 -vcfinput \
		 -polish\
		 -thread 20\
		 -maxgenethread 20



htsfile -h  Annotated.hg38_multianno.vcf > tmp_LGD.vcf
grep -E "stop|start|frameshift|splicing" Annotated.hg38_multianno.vcf >> tmp_LGD.vcf

bgzip -c -@ 20 tmp_LGD.vcf > LGD_$filename.gz
bcftools index LGD_$filename.gz

rm tmp*
rm Annotated*

## aim list all frameshift stopgains splice site 
## bcftools query -f '%ExonicFunc.refGene ' Annotated.hg38_multianno.vcf  | sort | uniq

## ExonicFunc.refGene
# .
# frameshift_deletion
# frameshift_insertion
# nonframeshift_deletion
# nonframeshift_insertion
# nonsynonymous_SNV
# startloss
# stopgain
# stoploss
# synonymous_SNV
# unknown

## Func.refGene
# downstream
# exonic
# exonic\x3bsplicing
# intergenic
# intronic
# ncRNA_exonic
# ncRNA_exonic\x3bsplicing
# ncRNA_intronic
# ncRNA_splicing
# ncRNA_UTR5
# splicing
# upstream
# upstream\x3bdownstream
# UTR3
# UTR5
# UTR5\x3bUTR3


echo "done"

cd $homedir



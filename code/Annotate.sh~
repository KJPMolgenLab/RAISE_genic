#!/bin/bash
#SBATCH --job-name=merger
#SBATCH --partition=fuchs
#SBATCH --nodes=2
#SBATCH --ntasks=20
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1200   
#SBATCH --time=12:00:00
#SBATCH --no-requeue
#SBATCH --mail-type=ALL

modus="notest"

reffasta="/scratch/fuchs/agmisc/chiocchetti/ReferenceGenomes/hg38.fa"
refdb="/scratch/fuchs/agmisc/chiocchetti/annovar/humandb/"

filename=$1

echo "the length of filename is ${#filename}"

if [[ ${#filename} == 0  ]] | [[ ! -f $1 ]];
then
    echo "file does not exist or is not specified";
    exit;
else
    echo "$1 is used"
    filename=$(basename $1)
    filedir=$(dirname $1)
fi

homedir=$(pwd)
cd $filedir

if [ "$modus" == "test" ]
then
    echo "modus is test: subset is used"
    if [ -f testset.vcf.gz ]
    then
	echo "using existing testset file; if you want to use an otherone delete or rename testset.vcf.gz in $filedir"
	filename="testset.vcf.gz"
    else
	bcftools query -l $filename | grep 5 > samples.tmp
	bcftools view $filename -S samples.tmp -Oz -o testset.vcf.gz
	filename="testset.vcf.gz"
    fi
fi

echo "run normalization step1"
bcftools norm -m-both --threads 20 -Oz -o 01_merged.vcf.gz $filename

echo "run normalization step2"
bcftools norm -f $reffasta --threads 20 -Oz -o 02_merged.vcf.gz 01_merged.vcf.gz

## todo hg 38
echo "run annotation"
## annotate_variation.pl -vcfinput -out Annotated -build hg19 02_merged.vcf.gz humandb/

table_annovar.pl 02_merged.vcf.gz  $refdb \
		 -buildver hg38 \
		 -out Annotated \
		 -remove \
		 -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a\
		 -operation gx,r,f,f,f\
		 -nastring . \
		 -vcfinput \
		 -polish\
		 -thread 20\
		 -maxgenethread 20


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


bcftools query -e 'GT ="."' \
	 -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%Gene.refGene\t%GeneDetail.refGene\t%Func.refGene\t%ExonicFunc.refGene\t%AAChange.refGene\n' \
	 Annotated.hg38_multianno.vcf -threads 20 | \
    grep -E "stop|start|frameshift|splicing" > Annotated.hg38_LDG.txt

echo "done"

cd $homedir



#Read filter
for filename in *.fastq
do
# first, make the base by removing fastq.gz
        base=$(basename $filename .fastq)
        echo $base

       trimmomatic SE ${base}.fastq \
                ${base}.qc.fq \
                ILLUMINACLIP:/home/afsheenyousaf/anaconda3/pkgs/trimmomatic-0.39-1/share/trimmomatic-0.39-1/adapters/TruSeq3-SE.fa:2:0:15 \
                LEADING:3 TRAILING:3 \
                SLIDINGWINDOW:4:2 \
                MINLEN:36
done              
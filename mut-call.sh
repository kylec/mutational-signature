T=02
N=1
BAM="*.bam"
SNVDIR=~/data/sol-snv
COVDIR=~/data/sol-cov
COSMIC=~/references/hg18_cosmic_v54_120711.vcf
REF=~/references/hsap_36.1_hg18.fa
DBSNP=~/references/bcm_dbsnp_132.hg18.vcf

for PAT in `ls $BAM | cut -d- -f1-3 | sort -u`; do
TUM=`ls $BAM | grep $PAT-$T`; NRM=`ls $BAM | grep $PAT-$N`;
q "java -Xmx4g -jar ~/src/jar/muTect-1.1.4.jar --analysis_type MuTect --reference_sequence $REF --dbsnp $DBSNP --cosmic $COSMIC --input_file:normal $NRM --input_file:tumor $TUM --out $SNVDIR/$PAT-$T.out.txt --coverage_file $COVDIR/$PAT-$T.wig.txt" $PAT-$T $PAT-$T.log;
done
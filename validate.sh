# query all rna bams
cgquery "disease_abbr=OV&analyte_code=R&sample_type=01&library_strategy=RNA-Seq&state=live&filetype=bam" -o prim-output.xml

cgquery "disease_abbr=OV&analyte_code=R&sample_type=02&library_strategy=RNA-Seq&state=live&filetype=bam" -o recur-output.xml

# construct cgquery
xml=recur-output.xml
for a in `cat solid-download.txt`; do file=`grep $a $xml  | grep -v -P "bai|GRC" | sed 's/.*<filename>//' | sed 's/<\/filename>//'` ; if [ $file ]; then cgquery "filename=$file" -o $file.xml; fi; done

# download
for a in TCGA*02A*.xml; do 
q "/apps/x86_64/GT302/bin/GeneTorrent -C /apps/x86_64/GT302/share/GeneTorrent -v -c ~/references/cghub.key -d $a --maxChildren 1" $a $a.log
done

#download by date
for a in `find . -type f -newermt 2013-12-23 -name "TCGA*.xml"`; do q "/apps/x86_64/GT302/bin/GeneTorrent -C /apps/x86_64/GT302/share/GeneTorrent -v -c ~/references/cghub.key -d $a --maxChildren 1" $a $a.log; done

# keep snvs
for a in *.out.txt; do grep -P "judgement|KEEP" $a > $a.keep; done

# get cov
for a in *.keep; do 
q "perl ~/src/scripts/getTumCov.pl -i $a -bamdir ~/data/sol-bam -short > $a.cov" $a $a.log
done

# get rna cov
for a in *.cov; do q "perl ~/src/scripts/getTumCov.pl -i $a -bamdir ~/rna/ -short -label RNA " $a $a.log; done

# add ttot tvar seq cov at the end
for a in *.cov.cov; do head -1 $a | sed 's/$/\tTTotSeq\tTVarSeq/' > $a.fmt; sed '1d' $a | awk '{FS=OFS="\t"; print $0, $21+$22, $22}' >> $a.fmt; done

##### filter mutation call
for a in *.fmt; do perl ~/src/scripts/filter.pl $a 20 2 20 2 20 2 > $a.filt; done
for a in *.filt; do py2 ~/src/scripts/addRefBase.py $a; done

### count novel/dbsnp
for a in *.annotated; do sample=`echo $a | cut -d. -f1 | sed 's/TCGA-//g'`; novel=`cut -f42 $a | grep novel| wc -l`; total=`sed '1d' $a | wc -l`; echo -e "$sample\t$novel\t$total"; done

#### annotation count
echo -e "sample\texonic\trna\tutr\tintron\tintergenic"
for a in *.annotated; do sample=`echo $a | cut -d. -f1 | sed 's/TCGA-//g'`; code=`cut -f43 $a | grep -P "SNV|splicing" | wc -l`; rna=`cut -f43 $a | grep RNA | wc -l`; utr=`cut -f43 $a | grep UTR | wc -l`; intron=`cut -f43 $a | grep intronic | wc -l`; int=`cut -f43 $a | grep intergenic | wc -l`; echo -e "$sample\t$code\t$rna\t$utr\t$intron\t$int"; done

### recurrent gene
 cut -f43,44 *.annotated | awk -F"\t" '$1 ~ /SNV|splicing|exonic/'  | cut -f2 | sort | uniq -c | awk '$1>1'

#### clinical info (patient, age, days to death)
grep -f<(cut -d- -f1-2 ~/data/ill-snv/barcode.txt) nationwidechildrens.org_clinical_patient_ov.txt  | cut -f1,2,12

# create matlab file for mutation signature
paste *filt.subtypes > originalGenomes.txt
paste TCGA*01.out*filt.subtypes > prim-originalGenomes.txt
paste TCGA*02.out*filt.subtypes > recur-originalGenomes.txt
matlab -nodesktop -nosplash -nodisplay -r "run /home/kchang3/src/scripts/makeMutClassInput.m;quit;"
scp *.mat kchang3@mdarisngc03:~

##### run matlab ######
q "matlab -nodesktop -nodisplay -r \"test('all');quit;\"" all all.log
q "matlab -nodesktop -nodisplay -r \"test('prim');quit;\"" prim prim.log
q "matlab -nodesktop -nodisplay -r \"test('recur');quit;\"" recur recur.log

# cov columns
 cut -f21,22,36- *.cov

# solid 
head -1 ova-triplets-primary-somatic-wildtype.val.mafplus > ova-solid-maf.txt
awk -F"\t" '$10=="SNP" && $26=="Somatic"' *.mafplus  >> ova-solid-maf.txt
perl ~/src/scripts/addSeqContext.pl ova-solid-maf.txt > ova-solid-maf-seq.txt

# split file
IN=ova-solid-maf-seq.txt; for a in `cut -f16 ova-solid-maf-seq.txt | sort -u | cut -d- -f1-4 | sed 's/A$//g' | sed 's/B$//g' | grep -v Tumor `; do OUT=$a.txt; head -1 $IN > $OUT; grep $a $IN >> $OUT; done

# count mut
for a in TCGA*.txt; do py2 ~/src/scripts/solid_addRefBase.py $a; done

# clincial
grep -f<(cut -d- -f1-2 patients.txt) nationwidechildrens.org_clinical_patient_ov.txt  | cut -f1,2,9,12,50
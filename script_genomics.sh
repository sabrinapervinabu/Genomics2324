#!/bin/bash

if [ ! -d "project"]; then 
	mkdir project
cd project

#create a directory for each assigned case
cases=(590 593 601 631 669 672 681 708 710 745)
for num in "${cases[@]}"; do
	if [ ! -d "case$num" ]; then
		mkdir case$num
	fi
done

#check the FASTA files
names=("mother" "father" "child")
for num in "${cases[@]}"; do
	for name in "${names[@]}"; do
		fastqc /home/BCG2024_genomics_exam/case${num}_${name}.fq.gz -o case${num}/
	done
done

#create the BAM files
for num in "${cases[@]}"; do
	for name in "${names[@]}"; do
		bowtie2 -U /home/BCG2024_genomics_exam/case${num}_${name}.fq.gz \
			-p 8 -x /home/BCG2024_genomics_exam/uni --threads 4 \
			--rg-id "$name" --rg "SM:$name" \
			| samtools view -Sb -@ 4 \
			| samtools sort -@ 4 -o case${num}/case${num}_${name}.bam
	done
done

#create the BG files from each BAM file
for num in "${cases[@]}"; do
	for name in "${names[@]}"; do
		bedtools genomecov -ibam case${num}/case${num}_${name}.bam \
			-bg -trackline -trackopts "name="$name"" -max 100 \
			> case${num}/${name}Cov.bg
	done
done

#do the quality check for each BAM file (with target exome)
for num in "${cases[@]}"; do
	for name in "${names[@]}"; do
		qualimap bamqc -bam case${num}/case${num}_${name}.bam \
			-gff /home/BCG2024_genomics_exam/exons16Padded_sorted.bed \
			-outdir case${num}/case${num}_${name}
	done
done

#do multiqc
for num in "${cases[@]}"; do
	multiqc case${num}/ -o case${num}/
done

#create the multi-sample VCF files
mkdir results
mkdir final

for num in "${cases[@]}"; do
	freebayes -f /home/BCG2024_genomics_exam/universe.fasta \
		-m 20 -C 5 -Q 10 --min-coverage 10 case${num}/case${num}_mother.bam \
		case${num}/case${num}_father.bam case${num}/case${num}_child.bam \
		> results/case${num}.vcf
done

#sort the sample’s names in the VCF files
for num in "${cases[@]}"; do
	bcftools query -l results/case${num}.vcf \
		| sort > results/samples_${num}.txt
done

for num in "${cases[@]}"; do
	bcftools view -S results/samples_${num}.txt results/case${num}.vcf \
		> results/case${num}.sorted.vcf
done

rm results/*.txt

#variant prioritization
AR=(601 708 745)
AD=(590 593 631 669 672 681 710)

AD_pipeline () {
	num=$1
	grep "#" results/case${num}.sorted.vcf > results/candilist${num}.vcf
	grep "0/1.*0/0.*0/1" results/case${num}.sorted.vcf >> results/candilist${num}.vcf
	#grep "0/1.*0/0.*1/0" results/case${num}.sorted.vcf >> results/candilist${num}.vcf
	#grep "1/0.*0/0.*0/1" results/case${num}.sorted.vcf >> results/candilist${num}.vcf
	#grep "1/0.*0/0.*1/0" results/case${num}.sorted.vcf >> results/candilist${num}.vcf
	#grep "1/0.*0/1.*0/0" results/case${num}.sorted.vcf >> results/candilist${num}.vcf
	#grep "1/0.*1/0.*0/0" results/case${num}.sorted.vcf >> results/candilist${num}.vcf
	grep "0/1.*0/1.*0/0" results/case${num}.sorted.vcf >> results/candilist${num}.vcf
	#grep "0/1.*1/0.*0/0" results/case${num}.sorted.vcf >> results/candilist${num}.vcf
	#grep "1/0.*0/0.*0/0" results/case${num}.sorted.vcf >> results/candilist${num}.vcf
	grep "0/1.*0/0.*0/0" results/case${num}.sorted.vcf >> results/candilist${num}.vcf
	bedtools intersect -a results/candilist${num}.vcf \
		-b /home/BCG2024_genomics_exam/exons16Padded_sorted.bed -u > final/${num}TG.vcf
}

for num in "${AD[@]}"; do
	AD_pipeline "$num" &
done


AR_pipeline () {
	num=$1
	grep "#" results/case${num}.sorted.vcf > results/candilist${num}.vcf
	grep "1/1.*0/1.*0/1" results/case${num}.sorted.vcf >> results/candilist${num}.vcf
	#grep "1/1.*0/1.*1/0" results/case${num}.sorted.vcf >> results/candilist${num}.vcf
	#grep "1/1.*1/0.*0/1" results/case${num}.sorted.vcf >> results/candilist${num}.vcf
	#grep "1/1.*1/0.*1/0" results/case${num}.sorted.vcf >> results/candilist${num}.vcf
	bedtools intersect -a results/candilist${num}.vcf \
		-b /home/BCG2024_genomics_exam/exons16Padded_sorted.bed -u > final/${num}TG.vcf
}

for num in "${AR[@]}"; do
	AR_pipeline "$num" &
done


#END OF THE SCRIPT


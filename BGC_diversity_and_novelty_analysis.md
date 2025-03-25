## BGC analysis of CRBC and published genomes

### 1. The genomic data processing and AI-2-related proteins prospecting
```shell
# run checkm
time checkm lineage_wf -t 40 -x fna $source/  $wd/checkm/
```


```python
#!/bin/python3

'''
This program is designed to process CheckM output file (e.g. bin_stats_ext.tsv) to easy-to-read text format (.txt)
Usage: python3 checkm_summary.py <inputfile> <outputfile>
Written by: Emmett Peng
'''
import json
import sys
input_="D:/bin_stats_ext.tsv"
out_="D:/checkout.txt"

with open(input_ ,'r') as f:
	Load = {}
	for line in f:
		line = line.replace('\'','\"') #替换
		line = line.split('\t')
		line[0] = 'bin_' + line[0].replace('.bin.', '_')		
		Load[line[0]] = json.loads(line[1])		

with open(out_, 'w+') as output:
	output.write('Bin Id\tMarker Lineage\tGenomes\tMarkers\tMarker Sets\t0\t1\t2\t3\t4\t5+\tCompleteness\tContamination\tGC\tGC std\tGenome Size\tAmbiguous Bases\tScaffolds\tContigs\tTranslation Table\tPredicted Genes\n')
	for key in Load:
		output.write(key + '\t')
		output.write(Load[key]['marker lineage'] + '\t')
		output.write(str(Load[key]['# genomes']) + '\t')
		output.write(str(Load[key]['# markers']) + '\t')
		output.write(str(Load[key]['# marker sets']) + '\t')
		output.write(str(Load[key]['0']) + '\t')
		output.write(str(Load[key]['1']) + '\t')
		output.write(str(Load[key]['2']) + '\t')
		output.write(str(Load[key]['3']) + '\t')
		output.write(str(Load[key]['4']) + '\t')
		output.write(str(Load[key]['5+']) + '\t')
		output.write(str(Load[key]['Completeness']) + '\t')
		output.write(str(Load[key]['Contamination']) + '\t')
		output.write(str(Load[key]['GC']) + '\t')
		output.write(str(Load[key]['GC std']) + '\t')
		output.write(str(Load[key]['Genome size']) + '\t')
		output.write(str(Load[key]['# ambiguous bases']) + '\t')
		output.write(str(Load[key]['# scaffolds']) + '\t')
		output.write(str(Load[key]['# contigs']) + '\t')
		output.write(str(Load[key]['Translation table']) + '\t')
		output.write(str(Load[key]['# predicted genes']) + '\t')
		output.write('\n')
print("OK")
```
```shell
#run prokka
ls *.fa | parallel --verbose "prokka --cpu 40 {} --prefix {.}" &
#
cat prokka/*.faa > prokka.faa
#run hummsearch
hmmsearch -o PF02664.txt --cut_ga PF02664.hmm prokka.faa --domE 1e-5
hmmsearch -o PF13407.txt --cut_ga PF13407.hmm prokka.faa --domE 1e-5
hmmsearch -o PF17155.txt --cut_ga PF17155.hmm prokka.faa --domE 1e-5
hmmsearch -o PF02743.txt --cut_ga PF02743.hmm prokka.faa --domE 1e-5
```
### 2. The construction of phylogenetic trees based on microbial genomes 
```shell
#run gtdb
gtdbtk classify_wf --genome_dir seq/ --out_dir gtdb_classify --extension fa --prefix tax --cpus 50 --skip_ani_screen 
#run fasttree
FastTree protein_alignment > tree
```
### 3. The signal transduction modes architecture of AI-2 
```shell
#run pfam
pfam_scan.pl -fasta pfam/domainseq.txt -dir /pfam/ -outfile results_0.5.faa -as 
```


### 4. The gene expression within mouse gut microbial metatranscriptomes
```shell
#run  Trimmomatic 
for i in `cat name.txt`;do
    trimmomatic PE -phred33 -threads 40 $i_1.fastq.gz $i_2.fastq.gz $i_1.clean.fq.gz $i_1.unpaired.fq.gz $i_2.clean.fq.gz $i_2.unpaired.fq.gz ILLUMINACLIP:./adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:5:20 LEADING:5 TRAILING:5 MINLEN:30 
done
#run fastqc
time fastqc qc/*clean.fq.gz -t 40
multiqc -d qc/ -o mqc/
#run STAR v2.7.10a50
STAR --runThreadN 40 --runMode genomeGenerate \
--genomeDir arab_STAR_genome \
--genomeFastaFiles mou933.fas \
--sjdbGTFfile gtfmouse933new.gtf \
--limitGenomeGenerateRAM 92048356277 \
--sjdbOverhang 149 

STAR --runThreadN 20 --genomeDir arab_STAR_genome \
--readFilesIn ~/meta/clean/mouse933clean/clean/ERR1146183_1.clean.fq ~/meta/clean/mouse933clean/clean/ERR1146183_2.clean.fq  \
--outFileNamePrefix 03align_out/mou \
--outSAMtype BAM SortedByCoordinate \
--outBAMsortingThreadN 10 \
--quantMode TranscriptomeSAM GeneCounts 
#run RSEM v1.351
rsem-prepare-reference --gtf gtfmouse933new.gtf \
mou933.fas \
arab_RSEM/arab_rsem

rsem-calculate-expression --paired-end --no-bam-output --alignments -p 35 -q 03align_out/mouse33Aligned.toTranscriptome.out.bam arab_RSEM/arab_rsem 04rsem_out/mouse33
```



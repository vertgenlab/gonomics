# Generating Testing Files for *samAssembler*

  ### Riley J. Mangan

In this file, I'll explain how I generated the test files used to evaluate *samAssembler*.

 First, we'll need a reference genome sequence. In this case, I generated two 10kb random sequences with *randSeq*, corresponding to a 20kb genome with 2 chromosomes, each 10kb in length.

  ```
~/go/bin/randSeq -numSeq 2 -lenSeq 10000 -setSeed 23 ref.fa
  ```

Next, I generated two evolved sequences with branch length 0.05 with simulateEvol for each chromosome of the reference.

  

```
~/go/bin/faFilter -name Sequence_0 ref.fa Sequence_0.fa
~/go/bin/simulateEvol -qName evol1 -branchLength 0.05 -propIndel 0.1 -setSeed 19 -transitionBias 3 Sequence_0.fa Sequence_0.evol1.fa
~/go/bin/simulateEvol -qName evol2 -branchLength 0.05 -propIndel 0.1 -setSeed 17 -transitionBias 3 Sequence_0.fa Sequence_0.evol2.fa

~/go/bin/faFilter -name Sequence_1 ref.fa Sequence_1.fa
~/go/bin/simulateEvol -qName evol1 -branchLength 0.05 -propIndel 0.1 -setSeed 27 -transitionBias 3 Sequence_1.fa Sequence_1.evol1.fa
~/go/bin/simulateEvol -qName evol2 -branchLength 0.05 -propIndel 0.1 -setSeed 29 -transitionBias 3 Sequence_1.fa Sequence_1.evol2.fa
```
  
I then filtered out the evolved sequences from the resulting multiFa alignments and trimmmed out gaps:

 ``` 
~/go/bin/faFilter -notName Sequence_0 Sequence_0.evol1.fa tmp.Sequence_0.evol1.fa
~/go/bin/faFormat -noGaps tmp.Sequence_0.evol1.fa Sequence_0.evol1.diverge.noGaps.fa
~/go/bin/faFilter -notName Sequence_1 Sequence_1.evol1.fa tmp.Sequence_1.evol1.fa
~/go/bin/faFormat -noGaps tmp.Sequence_1.evol1.fa Sequence_1.evol1.diverge.noGaps.fa
```
  
```
~/go/bin/faFilter -notName Sequence_0 Sequence_0.evol2.fa tmp.Sequence_0.evol2.fa
~/go/bin/faFormat -noGaps tmp.Sequence_0.evol2.fa Sequence_0.evol2.diverge.noGaps.fa
~/go/bin/faFilter -notName Sequence_1 Sequence_1.evol2.fa tmp.Sequence_1.evol2.fa
~/go/bin/faFormat -noGaps tmp.Sequence_1.evol2.fa Sequence_1.evol2.diverge.noGaps.fa
```
 ```
cat Sequence_0.evol1.diverge.noGaps.fa Sequence_1.evol1.diverge.noGaps.fa > evol1.diverge.noGaps.fa
cat Sequence_0.evol2.diverge.noGaps.fa Sequence_1.evol2.diverge.noGaps.fa > evol2.diverge.noGaps.fa
```

We'll merge the diverged sequences and place them in the multiFa folder for validation later
```bash
~/go/bin/mergeMultiFa Sequence_0.evol1.fa Sequence_0.evol2.fa multiFa/Sequence_0.evol.fa
~/go/bin/mergeMultiFa Sequence_1.evol1.fa Sequence_1.evol2.fa multiFa/Sequence_1.evol.fa
```

Next I generated 60x coverage short reads.
(# of reads per chromosome = 2000 = 60xCoverage * 10kb Chrom Length / 2x150bp paired read length) 

```
~/go/bin/simulateSam -n 2000 -setSeed 9 -readLength 150 -fragmentLength 400 evol1.diverge.noGaps.fa evol1.Reads.bam
~/go/bin/simulateSam -n 2000 -setSeed 9 -readLength 150 -fragmentLength 400 evol2.diverge.noGaps.fa evol2.Reads.bam
rm evol1.diverge.noGaps.fa evol2.diverge.noGaps.fa
```
  
Next I recovered fastq reads from the simulated bam files and merged the fastq files together

```
~/Software/bedtools2/bin/bedtools bamtofastq -i evol1.Reads.bam -fq evol1.readA.fastq -fq2 evol1.readB.fastq
~/Software/bedtools2/bin/bedtools bamtofastq -i evol2.Reads.bam -fq evol2.readA.fastq -fq2 evol2.readB.fastq
cat evol1.readA.fastq evol2.readA.fastq > merged.readA.fastq
cat evol1.readB.fastq evol2.readB.fastq > merged.readB.fastq
```
  
After this, I aligned these reads to the reference genome and sorted the reads

 ```bash
~/Software/bwa-0.7.17/bwa index ref.fa
~/Software/bwa-0.7.17/bwa mem -t 8 ref.fa merged.readA.fastq merged.readB.fastq | ~/Software/samtools-1.9/samtools view -bh - > diverged.RefAln.bam
~/Software/samtools-1.9/samtools sort diverged.RefAln.bam -o diverged.RefAln.sorted.bam
```

We now have diploid testing data generated from a synthetic sequence to be used with samAssembler
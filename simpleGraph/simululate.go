package simpleGraph

import (
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
	//"os"
	"fmt"
	"sort"
	"strings"
)

func checkAlignment(aln *sam.SamAln) bool {
	var answer bool = false
	words := strings.Split(aln.QName, "_")
	var ref int64
	//var query int64
	var blastScore int64
	alignedPos := common.StringToInt64(words[1])
	if alignedPos <= ref-10 || (alignedPos > ref+100) {
		words = strings.Split(aln.Extra, "\t")
		blastScore = common.StringToInt64(words[0][5:])
		if blastScore > 5000 {
			answer = true
		}
		log.Printf("\t%s\t%s\t%d\t%s\t%s\n", cigar.ToString(aln.Cigar), aln.QName, aln.Pos, aln.RName, aln.Extra)
	}
	return answer
}

func CheckAnswers(query []*sam.SamAln) {
	var yes, no int64 = 0, 0
	for i := 0; i < len(query); i++ {
		if checkAlignment(query[i]) {
			yes++
			//log.Printf(sam.SamAlnToString(query[i]))
		} else {
			no++
			//log.Printf("This did not map:\n%s\n", sam.SamAlnToString(query[i]))
		}
	}
	log.Printf("Total number of reads aligned: %d...", len(query))
	log.Printf("Number of reads correctly aligned: %d...\n", yes)
	log.Printf("Number of reads mismapped: %d...\n", no)
}

type Location struct {
	Chr string
	Pos int64
}

//will only look at lines that include chr name
func samToGenomeNotes(samFileName string, chr *fasta.Fasta, chromIdx int) map[uint64][]string {
	var i int
	fasta.ToUpper(chr)
	samFile := fileio.EasyOpen(samFileName)
	defer samFile.Close()
	var done bool = false
	var RefIndex, SeqIndex, subSeqIdx int64
	var currentSeq []dna.Base
	var aln *sam.SamAln
	sam.ReadHeader(samFile)

	votes := make(map[uint64][]string)
	var seqCode uint64

	//var basesToAdd []dna.Base
	var progress int
	log.Printf("Reading in sam alignments...")
	//var currLoc int64
	//var currLoc Location = Location{Chr: "", Pos: 0}
	for aln, done = sam.NextAlignment(samFile); done != true; aln, done = sam.NextAlignment(samFile) {
		progress++
		SeqIndex = 0
		RefIndex = aln.Pos - 1
		if aln.Cigar[0].Op != '*' {
			for i = 0; i < len(aln.Cigar); i++ {
				currentSeq = aln.Seq
				if aln.Cigar[i].Op == 'D' {
					gapSeq := make([]dna.Base, aln.Cigar[i].RunLength)
					for subSeqIdx = 0; subSeqIdx < aln.Cigar[i].RunLength; subSeqIdx++ {
						gapSeq[subSeqIdx] = dna.Gap
					}
					seqCode = chromAndPosToNumber(chromIdx, int(RefIndex))
					votes[seqCode] = append(votes[seqCode], dna.BasesToString(gapSeq))
					//RefIndex, SeqIndex = cigar.UpdateIndices(aln.Cigar[i], RefIndex, SeqIndex)
				} else if aln.Cigar[i].Op == 'M' {
					for subSeqIdx = 0; subSeqIdx < aln.Cigar[i].RunLength; subSeqIdx++ {
						seqCode = chromAndPosToNumber(chromIdx, int(RefIndex+subSeqIdx))
						votes[seqCode] = append(votes[seqCode], dna.BaseToString(currentSeq[SeqIndex+subSeqIdx]))
					}
					//RefIndex++
					//SeqIndex++
					//RefIndex, SeqIndex = cigar.UpdateIndices(aln.Cigar[i], RefIndex, SeqIndex)
				} else if aln.Cigar[i].Op == 'I' {
					if i > 0 {
						seqCode = chromAndPosToNumber(chromIdx, int(RefIndex))
						votes[seqCode] = append(votes[seqCode], dna.BasesToString(currentSeq[SeqIndex-1:SeqIndex-1+aln.Cigar[i].RunLength]))
					
					}
				}else {
					//log.Printf("Skipped parts of this cigar: %s\n", cigar.ToString(aln.Cigar))
				}
				RefIndex, SeqIndex = cigar.UpdateIndices(aln.Cigar[i], RefIndex, SeqIndex)
			}
		}
	}
	log.Printf("Finished analyzing %d alignments...", progress)
	return votes
}
func toVcf(samFileName string, chr *fasta.Fasta, chromIdx int, votes map[uint64][]string) []*vcf.Vcf {
	var concensus []dna.Base
	var seqCode uint64
	var vcfs []*vcf.Vcf
	var base int64
	for base = 0; base < int64(len(chr.Seq)); base++ {
		seqCode = chromAndPosToNumber(chromIdx, int(base))
		_, ok := votes[seqCode]
		if ok {

			concensus = dna.StringToBases(MostOccuringSeq(votes[seqCode]))
			if len(concensus) >= 1 {
				if len(concensus) == 1 && concensus[0] != chr.Seq[seqCode] {
					snp := vcf.Vcf{Chr: chr.Name, Pos: base + 1, Id: "", Ref: dna.BaseToString(chr.Seq[base]), Alt: dna.BaseToString(concensus[0]), Qual: 0, Filter: "", Info: "", Notes: "SVTYPE=SNP"}
					vcfs = append(vcfs, &snp)
				}
				if len(concensus) > 1 {
					numGaps := dna.CountBaseInterval(concensus, dna.Gap, 0, len(concensus))
					if numGaps > 0 {
						delSeq := vcf.Vcf{Chr: chr.Name, Pos: base + 1, Id: "", Ref: dna.BasesToString(chr.Seq[base : base+int64(len(concensus))]), Alt: dna.BaseToString(chr.Seq[base]), Qual: 0, Filter: "", Info: "", Notes: "SVTYPE=DEL"}
						vcfs = append(vcfs, &delSeq)
					} else {
						ins := vcf.Vcf{Chr: chr.Name, Pos: base + 1, Id: "", Ref: dna.BaseToString(chr.Seq[base]), Alt: dna.BasesToString(concensus), Qual: 0, Filter: "", Info: "", Notes: "SVTYPE=INS"}
						vcfs = append(vcfs, &ins)
					}
				}
			}
		}
	}
	return vcfs
}

func EditGenome(samFileName string, chr *fasta.Fasta, chromIdx int, votes map[uint64][]string) (*fasta.Fasta, []*vcf.Vcf) {
	var editGenome *fasta.Fasta = chr
	var concensus []dna.Base
	var seqCode uint64
	var vcfs []*vcf.Vcf
	var base int64
	for base = int64(len(chr.Seq) - 1); 0 <= base; base-- {
		seqCode = chromAndPosToNumber(chromIdx, int(base))
		_, ok := votes[seqCode]
		if ok {
			concensus = dna.StringToBases(MostOccuringSeq(votes[seqCode]))
			if len(concensus) >= 1 {
				if len(concensus) == 1 && concensus[0] != chr.Seq[seqCode] {
					snp := vcf.Vcf{Chr: chr.Name, Pos: base + 1, Id: "", Ref: dna.BaseToString(chr.Seq[base]), Alt: dna.BaseToString(concensus[0]), Qual: 0, Filter: "", Info: "", Notes: "SVTYPE=SNP"}
					vcfs = append(vcfs, &snp)
					editGenome.Seq[base] = concensus[0]
				}
				if len(concensus) > 1 {
					numGaps := dna.CountBaseInterval(concensus, dna.Gap, 0, len(concensus))
					if numGaps > 0 {
						delSeq := vcf.Vcf{Chr: chr.Name, Pos: base + 1, Id: "", Ref: dna.BasesToString(chr.Seq[base : base+int64(len(concensus))]), Alt: dna.BaseToString(chr.Seq[base]), Qual: 0, Filter: "", Info: "", Notes: "SVTYPE=DEL"}
						vcfs = append(vcfs, &delSeq)
						editGenome.Seq = dna.Delete(editGenome.Seq, base, base + int64(len(concensus)))
					} else  {
						ins := vcf.Vcf{Chr: chr.Name, Pos: base + 1, Id: "", Ref: dna.BaseToString(chr.Seq[base]), Alt: dna.BasesToString(concensus), Qual: 0, Filter: "", Info: "", Notes: "SVTYPE=INS"}
						vcfs = append(vcfs, &ins)
						editGenome.Seq = dna.Insert(editGenome.Seq, base+1, concensus[1:])
					}
				}
			}
		}
	}
	vcf.Sort(vcfs)
	return editGenome, vcfs
}

func SimulateVcfGenomeWide(genome []*fasta.Fasta, numChanges int) []*vcf.Vcf {

	var vcfs []*vcf.Vcf = make([]*vcf.Vcf, 0)
	for i := 0; i < len(genome); i++ {
		chrVcf := make([]*vcf.Vcf, numChanges)
		size := randIntInRange(2, 8)
		locations := make([]int, numChanges)
		var j int
		for j = 0; j < numChanges; {

			loc := randIntInRange(0, len(genome[i].Seq)-size)
			if genome[i].Seq[loc] != dna.N {
				locations[j] = loc
				j++
			}
		}
		sort.Sort(sort.Reverse(sort.IntSlice(locations)))
		for j = 0; j < numChanges; j++ {
			currVcf := vcf.Vcf{}
			currVcf = SimulateVcf(genome[i], int64(locations[j]), size, i)
			chrVcf[j] = &currVcf
		}
		vcfs = append(vcfs, chrVcf...)
	}
	return vcfs
}

func SimulateVcf(chr *fasta.Fasta, location int64, size int, chrom int) vcf.Vcf {
	editGenome := randIntInRange(0, 3)
	possibleBases := []dna.Base{dna.A, dna.C, dna.G, dna.T}
	extra := make([]dna.Base, size)
	if editGenome == 0 {
		var before dna.Base = chr.Seq[location]
		mutatePos(chr.Seq, int(location))
		snp := vcf.Vcf{Chr: chr.Name, Pos: location + 1, Id: fmt.Sprint(chrom), Ref: dna.BaseToString(before), Alt: dna.BaseToString(chr.Seq[location]), Qual: 0, Filter: "", Info: fmt.Sprintf("%s:SNP:%sto%s", chr.Name, dna.BaseToString(before), dna.BaseToString(chr.Seq[location])), Format: "SVTYPE=SNP"}
		return snp
	} else if editGenome == 1 {
		for base := 0; base < size; base++ {
			extra[base] = possibleBases[randIntInRange(0, len(possibleBases))]
		}
		dna.Insert(chr.Seq, location, extra)
		ins := vcf.Vcf{Chr: chr.Name, Pos: location + 1, Id: fmt.Sprint(chrom), Ref: dna.BaseToString(chr.Seq[location-1]), Alt: dna.BaseToString(chr.Seq[location-1]) + dna.BasesToString(extra), Qual: 0, Filter: "", Info: fmt.Sprintf("%s:insertion:%dto%d", chr.Name, location, location+int64(size)), Format: "SVTYPE=INS"}
		return ins
	} else {
		del := vcf.Vcf{Chr: chr.Name, Pos: location + 1, Id: fmt.Sprint(chrom), Ref: dna.BasesToString(chr.Seq[location-1 : location+int64(size)-1]), Alt: dna.BaseToString(chr.Seq[location-1]), Qual: 0, Filter: "", Info: fmt.Sprintf("%s:deletion:%dto%d", chr.Name, location, location+int64(size)), Format: "SVTYPE=DEL"}
		dna.Delete(chr.Seq, location, location+int64(size))
		return del
	}
}

func SimulateVcfFastq(genome []*fasta.Fasta, numChanges int, readLength int) ([]*vcf.Vcf, []*fastq.Fastq) {

	var vcfs []*vcf.Vcf = make([]*vcf.Vcf, 0)
	var simReads []*fastq.Fastq = make([]*fastq.Fastq, 0)
	for i := 0; i < len(genome); i++ {
		chrVcf := make([]*vcf.Vcf, numChanges)
		fqs := make([]*fastq.Fastq, numChanges)
		size := randIntInRange(2, 8)
		locations := make([]int, numChanges)
		var j int
		for j = 0; j < numChanges; {

			loc := randIntInRange(0, len(genome[i].Seq)-size)
			if genome[i].Seq[loc] != dna.N {
				locations[j] = loc
				j++
			}
		}
		sort.Sort(sort.Reverse(sort.IntSlice(locations)))
		for j = 0; j < numChanges; j++ {
			currVcf := vcf.Vcf{}
			currFastq := fastq.Fastq{}
			currVcf, currFastq = vcfFastq(genome[i], int64(locations[j]), size, i, readLength)
			chrVcf[j], fqs[i] = &currVcf, &currFastq
		}
		vcfs = append(vcfs, chrVcf...)
		simReads = append(simReads, fqs...)
	}
	return vcfs, simReads
}

func vcfFastq(chr *fasta.Fasta, location int64, size int, chrom int, readLength int) (vcf.Vcf, fastq.Fastq) {
	editGenome := randIntInRange(0, 3)
	possibleBases := []dna.Base{dna.A, dna.C, dna.G, dna.T}
	extra := make([]dna.Base, size)
	var half int64 = int64(readLength / 2)
	if editGenome == 0 {
		var before dna.Base = chr.Seq[location]
		mutatePos(chr.Seq, int(location))
		snp := vcf.Vcf{Chr: chr.Name, Pos: location + 1, Id: fmt.Sprint(chrom), Ref: dna.BaseToString(before), Alt: dna.BaseToString(chr.Seq[location]), Qual: 0, Filter: "", Info: fmt.Sprintf("%s:SNP:%sto%s", chr.Name, dna.BaseToString(before), dna.BaseToString(chr.Seq[location])), Format: "SVTYPE=SNP"}
		read := fastq.Fastq{Name: snp.Info, Seq: make([]dna.Base, readLength)}
		copy(read.Seq, chr.Seq[common.MaxInt64(location-half, 0):common.MaxInt64(location-half, 0)+int64(readLength)])
		return snp, read
	} else if editGenome == 1 {
		for base := 0; base < size; base++ {
			extra[base] = possibleBases[randIntInRange(0, len(possibleBases))]
		}
		dna.Insert(chr.Seq, location, extra)
		ins := vcf.Vcf{Chr: chr.Name, Pos: location + 1, Id: fmt.Sprint(chrom), Ref: dna.BaseToString(chr.Seq[location-1]), Alt: dna.BaseToString(chr.Seq[location-1]) + dna.BasesToString(extra), Qual: 0, Filter: "", Info: fmt.Sprintf("%s:insertion:%dto%d", chr.Name, location, location+int64(size)), Format: "SVTYPE=INS"}
		read := fastq.Fastq{Name: ins.Info, Seq: make([]dna.Base, readLength)}
		copy(read.Seq, chr.Seq[common.MaxInt64(location-half, 0):common.MaxInt64(location-half, 0)+int64(readLength)])

		return ins, read
	} else {
		del := vcf.Vcf{Chr: chr.Name, Pos: location + 1, Id: fmt.Sprint(chrom), Ref: dna.BasesToString(chr.Seq[location-1 : location+int64(size)-1]), Alt: dna.BaseToString(chr.Seq[location-1]), Qual: 0, Filter: "", Info: fmt.Sprintf("%s:deletion:%dto%d", chr.Name, location, location+int64(size)), Format: "SVTYPE=DEL"}
		dna.Delete(chr.Seq, location, location+int64(size))
		read := fastq.Fastq{Name: del.Info, Seq: make([]dna.Base, readLength)}
		copy(read.Seq, chr.Seq[common.MaxInt64(location-half, 0):common.MaxInt64(location-half, 0)+int64(readLength)])
		return del, read
	}
}

/*
//Generate fake vcf still work in progress
func makeFakeVcf(genome []*Node, readLength int, numReads int, numChanges int) ([]*fastq.Fastq, []*vcf.Vcf) {
	var vcfs []*vcf.Vcf = make([]*vcf.Vcf, numChanges)

	for i := 0; i < len(genome); i++ {
		locations := make([]int, numChanges)
		var j int
		for j = 0; j < numChanges;{
			loc := randIntInRange(0, len(genome[i].Seq))
			if genome[i].Seq[loc] != dna.N {
				locations[j] = loc
				j++
			}
		}
		sort.Sort(sort.Reverse(sort.IntSlice(locations)))
		for j = 0; j < numChanges; j++ {
			currVcf := vcf.Vcf{}
			genome[i], currVcf = SimulateMutations(genome[i], int64(locations[j]))
			vcfs[j]= &currVcf
		}
	}
	simReads := RandomReads(genome, readLength, numReads)
	return simReads, vcfs
}*/
/*
func EditGenome(chr *fasta.Fasta, votes map[Location][]string) []*vcf.Vcf {
	var concensus []dna.Base
	var currLoc Location = Location{Chr: chr.Name, Pos: 0}

	file, _ := os.Create("/dev/stdout")
	defer file.Close()
	var base int64
	for base = int64(len(chr.Seq) - 1); 0 <= base; base-- {
		currLoc.Pos = base
		concensus = dna.StringToBases(MostOccuringSeq(votes[currLoc]))
		variant := vcf.Vcf{Chr: chr.Name, Pos: 0, Id: "", Ref: "", Alt: "", Qual: 0, Filter: "", Info: ""}
		if len(concensus) != 0 {
			if chr.Seq[base] != concensus[0] {
				if len(concensus) == 1 {
					if concensus[0] == dna.Gap {


						variant.Notes = "SVTYPE=DEL"
					} else {
						variant.Notes = "SVTYPE=SNP"
					}
					chr.Seq[base] = concensus[0]
				} else {
					ins(chr.Seq, concensus, int(base))
					variant.Notes = "SVTYPE=INS"
				}
				fmt.Fprintf(file, "%s\t%v\t%s\t%s\t%s\t%v\t%s\t%s\t%s\n", variant.Chr, variant.Pos, variant.Id, variant.Ref, variant.Alt, variant.Qual, variant.Filter, variant.Info, variant.Notes)
			}
		}
	}
	return chr
}*/
/*
//will only look at lines that include chr name
func TakeNotesGenome(samFileName string, chr *fasta.Fasta, chromIdx int) map[uint64][]string {
	var i, k int
	samFile := fileio.EasyOpen(samFileName)
	defer samFile.Close()
	var done bool = false
	var RefIndex, SeqIndex int64
	var currentSeq []dna.Base
	var aln *sam.SamAln
	sam.ReadHeader(samFile)
	votes := make(map[uint64][]string)
	var seqCode uint64

	//var basesToAdd []dna.Base
	var progress int
	log.Printf("Reading in sam alignments...")
	for aln, done = sam.NextAlignment(samFile); done != true; aln, done = sam.NextAlignment(samFile) {
		progress++
		if aln.Cigar[0].Op != '*' {
			SeqIndex = 0
			RefIndex = aln.Pos - 1
			for i = 0; i < len(aln.Cigar); i++ {
				currentSeq = aln.Seq
				if aln.Cigar[i].Op == 'D' {
					seqCode = chromAndPosToNumber(chromIdx, int(RefIndex))
					votes[seqCode] = append(votes[seqCode], dna.BasesToString(chr.Seq[RefIndex:RefIndex+aln.Cigar[i].RunLength]))
					cigar.UpdateIndices(RefIndex, SeqIndex, aln.Cigar[i])
				} else if aln.Cigar[i].Op == 'M' {
					for k = 0; k < int(aln.Cigar[i].RunLength); k++ {
						seqCode = chromAndPosToNumber(chromIdx, int(RefIndex))
						votes[seqCode] = append(votes[seqCode], dna.BasesToString([]dna.Base{currentSeq[SeqIndex]}))
						RefIndex++
						SeqIndex++
					}
				} else if aln.Cigar[i].Op == 'I' && (SeqIndex-1 > 0) {
					seqCode = chromAndPosToNumber(chromIdx, int(RefIndex))
					votes[seqCode] = append(votes[seqCode], dna.BasesToString(currentSeq[SeqIndex:SeqIndex+aln.Cigar[i].RunLength]))
					cigar.UpdateIndices(RefIndex, SeqIndex, aln.Cigar[i])
				} else {
					cigar.UpdateIndices(RefIndex, SeqIndex, aln.Cigar[i])
				}
			}
		}
	}
	log.Printf("Finished analyzing %d alignments...", progress)
	return votes
}

func EditGenome(chr *fasta.Fasta, votes map[uint64][]string, chromIdx int) *fasta.Fasta {
	var concensus []dna.Base
	//var seqCode uint64
	var base int
	file, _ := os.Create("/dev/stdout")
	defer file.Close()
	variant := vcf.Vcf{Chr: chr.Name, Pos: 0, Id: "", Ref: "", Alt: "", Qual: 0, Filter: "", Info: ""}
	for base = len(chr.Seq) ; 0 <= base; base-- {
		seqCode := chromAndPosToNumber(chromIdx, base)
		concensus = dna.StringToBases(MostOccuringSeq(votes[seqCode]))
		if len(concensus) != 0 {
			if chr.Seq[base] != concensus[0] {
				if len(concensus) == 1 {
					if concensus[0] == dna.Gap {
						variant.Notes = "SVTYPE=DEL"
					} else {
						variant.Notes = "SVTYPE=SNP"
					}
					chr.Seq[base] = concensus[0]
				} else {
					ins(chr.Seq, concensus, base)
					variant.Notes = "SVTYPE=INS"
				}
				fmt.Fprintf(file, "%s\t%v\t%s\t%s\t%s\t%v\t%s\t%s\t%s\n", variant.Chr, variant.Pos, variant.Id, variant.Ref, variant.Alt, variant.Qual, variant.Filter, variant.Info, variant.Notes)
			}


		}

	}
	return chr
}*/

func ins(ref []dna.Base, ins []dna.Base, index int) {
	if index == len(ref) {
		ref = append(ref, ins...)
	} else if len(ins) == 1 {
		if ref[index] != ins[0] {
			ref[index] = ins[0]
		}
	} else {
		var tail []dna.Base = ref[index:]
		ref = append(ref[:index-1], ins...)
		ref = append(ref, tail...)
	}
}

func MostOccuringSeq(alleles []string) string {
	votes := make(map[string]int64)
	var mostCount int64 = -1
	var refSeq string
	for _, allele := range alleles {
		_, ok := votes[allele]
		if !ok {
			votes[allele] = 1
		} else {
			votes[allele]++
			if votes[allele] > mostCount {
				mostCount = votes[allele]
				refSeq = allele
			}
		}
	}
	return refSeq
}

/*
func WriteReadsToFile(genome []*Node, readLength int, numReads int, output string) {

	var numberOfReads int = 10000
	var readLength int = 150
	var mutations int = 1
	log.Printf("Simulating %d reads...\n", numberOfReads)
	simReads := RandomReads(genome, readLength, numberOfReads, mutations, true)
	fastq.Write("genomeDiversity.fq", simReads)
	log.Printf("Finished writing to file %d reads...\n", numberOfReads)
}*/

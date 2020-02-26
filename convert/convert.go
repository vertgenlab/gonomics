package convert

import (
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	"github.com/vertgenlab/gonomics/wig"
	"log"
	"strings"
)

func singleBedToFasta(b *bed.Bed, ref []*fasta.Fasta) *fasta.Fasta {
	for i := 0; i < len(ref); i++ {
		if b.Chrom == ref[i].Name {
			return fasta.Extract(ref[i], b.ChromStart, b.ChromEnd, b.Name)
		}
	}
	log.Fatalf("Chrom not found in fasta")
	return nil
}

func BedToFasta(b []*bed.Bed, ref []*fasta.Fasta) []*fasta.Fasta {
	outlist := make([]*fasta.Fasta, len(b))
	for i := 0; i < len(b); i++ {
		outlist[i] = singleBedToFasta(b[i], ref)
	}
	return outlist
}

func SamToBed(s *sam.SamAln) *bed.Bed {
	if s.Cigar[0].Op == '*' {
		return nil
	} else {
		return &bed.Bed{Chrom: s.RName, ChromStart: s.Pos - 1, ChromEnd: s.Pos -1 + cigar.ReferenceLength(s.Cigar), Name: s.QName}
	}
}

/* TODO: Write Sam to Bed conversion for paired reads.

func SamToBedPaired(s *sam.Sam) []*bed.Bed {
	//sort sam by QName
	//check for "properly aligned" flag
	//grab two properly paired samAln (sanme QName with strings.suffix removed), feed into helper function for bed conversion
	//add output to bedlist
} */

func SamToBedFrag(s *sam.SamAln, fragLength int64, reference map[string]*chromInfo.ChromInfo) *bed.Bed {
	var answer *bed.Bed

	if s.Cigar[0].Op == '*' {
		return nil
	} else {
		answer = &bed.Bed{Chrom: s.RName, Name: s.QName}
		if sam.IsPosStrand(s) {
			answer.ChromStart = s.Pos - 1
			answer.ChromEnd = common.MinInt64(answer.ChromStart+fragLength-cigar.NumInsertions(s.Cigar)+cigar.NumDeletions(s.Cigar), reference[answer.Chrom].Size)
			answer.Strand = true
		} else {
			answer.ChromEnd = s.Pos - 1 + cigar.ReferenceLength(s.Cigar)
			answer.Strand = false
			answer.ChromStart = common.MaxInt64(answer.ChromEnd-(fragLength-cigar.NumInsertions(s.Cigar)+cigar.NumDeletions(s.Cigar)), 0)
		}
		return answer
	}
}

func BedScoreToWig(infile string, reference map[string]*chromInfo.ChromInfo) []*wig.Wig {
	wigSlice := make([]*wig.Wig, len(reference))
	var line string
	var chromIndex int
	var midpoint int
	var startNum, endNum, x int64
	var i int = 0
	var doneReading bool = false
	var current *bed.Bed

	//generate Wig skeleton from reference
	for _, v := range reference {
		currentWig := wig.Wig{StepType: "fixedStep", Chrom: v.Name, Start: 1, Step: 1}
		currentWig.Values = make([]*wig.WigValue, v.Size)
		for x = 0; x < v.Size; x++ {
			currentWig.Values[x] = &wig.WigValue{Position: x, Value: 0}
		}
		wigSlice[i] = &currentWig
		i++
	}

	log.Println("Completed wig skeleton, looping through bed.")

	//loop through bed line at a time
	file := fileio.EasyOpen(infile)
	defer file.Close()

	for line, doneReading = fileio.EasyNextRealLine(file); !doneReading; line, doneReading = fileio.EasyNextRealLine(file) {
		words := strings.Split(line, "\t")
		startNum = common.StringToInt64(words[1])
		endNum = common.StringToInt64(words[2])
		current = &bed.Bed{Chrom: words[0], ChromStart: startNum, ChromEnd: endNum}
		if len(words) >= 4 {
			current.Name = words[3]
		}
		if len(words) >= 5 {
			current.Score = common.StringToInt64(words[4])
		}
		chromIndex = getWigChromIndex(current.Chrom, wigSlice)
		midpoint = bedMidpoint(current)
		if wigSlice[chromIndex].Values[midpoint].Value != 0 {
			log.Fatalf("Multiple scores for one position.")
		}

		wigSlice[chromIndex].Values[midpoint].Value = float64(current.Score)

	}
	return wigSlice
}

func BedScoreToWigRange(infile string, reference map[string]*chromInfo.ChromInfo) []*wig.Wig {
	wigSlice := make([]*wig.Wig, len(reference))
	var line string
	var chromIndex int
	var midpoint int
	var startNum, endNum, x int64
	var i int = 0
	var doneReading bool = false
	var current *bed.Bed

	//generate Wig skeleton from reference
	for _, v := range reference {
		currentWig := wig.Wig{StepType: "fixedStep", Chrom: v.Name, Start: 1, Step: 1}
		currentWig.Values = make([]*wig.WigValue, v.Size)
		for x = 0; x < v.Size; x++ {
			currentWig.Values[x] = &wig.WigValue{Position: x, Value: 0}
		}
		wigSlice[i] = &currentWig
		i++
	}

	log.Println("Completed wig skeleton, looping through bed.")

	//loop through bed line at a time
	file := fileio.EasyOpen(infile)
	defer file.Close()

	for line, doneReading = fileio.EasyNextRealLine(file); !doneReading; line, doneReading = fileio.EasyNextRealLine(file) {
		words := strings.Split(line, "\t")
		startNum = common.StringToInt64(words[1])
		endNum = common.StringToInt64(words[2])
		current = &bed.Bed{Chrom: words[0], ChromStart: startNum, ChromEnd: endNum}
		if len(words) >= 4 {
			current.Name = words[3]
		}
		if len(words) >= 5 {
			current.Score = common.StringToInt64(words[4])
		}
		chromIndex = getWigChromIndex(current.Chrom, wigSlice)
		if wigSlice[chromIndex].Values[midpoint].Value != 0 {
			log.Fatalf("Multiple scores for one position.")
		}
		for k := current.ChromStart; k < current.ChromEnd; k++ {
			//fmt.Printf("b[j].Chrom: %s, b[j].ChromStart: %d, b[j].ChromEnd: %d, k: %d, len(wigSlice[chromIndex].Values) %d\n", b[j].Chrom, b[j].ChromStart, b[j].ChromEnd, k, len(wigSlice[chromIndex].Values))
			wigSlice[chromIndex].Values[k+1].Value = float64(current.Score)
		}
	}
	return wigSlice
}

func BedReadsToWig(b []*bed.Bed, reference map[string]*chromInfo.ChromInfo) []*wig.Wig {
	wigSlice := make([]*wig.Wig, len(reference))
	var chromIndex int
	var i int = 0
	var x int64 = 0
	//generate Wig skeleton from reference
	for _, v := range reference {
		currentWig := wig.Wig{StepType: "fixedStep", Chrom: v.Name, Start: 1, Step: 1}
		currentWig.Values = make([]*wig.WigValue, v.Size)
		for x = 0; x < v.Size; x++ {
			currentWig.Values[x] = &wig.WigValue{Position: x, Value: 0}
		}
		wigSlice[i] = &currentWig
		i++
	}

	for j := 0; j < len(b); j++ {
		chromIndex = getWigChromIndex(b[j].Chrom, wigSlice)
		fmt.Printf("b[j].Chrom: %s, b[j].ChromStart: %d, b[j].ChromEnd: %d, j: %d, len(wigSlice[chromIndex].Values) %d\n", b[j].Chrom, b[j].ChromStart, b[j].ChromEnd, j, len(wigSlice[chromIndex].Values))
		for k := b[j].ChromStart; k < b[j].ChromEnd; k++ {
			//fmt.Printf("b[j].Chrom: %s, b[j].ChromStart: %d, b[j].ChromEnd: %d, k: %d, len(wigSlice[chromIndex].Values) %d\n", b[j].Chrom, b[j].ChromStart, b[j].ChromEnd, k, len(wigSlice[chromIndex].Values))
			wigSlice[chromIndex].Values[k].Value++
		}
		//fmt.Printf("b[j].Chrom: %s, b[j].ChromStart: %d, b[j].ChromEnd: %d, j: %d, len(wigSlice[chromIndex].Values) %d\n", b[j].Chrom, b[j].ChromStart, b[j].ChromEnd, j, len(wigSlice[chromIndex].Values))
	}
	return wigSlice
}

func bedMidpoint(b *bed.Bed) int {
	return int(b.ChromEnd+b.ChromStart) / 2
}

func getWigChromIndex(s string, wigSlice []*wig.Wig) int {
	for i := 0; i < len(wigSlice); i++ {
		if s == wigSlice[i].Chrom {
			return i
		}
	}
	log.Fatalf("Bed Chromosome, %s, not in reference genome.", s)
	return -1
}

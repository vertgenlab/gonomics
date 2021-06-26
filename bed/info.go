package bed

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
)

//UngappedRegionsThresholdFromFa finds all regions with at least 'threshold' bases out of 'windowSize' that are ungapped (not '-' or 'N')
func UngappedRegionsThresholdFromFa(fa fasta.Fasta, windowSize int, threshold int) []Bed {
	var inUnGappedRegion bool = false
	var answer []Bed
	var curr Bed
	var refCounter int = 0
	var regionCounter int = 0
	for alnCounter := 0; alnCounter < len(fa.Seq) - windowSize; alnCounter++ {
		if fa.Seq[alnCounter] == '-'{
			continue
		} else if inUnGappedRegion {
			if unGappedThresholdCheck(fa, alnCounter, windowSize, threshold) {
				curr.ChromEnd = refCounter + windowSize
			} else {
				inUnGappedRegion = false
				answer = append(answer, curr)
			}
			refCounter++
		} else {
			if unGappedThresholdCheck(fa, alnCounter, windowSize, threshold) {
				curr = Bed{fa.Name, refCounter, refCounter + windowSize, fmt.Sprintf("Region_%d", regionCounter), 0, None, 7, nil}
				regionCounter++
				inUnGappedRegion = true
			}
			refCounter++
		}
	}
	if inUnGappedRegion { //if we ended in an ungapped region, we need to add the last one to the output slice.
		answer = append(answer, curr)
	}
	return answer
}

func unGappedThresholdCheck(fa fasta.Fasta, alnCounter int, windowSize int, threshold int) bool {
	var unGappedCounter int = 0
	for i := alnCounter; i < alnCounter + windowSize; i++ {
		if !dna.DefineBase(fa.Seq[i]) {
			unGappedCounter++
		}
	}
	return unGappedCounter > threshold
}

//UngappedRegionsFromFa: finds all regions outside gaps in a given fasta record
func UngappedRegionsFromFa(fa fasta.Fasta) []Bed {
	var answer []Bed
	var inRegion bool = false
	var startIndex, index int = 0, 0
	for index = range fa.Seq {
		if dna.DefineBase(fa.Seq[index]) && inRegion == false {
			inRegion = true
			startIndex = index
		} else if !(dna.DefineBase(fa.Seq[index])) && inRegion == true {
			answer = append(answer, Bed{Chrom: fa.Name, ChromStart: startIndex, ChromEnd: index, FieldsInitialized: 3})
			inRegion = false
		}
	}
	if inRegion == true {
		answer = append(answer, Bed{Chrom: fa.Name, ChromStart: startIndex, ChromEnd: len(fa.Seq), FieldsInitialized: 3})
	}
	return answer
}

//UngappedRegionsAllFromFa: Finds ungapped regions or bases that do not contain Ns. Returns a slice of bed records.
func UngappedRegionsAllFromFa(records []fasta.Fasta) []Bed {
	var answer []Bed
	var idx int = 0
	for idx = range records {
		answer = append(answer, UngappedRegionsFromFa(records[idx])...)
	}
	return answer
}

//TotalSize gives back to total region covered by bed entry.
func TotalSize(b []Bed) int {
	var ans, curLen int
	for i := 0; i < len(b); i++ {
		curLen = b[i].ChromEnd - b[i].ChromStart
		ans += curLen
	}
	return ans
}

//Splits fasta regions by using bed regions and concatenate fasta sequences by filling 100 Ns in between
func MakeContigFromBed(fa *fasta.Fasta, beds []Bed) *fasta.Fasta {
	var ans *fasta.Fasta = &fasta.Fasta{Name: fa.Name, Seq: make([]dna.Base, 0)}
	for i, b := range beds {
		ans.Seq = append(ans.Seq, fa.Seq[b.ChromStart:b.ChromEnd]...)
		//adds 100n in between bed regions
		if i < len(beds)-2 {
			ans.Seq = append(ans.Seq, dna.CreateAllNs(100)...)
		}
	}
	return ans
}

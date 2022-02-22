package fasta

import (
	"github.com/vertgenlab/gonomics/dna"
	"log"
)

//RefPosToAlnPos returns the alignment position associated with a given reference position for an input MultiFa. 0 based.
func RefPosToAlnPos(record Fasta, RefPos int) int {
	return RefPosToAlnPosCounter(record, RefPos, 0, 0)
}

//RefPosToAlnPosCounter is like RefPosToAlnPos, but can begin midway through a chromosome at a refPosition/alnPosition pair, defined by the input variables refStart and alnStart.
func RefPosToAlnPosCounter(record Fasta, RefPos int, refStart int, alnStart int) int {
	if refStart > RefPos {
		refStart, alnStart = 0, 0 //in case the refStart was improperly set (greater than the desired position, we reset these counters to 0.
	}
	for t := alnStart; refStart < RefPos; alnStart++ {
		t++
		if t == len(record.Seq) {
			log.Fatalf("Ran out of chromosome.")
		} else if record.Seq[t] != dna.Gap {
			refStart++
		}
	}
	return alnStart
}

//AlnPosToRefPos returns the reference position associated with a given AlnPos for an input Fasta. If the AlnPos corresponds to a gap, it gives the preceeding reference position.
//0 based.
func AlnPosToRefPos(record Fasta, AlnPos int) int {
	var RefPos int = 0
	for t := 0; t < AlnPos; t++ {
		if t == len(record.Seq) {
			log.Fatalf("Ran out of chromosome.")
		} else if record.Seq[t] != dna.Gap {
			RefPos++
		}
	}
	return RefPos
}

//CopySubset returns a copy of a multiFa from a specified start and end position.
func CopySubset(records []Fasta, start int, end int) []Fasta {
	c := make([]Fasta, len(records))
	length := end - start
	for i := 0; i < len(records); i++ {
		c[i] = Fasta{Name: records[i].Name}
		c[i].Seq = make([]dna.Base, length)
		copy(c[i].Seq, records[i].Seq[start:end])
	}
	return c
}

//RemoveMissingMult removes any entries comprised only of gaps in a multiple alignment block,.
func RemoveMissingMult(records []Fasta) []Fasta {
	var answer []Fasta
	var missing bool = true

	for i := 0; i < len(records); i++ {
		missing = true
		for j := 0; j < len(records[i].Seq) && missing; j++ {
			if records[i].Seq[j] != dna.Gap {
				missing = false
			}
		}
		if !missing {
			answer = append(answer, records[i])
		}
	}
	return answer
}

//returns alignment columns with no gaps or lowercase letters
func DistColumn(records []Fasta) []Fasta {
	var subFa = make([]Fasta, len(records))
	for i := 0; i < len(records); i++ {
		subFa[i] = Fasta{Name: records[i].Name, Seq: make([]dna.Base, 0)}
	}

	for i := 0; i < len(records[0].Seq); i++ {
		currentBase := records[0].Seq[i]
		allValid := true
		if !(dna.IsLower(currentBase) || currentBase == dna.Gap) {
			for j := 1; j < len(records); j++ {
				if records[j].Seq[i] == dna.Gap || dna.IsLower(records[j].Seq[i]) {
					allValid = false
				}
			}
		} else {
			allValid = false
		}
		if allValid {
			for k := 0; k < len(records); k++ {
				subFa[k].Seq = append(subFa[k].Seq, records[k].Seq[i])
			}
		}
	}
	return subFa
}

//This function takes in a multiFa alignment block and returns only the columns that contain segregating sites.
func SegregatingSites(aln []Fasta) []Fasta {
	var answer = make([]Fasta, len(aln))
	for i := 0; i < len(aln); i++ {
		answer[i] = Fasta{Name: aln[i].Name, Seq: make([]dna.Base, 0)}
	}
	var current dna.Base
	var isSegregating bool
	for i := 0; i < len(aln[0].Seq); i++ {
		current = aln[0].Seq[i]
		isSegregating = false
		for j := 1; j < len(aln); j++ {
			if aln[j].Seq[i] != current {
				isSegregating = true
			}
		}
		if isSegregating {
			for k := 0; k < len(aln); k++ {
				answer[k].Seq = append(answer[k].Seq, aln[k].Seq[i])
			}
		}
	}
	return answer
}

//NumSegregatingSites returns the number of sites in an alignment block that are segregating.
func NumSegregatingSites(aln []Fasta) int {
	if len(SegregatingSites(aln)) == 0 {
		return 0
	}
	return len(SegregatingSites(aln)[0].Seq)
}

//PairwiseMutationDistanceReferenceWindow takes two input fasta sequences and calculates the number of mutations in a reference window of a given size. Segregating sites are counted as 1, as are INDELs regardless of length.
//alnStart indicates the beginning alignment column for distance evaluation, and windowSize is the number of references bases to compare.
//Three returns, first is the pairwise mutation distance, second is reachedEnd, a bool that is true for incomplete windows. The third return is alignmentEnd, or the last alignment column evaluated.
func PairwiseMutationDistanceReferenceWindow(seq1 Fasta, seq2 Fasta, alnStart int, windowSize int) (int, bool, int) {
	diff := 0
	baseCount := 0
	var seq1Indel bool = false
	var seq2Indel bool = false
	var reachedEnd bool = false
	var i int = 0

	for i = alnStart; baseCount < windowSize && i < len(seq1.Seq); i++ {
		if seq1.Seq[i] == seq2.Seq[i] {
			if seq1.Seq[i] != dna.Gap {
				seq1Indel = false
				seq2Indel = false
				baseCount++
			}
		} else if seq1.Seq[i] == dna.Gap {
			seq2Indel = false
			if !seq1Indel {
				seq1Indel = true
				diff++
			}
		} else if seq2.Seq[i] == dna.Gap {
			baseCount++
			seq1Indel = false
			if !seq2Indel {
				seq2Indel = true
				diff++
			}
		} else if seq1.Seq[i] != seq2.Seq[i] {
			seq1Indel = false
			seq2Indel = false
			baseCount++
			diff++
		} else {
			log.Fatalf("Something went horribly wrong.")
		}
	}

	if baseCount != windowSize {
		reachedEnd = true
	}
	return diff, reachedEnd, i
}

//PairwiseMutationDistanceInRange calculates the number of mutations between two Fasta sequences from a specified
//start and end alignment column. Segregating sites are counted as 1, as are INDELs regardless of length.
func PairwiseMutationDistanceInRange(seq1 Fasta, seq2 Fasta, alnStart int, alnEnd int) int {
	diff := 0
	var seq1Indel bool = false
	var seq2Indel bool = false

	if alnEnd >= len(seq1.Seq) {
		log.Fatalf("Error in PairwiseMutationDistanceInRange, alnEnd must be less than the length of the sequence.")
	}

	for i := alnStart; i < alnEnd; i++ {
		if seq1.Seq[i] == seq2.Seq[i] {
			if seq1.Seq[i] != dna.Gap {
				seq1Indel = false
				seq2Indel = false
			}
		} else if seq1.Seq[i] == dna.Gap {
			seq2Indel = false
			if !seq1Indel {
				seq1Indel = true
				diff++
			}
		} else if seq2.Seq[i] == dna.Gap {
			seq1Indel = false
			if !seq2Indel {
				seq2Indel = true
				diff++
			}
		} else if seq1.Seq[i] != seq2.Seq[i] {
			seq1Indel = false
			seq2Indel = false
			diff++
		} else {
			log.Fatalf("Something went horribly wrong.")
		}
	}

	return diff
}

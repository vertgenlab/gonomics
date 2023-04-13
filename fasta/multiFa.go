package fasta

import (
	"log"

	"github.com/vertgenlab/gonomics/dna"
)

// RefPosToAlnPos returns the alignment position associated with a given reference position for an input MultiFa. 0 based.
func RefPosToAlnPos(record Fasta, RefPos int) int {
	return RefPosToAlnPosCounter(record, RefPos, 0, 0)
}

// RefPosToAlnPosCounter is like RefPosToAlnPos, but can begin midway through a chromosome at a refPosition/alnPosition pair, defined by the input variables refStart and alnStart.
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

// AlnPosToRefPos returns the reference position associated with a given AlnPos for an input Fasta. If the AlnPos corresponds to a gap, it gives the preceding reference position.
// 0 based.
func AlnPosToRefPos(record Fasta, AlnPos int) int {
	return AlnPosToRefPosCounter(record, AlnPos, 0, 0)
}

// AlnPosToRefPosCounter is like AlnPosToRefPos, but can begin midway through a chromosome at a refPosition/alnPosition pair, defined with the input variables refStart and alnStart.
func AlnPosToRefPosCounter(record Fasta, AlnPos int, refStart int, alnStart int) int {
	if alnStart > AlnPos {
		refStart, alnStart = 0, 0 //in case the alnStart was improperly set (greater than the desired position, we reset the counters to 0.
	}
	for t := alnStart; t < AlnPos; t++ {
		if t == len(record.Seq) {
			log.Fatalf("Ran out of chromosome.")
		} else if record.Seq[t] != dna.Gap {
			refStart++
		}
	}
	return refStart
}

// CopySubset returns a copy of a multiFa from a specified start and end position.
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

// RemoveMissingMult removes any entries comprised only of gaps in a multiple alignment block,.
func RemoveMissingMult(records []Fasta) []Fasta {
	var answer []Fasta
	var missing bool
	var i, j int
	for i = 0; i < len(records); i++ {
		missing = true
		for j = 0; j < len(records[i].Seq) && missing; j++ {
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

// returns alignment columns with no gaps or lowercase letters
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

// emptyCopy returns a new alignment where the sequences have the same names as the input
// alignment, but empty sequences.
func emptyCopy(aln []Fasta) []Fasta {
	var answer []Fasta = make([]Fasta, len(aln))
	for i := range aln {
		answer[i].Name = aln[i].Name
	}
	return answer
}

// isSegregating returns false if the value of all bases in the column (colIdx) are of
// equal value, and true otherwise.
func isSegregating(aln []Fasta, colIdx int) bool {
	var i int
	var firstBase dna.Base

	firstBase = aln[0].Seq[colIdx]
	for i = 1; i < len(aln); i++ {
		if aln[i].Seq[colIdx] != firstBase {
			return true
		}
	}
	return false
}

// SegregatingSites takes in a multiFa alignment and returns a new alignment containing only the columns with segregating sites.
func SegregatingSites(aln []Fasta) []Fasta {
	var answer []Fasta = emptyCopy(aln)
	var i, k int
	for i = 0; i < len(aln[0].Seq); i++ {
		if isSegregating(aln, i) {
			for k = 0; k < len(aln); k++ {
				answer[k].Seq = append(answer[k].Seq, aln[k].Seq[i])
			}
		}
	}
	return answer
}

// NumSegregatingSites returns the number of sites in an alignment block that are segregating.
func NumSegregatingSites(aln []Fasta) int {
	if len(SegregatingSites(aln)) == 0 {
		return 0
	}
	return len(SegregatingSites(aln)[0].Seq)
}

// PairwiseMutationDistanceReferenceWindow takes two input fasta sequences and calculates the number of mutations in a reference window of a given size. Segregating sites are counted as 1, as are INDELs regardless of length.
// alnStart indicates the beginning alignment column for distance evaluation, and windowSize is the number of references bases to compare.
// Three returns, first is the pairwise mutation distance, second is reachedEnd, a bool that is true for incomplete windows. The third return is alignmentEnd, or the last alignment column evaluated.
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

// PairwiseMutationDistanceInRange calculates the number of mutations between two Fasta sequences from a specified
// start and end alignment column. Segregating sites are counted as 1, as are INDELs regardless of length.
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

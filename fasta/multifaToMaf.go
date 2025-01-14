package fasta

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/maf"
	"log"
)

// ToMaf converts a multi-fasta alignment to a slice of maf blocks
// referene must be in the alignment
func ToMaf(aln []fasta.Fasta, refName string) []*maf.Maf {
	var maf []*maf.Maf
	refIndex := findRefIndex(aln, ref)
	
	// iterate over reference to identify blocks
	var start, size int
	inBlock := false
	for i := 0; i < len(aln[refIndex].Seq); i++ {
		if aln[refIndex].Seq[i] != dna.Gap {
			if !inBlock { // start new block
				start = i
				inBlock = true
			}
			size++
		} else if inBlock { // is gap, end of block
			maf = append(maf, createMafBlock(aln, refIndex, start, size))
			inBlock = false
			size = 0
		}
	}

	// if no gap at end of alignment, add last block
	if inBlock {
		maf = append(maf, createMafBlock(aln, refIndex, start, size))
	}

	return maf
}

// findRefIndex returns the index of the reference in the multi-alignment
func findRefIndex(aln []fasta.Fasta, refName string) int {
	for i, record := range aln {
		if record.Name == refName {
			return i
		}
	}
	log.Fatalf("Error: reference name %s not found in the alignment\n", refName)
	return -1
}

// createMafBlock creates a maf from the multi-alignment and the given range.
func createMafBlock(aln []fasta.Fasta, refIndex int, start int, size int) *maf.Maf {
	block := maf.Maf{
		Score:   0.0, // Set score to 0.0 as a placeholder
		Species: make([]*maf.MafSpecies, len(aln)),
	}

	for i, record := range aln {
		subSeq := dna.Extract(record.Seq, start, start+size)
		sLine := maf.MafSLine{
			Src:     record.Name,
			Start:   start,
			Size:    size,
			Strand:  true, // Todo: should I assume this is true?
			SrcSize: len(record.Seq),
			Seq:     subSeq,
		}
		block.Species[i] = maf.MafSpecies{Src: record.Name, SLine: sLine}
	}

	return block
}
package maf

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
)

// MultifaToMaf converts a multi-fasta alignment to a slice of maf blocks
// reference must be in the alignment (can i get rid of this comment)
// speciesChr must be at 1 or length of alignment
// TODO: check at end of string if ends in .txt find the file else take as string for default chrom
// TODO: try to use the expected to run on trails
// TODO: fix expected file
func MultifaToMaf(aln []fasta.Fasta, speciesGiven bool, speciesChr []string) []*Maf {
	// version=1 
	var output []*Maf = make([]*Maf, 0)
	
	// TODO: iterate over reference to identify blocks
	// TODO: don't stop at the gaps, just assume that the whole thing is 1 giant maf block

	block := Maf{
		Score:   0.0, // Todo: set score to 0.0 as a placeholder, check that this is legal with trails
		Species: make([]*MafSpecies, len(aln)),
	}

	var unmaskedCount int
	var maskedCount int
	var srcName string
	for i, record := range aln {
		unmaskedCount, maskedCount, _ = dna.CountMask(record.Seq)
		if speciesGiven {
			srcName = record.Name + "." + speciesChr[i]
		} else {
			srcName = record.Name + ".defaultChr" // TODO: need default chromosome?
		}
		sLine := &MafSLine{
			Src:     srcName,
			Start:   0, // hardcoded
			Size:    unmaskedCount + maskedCount,
			Strand:  true,
			SrcSize: 1000, // TODO: do we want an argument to input the actual src size? len(record.Seq) or some default number 1000
			Seq:     record.Seq,
		}
		block.Species[i] = &MafSpecies{Src: record.Name, SLine: sLine}
	}

	output = append(output, &block)

	return output
}
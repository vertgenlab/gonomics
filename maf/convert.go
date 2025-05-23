package maf

import (
	"log"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
)

func setupBlankAln(ref fasta.Fasta, species []string) []fasta.Fasta {
	var answer []fasta.Fasta
	answer = append(answer, ref)
	answer[0].Name = species[0]
	for i := 1; i < len(species); i++ {
		answer = append(answer, fasta.CreateAllGaps(species[i], len(ref.Seq)))
	}
	return answer
}

func insertMafBlockIntoFasta(aln []fasta.Fasta, m *Maf) []fasta.Fasta {
	var speciesBlock *MafSpecies
	var replaceStart, replaceEnd, replaceAlnLength int
	if m == nil || len(m.Species) < 1 || len(aln) < 1 {
		log.Fatalf("Error: empty maf or fasta alignment passed into insertMafBlockIntoFasta\n")
	}

	//first deal with reference species to get replace coordinates
	fastaRefName := aln[0].Name
	mafRefAssembly, mafRefChrom := SrcToAssemblyAndChrom(m.Species[0].Src)
	if mafRefAssembly != fastaRefName && mafRefChrom != fastaRefName {
		log.Fatalf("Error: name of reference/first species in maf block is %s, yet name of first species in fasta is %s\n", m.Species[0].Src, fastaRefName)
	}
	refSLine := m.Species[0].SLine
	if refSLine == nil {
		log.Fatalf("Error: did not find a SLine for reference species in maf\n")
	}
	replaceStart = refSLine.Start
	replaceEnd = refSLine.Start + refSLine.Size
	replaceAlnLength = len(refSLine.Seq)

	//now replace what is in each fasta record for the maf interval
	for i := 0; i < len(aln); i++ {
		speciesBlock = FindSpeciesBeforeDot(m, aln[i].Name)
		if i == 0 {
			if dna.CompareSeqsIgnoreCaseAndGaps(speciesBlock.SLine.Seq, aln[i].Seq[replaceStart:replaceEnd]) != 0 {
				log.Fatalf("Error: reference sequence in maf does not match that in the fasta\n"+
					"%s\n"+
					"%s\n"+
					"%d %d\n",
					dna.BasesToString(speciesBlock.SLine.Seq),
					dna.BasesToString(aln[i].Seq[replaceStart:replaceEnd]),
					replaceStart, replaceEnd)
			}
		}
		if speciesBlock == nil || speciesBlock.SLine == nil {
			aln[i].Seq = dna.Replace(aln[i].Seq, replaceStart, replaceEnd, dna.CreateAllGaps(replaceAlnLength))
		} else {
			aln[i].Seq = dna.Replace(aln[i].Seq, replaceStart, replaceEnd, speciesBlock.SLine.Seq)
		}
	}
	return aln
}

// ToFasta converts a slice of maf blocks to a multi-fasta alignment.
// Along with the slice of maf blocks, it also needs the reference genome/chromosome
// as an input and a slice of the other species/assemblies to include
// in the multi-fasta output alignment
func ToFasta(m []*Maf, ref fasta.Fasta, species []string) []fasta.Fasta {
	if len(ref.Seq) != m[0].Species[0].SLine.SrcSize {
		log.Fatalf("Error: ref seq supplied as fasta should match the length of the first seq in the first maf block\n")
	}
	aln := setupBlankAln(ref, species)
	SortByPosRev(m)
	for i := 0; i < len(m); i++ {
		insertMafBlockIntoFasta(aln, m[i])
	}
	return aln
}

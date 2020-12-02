package goGene

import (
	"errors"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/vcf"
)

type Feature int32
type MutationType byte

const (
	Intron   Feature = -1
	UtrThree Feature = -3
	UtrFive  Feature = -5
	// All positive values refer to cDNA position

	Silent     MutationType = 0
	Missense   MutationType = 1
	Nonsense   MutationType = 2
	Frameshift MutationType = 3
	Intergenic MutationType = 4
	Intronic   MutationType = 5
	Splice     MutationType = 6 // +/- 1-2 of E-I boundary
	FarSplice  MutationType = 7 // +/- 3-10 of E-I boundary
)

// GoGene is a processed version of a gtf record that enables easy
// traversal and manipulation of genes on the genomic and mRNA levels.
type GoGene struct {
	id           string
	startPos     int
	strand       bool
	cdsStarts    []int
	genomeSeq    []dna.Base
	cdnaSeq      []dna.Base
	featureArray []Feature
}

type EffectPrediction struct {
	vcf.Vcf
	RefId       string
	Consequence MutationType
	CdnaPos     int
	CdnaOffset  int
	AaPos       int
	AaRef       []dna.AminoAcid
	AaAlt       []dna.AminoAcid
}

// GenomicToCdna converts genomic coordinates to cDNA coordinates. The return format is c.100+10.
// The first int return is the nearest position in the coding sequence.
// The second int return is the distance from the nearest coding exon
// (>0 if after cds; <0 if before cds, ==0 if inside coding sequence).
// Input and output positions are zero-based
func GenomicPosToCdna(g *GoGene, genomePos int) (int, int, error) {
	var queryPos int
	if g.strand { // Positive Strand
		queryPos = genomePos - g.startPos
	} else { // Negative Strand
		queryPos = g.startPos - genomePos
	}

	feature := g.featureArray[queryPos]
	switch feature {
	case -1:
		var forwardOffset, reverseOffset int = 1, -1
		for {
			if g.featureArray[queryPos+reverseOffset] > 0 {
				// Note: the offset values are returned with their sign flipped as we are returning
				// the offset distance FROM the cds, not the offset distance TO the cds
				return int(g.featureArray[queryPos+reverseOffset]), reverseOffset * -1, nil
			}
			if g.featureArray[queryPos+forwardOffset] > 0 {
				return int(g.featureArray[queryPos+forwardOffset]), forwardOffset * -1, nil
			}
			forwardOffset++
			reverseOffset--
			if queryPos+forwardOffset > len(g.featureArray) || queryPos+reverseOffset < 0 {
				return 0, 0, errors.New("ERROR: No coding sequence could be found")
			}
		}

	case -3:
		var reverseOffset int = -1
		for g.featureArray[queryPos+reverseOffset] < 0 {
			reverseOffset--
			if queryPos+reverseOffset < 0 {
				return 0, 0, errors.New("ERROR: No coding sequence found before 3'UTR")
			}
		}
		return int(g.featureArray[queryPos+reverseOffset]), reverseOffset * -1, nil

	case -5:
		var forwardOffset int = 1
		for g.featureArray[queryPos+forwardOffset] < 0 {
			forwardOffset++
			if queryPos+forwardOffset > len(g.featureArray) {
				return 0, 0, errors.New("ERROR: No coding sequence found after 5'UTR")
			}
		}
		return int(g.featureArray[queryPos+forwardOffset]), forwardOffset * -1, nil

	default:
		return int(feature), 0, nil
	}
}

// CdnaPosToGenomic converts cDna coordinates to genomic coordinates
// Input and output positions are zero-based
func CdnaPosToGenomic(g *GoGene, cdnaPos int) (int, error) {
	if cdnaPos < 0 {
		return 0, errors.New("ERROR: Input cDNA position must be positive")
	}
	if cdnaPos > len(g.cdnaSeq)-1 {
		return 0, errors.New("ERROR: Input position is greater than the length of the cDNA")
	}
	var searchStartPos int = g.cdsStarts[0]
	for _, val := range g.cdsStarts {
		if int(g.featureArray[val]) > cdnaPos {
			break
		}
		searchStartPos = val
	}
	if g.strand { // Positive Strand
		return searchStartPos + (cdnaPos - int(g.featureArray[searchStartPos])) + g.startPos, nil
	} else { // Negative Strand
		return g.startPos - (searchStartPos + (cdnaPos - int(g.featureArray[searchStartPos]))), nil
	}
}

//WIP
func VariantEffect(g *GoGene, v *vcf.Vcf) EffectPrediction {
	var answer EffectPrediction

	return answer
}

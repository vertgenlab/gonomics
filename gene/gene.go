package gene

import (
	"errors"
	"github.com/vertgenlab/gonomics/dna"
)

type Feature int32
type MutationType byte

const (
	Intron   Feature = -1
	UtrThree Feature = -3
	UtrFive  Feature = -5
	// All positive values refer to cDNA position

	Silent           MutationType = 0
	Missense         MutationType = 1
	Nonsense         MutationType = 2
	Frameshift       MutationType = 3
	Intergenic       MutationType = 4
	Intronic         MutationType = 5
	Splice           MutationType = 6 // +/- 1-2 of E-I boundary
	FarSplice        MutationType = 7 // +/- 3-10 of E-I boundary
	DisruptStart     MutationType = 8
	DisruptStop      MutationType = 9
	InFrameInsertion MutationType = 10
	InFrameDeletion  MutationType = 11
)

// Gene is a processed version of a gtf record that enables easy
// traversal and manipulation of genes on the genomic and mRNA levels.
type Gene struct {
	id           string          // Identifier for the transcript the Gene is derived from. In GTF this is the GeneID field.
	startPos     int             // Genomic start position of the gene. This should use the coordinate system of the reference fasta, rather than Gene internal genomic coordinates
	posStrand    bool            // True if gene is on the positive strand, False if on negative strand.
	cdsStarts    []int           // The start position of each CDS. This value is stored in Gene genomic coordinates (slice index of genomeSeq)
	cdsEnds      []int           // The end position of each CDS. This value is stored in Gene genomic coordinates (slice index of genomeSeq)
	genomeSeq    []dna.Base      // The genomic sequence of the gene from 5' to 3'. NOTE: This field is reverse complemented relative to reference file if the gene is on the negative strand
	cdnaSeq      []dna.Base      // The cDNA sequence of the gene from 5' to 3'.
	protSeq      []dna.AminoAcid // The polypeptide sequence resulting from the cDNA.
	featureArray []Feature       // FeatureArray is a slice with len(featureArray) = len(genomeSeq). The index of featureArray corresponds to the same index of genomeSeq. featureArray denotes the features listed above as negative values, or the cDNA pos using all values > 0.
	orig         goGeneBackup    // Copy of initial Gene state to enable the Reset() function.
	changeLog    []diff          // Log of any mutations that have been performed on the Gene to enable the Reset() function.
}

// goGeneBackup stores the initial state of a GoGene to enable Reset() functionality.
type goGeneBackup struct {
	startPos     int
	cdsStarts    []int
	cdsEnds      []int
	genomeSeq    []dna.Base
	cdnaSeq      []dna.Base
	featureArray []Feature
}

// diff acts as an entry in a changelog listing how the sequence has been manipulated.
type diff struct {
	genomePos int
	removed   []dna.Base
	added     []dna.Base
}

// EffectPrediction outputs the effects of a mutation on the cDNA and protein sequences.
type EffectPrediction struct {
	Consequence MutationType    // Classification of mutation (see above for values)
	CdnaPos     int             // Base-zero position in the cDNA
	CdnaDist    int             // Distance from nearest CDS. Zero if in a CDS, >0 if 3' of CDS, <0 if 5' of CDS
	AaPos       int             // Base-zero position of first changed amino acid
	AaRef       []dna.AminoAcid // Slice of Ref amino acids (removed from protein)
	AaAlt       []dna.AminoAcid // Slice of Alt amino acids (added to protein)
	StopDist    int             // Distance to stop codon. This value is filled if and only if it changes as a result of the mutation. Value is -1 if unchanged and -2 if the new stop codon is beyond the bounds of the original gene.
}

// GenomicToCdna converts genomic coordinates to cDNA coordinates. The return format is c.100+10 (HGVS).
// The first int return is the nearest position in the coding sequence in cDNA coordinates.
// The second int return is the distance from the nearest coding exon
// (>0 if 5' of cds; <0 if 3' of cds, ==0 if inside coding sequence, ties break to <0).
// Input and output positions are zero-based
func GenomicPosToCdna(g *Gene, genomePos int) (int, int, error) {
	var queryPos int
	if g.posStrand { // Positive Strand
		queryPos = genomePos - g.startPos
	} else { // Negative Strand
		queryPos = g.startPos - genomePos
	}

	feature := g.featureArray[queryPos]
	switch feature {
	case Intron:
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
				return 0, 0, errors.New("no coding sequence could be found")
			}
		}

	case UtrThree:
		var reverseOffset int = -1
		for g.featureArray[queryPos+reverseOffset] < 0 {
			reverseOffset--
			if queryPos+reverseOffset < 0 {
				return 0, 0, errors.New("no coding sequence found before 3'UTR")
			}
		}
		return int(g.featureArray[queryPos+reverseOffset]), reverseOffset * -1, nil

	case UtrFive:
		var forwardOffset int = 1
		for g.featureArray[queryPos+forwardOffset] < 0 {
			forwardOffset++
			if queryPos+forwardOffset > len(g.featureArray) {
				return 0, 0, errors.New("no coding sequence found after 5'UTR")
			}
		}
		return int(g.featureArray[queryPos+forwardOffset]), forwardOffset * -1, nil

	default:
		return int(feature), 0, nil
	}
}

// CdnaPosToGenomic converts cDna coordinates to genomic coordinates
// Input and output positions are zero-based
func CdnaPosToGenomic(g *Gene, cdnaPos int) (int, error) {
	if cdnaPos < 0 {
		return 0, errors.New("input cDNA position must be positive")
	}
	if cdnaPos > len(g.cdnaSeq)-1 {
		return 0, errors.New("input position is greater than the length of the cDNA")
	}
	var searchStartPos int = g.cdsStarts[0]
	for _, val := range g.cdsStarts {
		if int(g.featureArray[val]) > cdnaPos {
			break
		}
		searchStartPos = val
	}
	if g.posStrand { // Positive Strand
		return searchStartPos + (cdnaPos - int(g.featureArray[searchStartPos])) + g.startPos, nil
	} else { // Negative Strand
		return g.startPos - (searchStartPos + (cdnaPos - int(g.featureArray[searchStartPos]))), nil
	}
}

func CdnaPosToCodon(g *Gene, cdnaPos int) (dna.Codon, error) {
	var answer dna.Codon
	if cdnaPos < 0 {
		return answer, errors.New("input cDNA position must be positive")
	}
	if cdnaPos > len(g.cdnaSeq)-1 {
		return answer, errors.New("input position is greater than the length of the cDNA")
	}

	switch cdnaPos % 3 {
	case 0:
		return dna.Codon{Seq: g.cdnaSeq[cdnaPos : cdnaPos+3]}, nil
	case 1:
		return dna.Codon{Seq: g.cdnaSeq[cdnaPos-1 : cdnaPos+2]}, nil
	case 2:
		return dna.Codon{Seq: g.cdnaSeq[cdnaPos-2 : cdnaPos+1]}, nil
	default: // never used
		return answer, errors.New("problem determining frame")
	}
}

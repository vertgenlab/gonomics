package goGene

import (
	"errors"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/numbers"
)

// PointMutation changes a single nucleotide to the desired base, predicts the effect,
// and updates the GoGene struct to reflect the change.
// The position of the mutation should be given in genomic coordinates.
func PointMutation(g *GoGene, genomePos int, alt dna.Base) (EffectPrediction, error) {
	var answer EffectPrediction
	var log diff
	var err error

	if alt != dna.A && alt != dna.C && alt != dna.G && alt != dna.T {
		return answer, errors.New("alt base must be A, C, T, or G")
	}

	if genomePos < 0 {
		return answer, errors.New("genomePos must be positive")
	}
	if g.strand {
		if genomePos < g.startPos {
			return answer, errors.New("input genomePos is not in the gene")
		}
	} else {
		if genomePos > g.startPos {
			return answer, errors.New("input genomePos is not in the gene")
		}
	}

	genomeIndexPos := numbers.Abs(genomePos - g.startPos) // abs needed to handle negative strand
	log.genomeIndexPos = genomeIndexPos
	log.removed = []dna.Base{g.genomeSeq[genomeIndexPos]}
	log.added = []dna.Base{alt}

	if genomeIndexPos > len(g.genomeSeq)-1 {
		return answer, errors.New("input genomePos is not in the gene")
	}

	g.genomeSeq[genomeIndexPos] = alt
	g.changeLog = append(g.changeLog, log)

	cdnaIndexPos := int(g.featureArray[genomeIndexPos])

	if cdnaIndexPos >= 0 { // mutation is coding
		answer.CdnaPos = cdnaIndexPos
		answer.AaPos = cdnaIndexPos / 3
		codon, err := CdnaPosToCodon(g, cdnaIndexPos)
		common.ExitIfError(err)
		answer.AaRef = []dna.AminoAcid{dna.TranslateCodon(&codon)}
		g.cdnaSeq[cdnaIndexPos] = alt
		answer.AaAlt = []dna.AminoAcid{dna.TranslateCodon(&codon)}
		if answer.AaRef[0] == answer.AaAlt[0] {
			answer.Consequence = Silent
		} else {
			if answer.AaAlt[0] == dna.Stop {
				answer.Consequence = Nonsense
			} else if answer.AaRef[0] == dna.Stop {
				answer.Consequence = DisruptStop
			} else if answer.AaPos == 0 {
				answer.Consequence = DisruptStart
			} else {
				answer.Consequence = Missense
			}
		}
	} else { // mutation is noncoding
		answer.CdnaPos, answer.CdnaOffset, err = GenomicPosToCdna(g, genomePos)
		common.ExitIfError(err)
		answer.Consequence = checkSplice(g, genomeIndexPos)
	}

	return answer, nil
}

//TODO
// Insertion adds bases to the GoGene, predicts the effect, and updates the GoGene struct to reflect the change.
// The position should be the genomic coordinate of the base directly BEFORE the inserted bases.
func Insertion(g *GoGene, genomePos int, alt []dna.Base) (EffectPrediction, error) {
	var answer EffectPrediction

	return answer, nil
}

//TODO
// Deletion removes bases from the GoGene, predicts the effect, and updates the GoGene struct to reflect the change.
// The input interval should be half-open (BED style). All positions should be given as genomic coordinates.
func Deletion(g *GoGene, genomeStartPos int, genomeEndPos int) (EffectPrediction, error) {
	var answer EffectPrediction

	return answer, nil
}

// Reset reverts all mutations done to a GoGene.
func Reset(g *GoGene) (int, error) {
	var numChangesUndone int
	var done bool
	var err error

	for !done {
		done, err = Undo(g)
		if done {
			break
		}
		if err != nil {
			return numChangesUndone, err
		}
		numChangesUndone++
	}

	return numChangesUndone, nil
}

// Undo reverts the last mutation done to a GoGene.
// bool return is true if changeLog is empty
func Undo(g *GoGene) (bool, error) {
	if len(g.changeLog) == 0 {
		return true, nil
	}

	toUndo := g.changeLog[len(g.changeLog)-1]
	g.changeLog = g.changeLog[:len(g.changeLog)-1]

	switch {
	case len(toUndo.added) == 1 && len(toUndo.added) == len(toUndo.removed): // Point Mutation
		g.genomeSeq[toUndo.genomeIndexPos] = toUndo.removed[0]
		if g.featureArray[toUndo.genomeIndexPos] >= 0 {
			g.cdnaSeq[g.featureArray[toUndo.genomeIndexPos]] = toUndo.removed[0]
		}
	case len(toUndo.added) == 0 && len(toUndo.removed) > 0: // Deletion
		//TODO
	case len(toUndo.removed) == 0 && len(toUndo.added) > 0: // Insertion
		//TODO
	default:
		return false, errors.New("unrecognized diff in changelog")
	}
	return false, nil
}

// checkSplice inputs a non-coding position in GoGene-adjusted genomic coordinates and determines if a variant at this site may affect splicing.
func checkSplice(g *GoGene, genomeIndexPos int) MutationType {
	var nearestExon, idx int
	var val Feature

	// check forwards
	for nearestExon, val = range g.featureArray[genomeIndexPos:] {
		if val >= 0 || nearestExon > 11 {
			break
		}
	}

	if nearestExon <= 2 {
		return Splice
	} else if nearestExon <= 10 {
		return FarSplice
	}

	// check backwards
	nearestExon = 0
	var lowerBound int
	if genomeIndexPos-11 > 0 {
		lowerBound = genomeIndexPos - 11
	}
	for idx, val = range g.featureArray[lowerBound:genomeIndexPos] {
		if val >= 0 {
			nearestExon = idx
		}
	}

	if nearestExon >= 10 {
		return Splice
	} else if nearestExon >= 1 {
		return FarSplice
	}

	return Intronic
}

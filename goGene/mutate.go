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
	log.genomePos = genomePos
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
	var log diff

	for _, val := range alt {
		if val != dna.A && val != dna.C && val != dna.G && val != dna.T {
			return answer, errors.New("all alt bases must be A, C, T, or G")
		}
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

	if genomeIndexPos > len(g.genomeSeq)-1 {
		return answer, errors.New("input genomePos is not in the gene")
	}

	// Fill changeLog
	log.genomePos = genomePos
	log.added = alt
	g.changeLog = append(g.changeLog, log)

	// Update genome sequence
	g.genomeSeq = append(g.genomeSeq, alt...)                                     // just to increase cap, will get overwritten next line
	copy(g.genomeSeq[genomeIndexPos+1+len(alt):], g.genomeSeq[genomeIndexPos+1:]) // shift seq over len(alt) positions
	for idx, val := range alt {                                                   // copy over inserted bases
		g.genomeSeq[genomeIndexPos+1+idx] = val
	}

	// Update CDS starts and ends
	for idx, val := range g.cdsStarts {
		if val <= genomeIndexPos && g.cdsEnds[idx] > genomeIndexPos {
			g.cdsEnds[idx] += len(alt)
		}
		if val > genomeIndexPos {
			g.cdsStarts[idx] += len(alt)
			g.cdsEnds[idx] += len(alt)
		}
	}

	if g.featureArray[genomeIndexPos] >= 0 && g.featureArray[genomeIndexPos+1] >= 0 { // Coding
		// Update cDNA
		cdnaPos := int(g.featureArray[genomeIndexPos])
		g.cdnaSeq = append(g.cdnaSeq, alt...)
		copy(g.cdnaSeq[cdnaPos+1+len(alt):], g.cdnaSeq[cdnaPos+1:])
		for idx, val := range alt {
			g.cdnaSeq[cdnaPos+1+idx] = val
		}

		// Update featureArray for coding
		fillVal := Feature(cdnaPos + 1)
		tmpSlice := make([]Feature, len(alt))
		g.featureArray = append(g.featureArray, tmpSlice...)
		copy(g.featureArray[genomeIndexPos+1+len(alt):], g.featureArray[genomeIndexPos+1:])
		var i int
		for i = 0; i < len(alt); i++ { // fill the inserted elements with the appropriate cDNA pos value
			g.featureArray[genomeIndexPos+1+i] = fillVal
			fillVal++
		}
		for ; g.featureArray[genomeIndexPos+1+i] >= 0; i++ { // update the rest of the current exon
			g.featureArray[genomeIndexPos+1+i] = fillVal
			fillVal++
		}
		for _, cdsStart := range g.cdsStarts { // update the rest of the CDS exons
			if cdsStart > genomeIndexPos {
				for k := 0; g.featureArray[cdsStart+k] >= 0; k++ {
					g.featureArray[cdsStart+k] = fillVal
					fillVal++
				}
			}
		}
	} else { // NonCoding
		// Update featureArray for noncoding
		// figure out how the inserted bases should be labelled in the feature array. Priority: 5'UTR=3'Utr>Intron
		fillVal := Feature(numbers.Min(int(g.featureArray[genomeIndexPos]), int(g.featureArray[genomeIndexPos+1])))
		tmpSlice := make([]Feature, len(alt))
		g.featureArray = append(g.featureArray, tmpSlice...)
		copy(g.featureArray[genomeIndexPos+1+len(alt):], g.featureArray[genomeIndexPos+1:])
		for i := 0; i < len(alt); i++ {
			g.featureArray[genomeIndexPos+1+i] = fillVal
		}
	}

	//TODO EffectPrediction

	return answer, nil
}

//TODO
// Deletion removes bases from the GoGene, predicts the effect, and updates the GoGene struct to reflect the change.
// genomeStartPos should be the base directly BEFORE the deleted bases. genomeEndPos should be the base directly AFTER the deleted bases.
// All positions should be given as genomic coordinates.
func Deletion(g *GoGene, genomeStartPos int, genomeEndPos int) (EffectPrediction, error) {
	var answer EffectPrediction
	var log diff
	var err error

	if genomeStartPos < 0 || genomeEndPos < 0 {
		return answer, errors.New("genome positions must be positive")
	}

	if genomeStartPos >= genomeEndPos {
		return answer, errors.New("genomeStartPos must be less than genomeEndPos")
	}

	if g.strand {
		if genomeStartPos < g.startPos {
			if genomeEndPos > g.startPos {
				genomeStartPos = g.startPos - 1 // if deletion partially overlaps the gene, then trim the value
			} else {
				return answer, errors.New("deletion does not overlap the gene")
			}
		}
	} else {
		if genomeStartPos > g.startPos {
			if genomeEndPos < g.startPos {
				genomeStartPos = g.startPos // if deletion partially overlaps the gene, then trim the value
			} else {
				return answer, errors.New("deletion does not overlap the gene")
			}
		}
	}

	genomeIndexStartPos := numbers.Abs(genomeStartPos - g.startPos) // abs needed to handle negative strand
	genomeIndexEndPos := numbers.Abs(genomeEndPos - g.startPos)
	if genomeIndexEndPos > len(g.genomeSeq) {
		genomeIndexEndPos = len(g.genomeSeq) // if deletion goes past the end of the gene, then trim the value
	}
	if genomeIndexStartPos >= len(g.genomeSeq)-1 {
		return answer, errors.New("deletion does not overlap the gene")
	}

	// Fill changeLog
	log.genomePos = genomeStartPos
	log.removed = make([]dna.Base, genomeIndexEndPos-(genomeIndexStartPos+1))
	copy(log.removed, g.genomeSeq[genomeIndexStartPos+1:genomeIndexEndPos])
	g.changeLog = append(g.changeLog, log)

	// Delete from Genome Sequence
	copy(g.genomeSeq[genomeIndexStartPos+1:], g.genomeSeq[genomeIndexEndPos:])
	g.genomeSeq = g.genomeSeq[:len(g.genomeSeq)-(genomeIndexEndPos-(genomeIndexStartPos+1))]

	// Update startPos if necessary
	if genomeIndexStartPos == -1 {
		g.startPos += genomeIndexEndPos - (genomeIndexStartPos + 1)
	}

	// Update CDS starts and ends and calculate how many cDNA bases have been deleted
	var cdsIndexesToDelete []int
	var deletedCdnaBases int
	var cdnaDelStart, cdnaDelEnd int = -1, -1
	deletionLen := genomeIndexEndPos - (genomeIndexStartPos + 1)
	for i := 0; i < len(g.cdsStarts); i++ {

		switch {
		// if cds is before deletion
		case genomeIndexStartPos >= g.cdsEnds[i]:
			cdnaDelStart = int(g.featureArray[g.cdsEnds[i]])
			continue

		// If cds exon is fully deleted
		case genomeIndexStartPos < g.cdsStarts[i] && genomeIndexEndPos > g.cdsEnds[i]:
			// mark cds for deletion
			cdsIndexesToDelete = append(cdsIndexesToDelete, i)
			deletedCdnaBases += (g.cdsEnds[i] + 1) - g.cdsStarts[i]
			if cdnaDelStart == -1 {
				cdnaDelStart = int(g.featureArray[g.cdsStarts[i]]) - 1
			}
			cdnaDelEnd = int(g.featureArray[g.cdsEnds[i]]) + 1

		// If deletion falls within a cds
		case genomeIndexStartPos >= g.cdsStarts[i] && genomeIndexStartPos < g.cdsEnds[i] &&
			genomeIndexEndPos > g.cdsStarts[i] && genomeIndexEndPos <= g.cdsEnds[i]:
			g.cdsEnds[i] -= deletionLen
			deletedCdnaBases += deletionLen
			cdnaDelStart = int(g.featureArray[genomeIndexStartPos])
			cdnaDelEnd = int(g.featureArray[genomeIndexEndPos])

		// If cds exon is partially deleted on right
		case genomeIndexStartPos >= g.cdsStarts[i] && genomeIndexStartPos < g.cdsEnds[i]:
			deletedCdnaBases += g.cdsEnds[i] - genomeIndexStartPos
			g.cdsEnds[i] = genomeIndexStartPos
			cdnaDelStart = int(g.featureArray[genomeIndexStartPos])

		// If cds exon is partially deleted on left
		case genomeIndexEndPos > g.cdsStarts[i] && genomeIndexEndPos <= g.cdsEnds[i]:
			deletedCdnaBases += genomeIndexEndPos - g.cdsStarts[i]
			g.cdsStarts[i] = genomeIndexEndPos - deletionLen
			g.cdsEnds[i] -= deletionLen
			cdnaDelEnd = int(g.featureArray[genomeIndexEndPos])

		// if cds is after deletion
		case genomeIndexEndPos <= g.cdsStarts[i]:
			if cdnaDelEnd == -1 {
				cdnaDelEnd = g.cdsStarts[i]
			}

			g.cdsStarts[i] -= deletionLen
			g.cdsEnds[i] -= deletionLen
		}
	}

	// Delete removed CDS
	if len(cdsIndexesToDelete) > 0 {
		copy(g.cdsStarts[cdsIndexesToDelete[0]:], g.cdsStarts[cdsIndexesToDelete[len(cdsIndexesToDelete)-1]+1:])
		g.cdsStarts = g.cdsStarts[:len(g.cdsStarts)-len(cdsIndexesToDelete)]
		copy(g.cdsEnds[cdsIndexesToDelete[0]:], g.cdsEnds[cdsIndexesToDelete[len(cdsIndexesToDelete)-1]+1:])
		g.cdsEnds = g.cdsEnds[:len(g.cdsEnds)-len(cdsIndexesToDelete)]
	}

	// Delete cDNA
	if deletedCdnaBases > 0 {
		copy(g.cdnaSeq[cdnaDelStart+1:], g.cdnaSeq[cdnaDelEnd:])
		g.cdnaSeq = g.cdnaSeq[:len(g.cdnaSeq)-deletedCdnaBases]
	}

	// Update Feature Array
	copy(g.featureArray[genomeIndexStartPos+1:], g.featureArray[genomeIndexEndPos:])
	g.featureArray = g.featureArray[:len(g.featureArray)-(genomeIndexEndPos-(genomeIndexStartPos+1))]

	var j, l int
	if g.featureArray[genomeIndexStartPos] >= 0 {
		for j = genomeIndexStartPos + 1; g.featureArray[j] >= 0; j++ {
			g.featureArray[j] -= Feature(deletedCdnaBases)
		}
	}
	for _, val := range g.cdsStarts {
		if val > j {
			for l = val; g.featureArray[l] >= 0; l++ {
				g.featureArray[l] -= Feature(deletedCdnaBases)
			}
		}
	}

	//TODO EffectPrediction

	return answer, err
}

// Reset reverts all mutations done to a GoGene.
func Reset(g *GoGene) {
	var hasDel bool
	var err error
	for _, val := range g.changeLog {
		if len(val.added) == 0 && len(val.removed) > 0 {
			hasDel = true
			break
		}
	}

	if !hasDel {
		for i := len(g.changeLog) - 1; i >= 0; i-- {
			if len(g.changeLog[i].added) == 1 && len(g.changeLog[i].removed) == 1 {
				_, err = PointMutation(g, g.changeLog[i].genomePos, g.changeLog[i].removed[0])
				g.changeLog = g.changeLog[:len(g.changeLog)-2]
				if err != nil {
					hasDel = true
					break
				}
			} else if len(g.changeLog[i].added) > 1 && len(g.changeLog[i].removed) == 0 {
				_, err = Deletion(g, g.changeLog[i].genomePos, g.changeLog[i].genomePos+len(g.changeLog[i].added)+1)
				g.changeLog = g.changeLog[:len(g.changeLog)-2]
				if err != nil {
					hasDel = true
					break
				}
			}
		}
	}

	if hasDel {
		g.startPos = g.orig.startPos
		copy(g.cdsStarts, g.orig.cdsStarts)
		g.cdsStarts = g.cdsStarts[:len(g.orig.cdsStarts)]
		copy(g.cdsEnds, g.orig.cdsEnds)
		g.cdsEnds = g.cdsEnds[:len(g.orig.cdsEnds)]
		copy(g.genomeSeq, g.orig.genomeSeq)
		g.genomeSeq = g.genomeSeq[:len(g.orig.genomeSeq)]
		copy(g.cdnaSeq, g.orig.cdnaSeq)
		g.cdnaSeq = g.cdnaSeq[:len(g.orig.cdnaSeq)]
		copy(g.featureArray, g.orig.featureArray)
		g.featureArray = g.featureArray[:len(g.orig.featureArray)]
	}

	g.changeLog = nil
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

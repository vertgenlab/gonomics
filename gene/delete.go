package gene

import (
	"errors"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/numbers"
)

var (
	ErrNoStopFound = errors.New("frameshift was unable to find a stop before the end of the gene")
)

// Deletion removes bases from the Gene, predicts the effect, and updates the Gene struct to reflect the change.
// genomeStartPos should be the first deleted base. genomeEndPos should be the base directly AFTER the last deleted base. i.e. left closed, right open.
// All positions should be given as base-zero genomic coordinates.
func Deletion(g *Gene, genomeStartPos int, genomeEndPos int) (EffectPrediction, error) {
	var answer EffectPrediction
	answer.StopDist = -1
	var err error

	// The interval genomeIndexStartPos:genomeIndexEndPos is zero base left closed relative to the start
	// of the genomic sequence of the gene.
	genomeIndexStartPos, genomeIndexEndPos, err := deletionPreRunChecks(g, genomeStartPos, genomeEndPos)
	if err != nil {
		return EffectPrediction{}, err
	}

	// determine cDNA position and distance from nearest CDS for effect pred
	var cdnaDistFromDelStart, cdnaDistFromDelEnd int                   // for effect prediction later
	_, cdnaDistFromDelStart, err = GenomicPosToCdna(g, genomeStartPos) // from deletion start
	if g.featureArray[genomeIndexStartPos] > 0 {
		answer.CdnaPos = int(g.featureArray[genomeIndexStartPos])
	} else {
		answer.CdnaDist = numbers.Min(cdnaDistFromDelStart, cdnaDistFromDelEnd)
	}
	_, cdnaDistFromDelEnd, err = GenomicPosToCdna(g, genomeEndPos-1) // from deletion end

	// Save CDS start and stop for later
	origCdsStart := g.cdsStarts[0]
	origCdsEnd := g.cdsEnds[len(g.cdsEnds)-1]

	// Delete from Genome Sequence
	deleteUpdateGenome(g, genomeIndexStartPos, genomeIndexEndPos)

	// Update CDS starts and ends and calculate how many cDNA bases have been deleted
	deletedCodingBases, err := deleteUpdateCds(g, genomeIndexStartPos, genomeIndexEndPos)
	if err != nil {
		return EffectPrediction{}, err
	}

	// Determine UTR bases to delete
	err = deleteUpdateUtr(g, genomeIndexStartPos, genomeIndexEndPos, origCdsStart, origCdsEnd)
	if err != nil {
		return EffectPrediction{}, err
	}

	// Update Feature Array
	deleteUpdateFeatureArray(g, genomeIndexStartPos, genomeIndexEndPos, deletedCodingBases)

	answer, err = deleteEffectPrediction(g, deletedCodingBases, answer)
	if err != nil {
		return EffectPrediction{}, err
	}

	g.protSeq = dna.TranslateSeqToTer(g.codingSeq.seq) // TODO improve efficiency by updating the protein as changes are made
	return answer, err
}

// deletionPreRunChecks performs error checking and appends the diff log prior to making a deletion.
func deletionPreRunChecks(g *Gene, genomeStartPos int, genomeEndPos int) (int, int, error) {
	var log diff

	if genomeStartPos < 0 || genomeEndPos < 0 {
		return -1, -1, ErrNegativeInputValue
	}

	if genomeStartPos >= genomeEndPos {
		return -1, -1, ErrInvalidRange
	}

	if g.posStrand {
		if genomeStartPos < g.startPos {
			if genomeEndPos > g.startPos {
				genomeStartPos = g.startPos - 1 // if deletion partially overlaps the gene, then trim the value
			} else {
				return -1, -1, ErrInputPosNotInGene
			}
		}
	} else {
		if genomeStartPos > g.startPos {
			if genomeEndPos < g.startPos {
				genomeStartPos = g.startPos // if deletion partially overlaps the gene, then trim the value
			} else {
				return -1, -1, ErrInputPosNotInGene
			}
		}
	}

	var genomeIndexStartPos, genomeIndexEndPos int // left closed, right open

	if g.posStrand {
		genomeIndexStartPos = genomeStartPos - g.startPos
		genomeIndexEndPos = genomeEndPos - g.startPos
	} else { // flipped for neg posStrand. start of deletion, and end of deletion must be flipped
		genomeIndexStartPos = g.startPos - (genomeEndPos - 1)
		genomeIndexEndPos = g.startPos - (genomeStartPos - 1)
	}

	if genomeIndexEndPos > len(g.genomeSeq) {
		genomeIndexEndPos = len(g.genomeSeq) // if deletion goes past the end of the gene, then trim the value
	}
	if genomeIndexStartPos > len(g.genomeSeq)-1 {
		return -1, -1, ErrInputPosNotInGene
	}

	// Fill changeLog
	log.genomePos = genomeStartPos
	log.removed = make([]dna.Base, genomeIndexEndPos-genomeIndexStartPos)
	copy(log.removed, g.genomeSeq[genomeIndexStartPos:genomeIndexEndPos])
	g.changeLog = append(g.changeLog, log)
	return genomeIndexStartPos, genomeIndexEndPos, nil
}

// safeDelete performs a deletion using the deleteStable function and adjust slide coordinates for all affected slices
func safeDelete(g *Gene, seq *[]dna.Base, delStart int, delEnd int, offset int) error {
	var err error
	if delStart == delEnd {
		return nil
	}
	err = deleteStable(seq, delStart, delEnd)
	delLen := delEnd - delStart

	// shift delete bases over so they are relative to cdnaSeq
	delStart += offset
	delEnd += offset

	// Fix 5' UTR
	switch {
	case delEnd <= g.utrFive.end: // internal deletion
		g.utrFive.end -= delLen
	case delStart <= g.utrFive.end: // right deletion of 5' UTR
		g.utrFive.end = delStart
	}
	if delStart <= g.utrFive.end { // check if deletion overlaps 5' UTR
		g.utrFive.end = delStart
		g.utrFive.seq = g.cdnaSeq[g.utrFive.start:g.utrFive.end]
	}

	// Fix CDS
	switch {
	case delEnd < g.codingSeq.start: // deletion before CDS and does not overlap
		g.codingSeq.start -= delLen
		g.codingSeq.end -= delLen

	case delStart > g.codingSeq.end: // deletion after CDS and does not overlap
		break // no changes needed

	case delStart <= g.codingSeq.start && delEnd >= g.codingSeq.end: // delete entire CDS
		g.codingSeq.start = delStart
		g.codingSeq.end = delStart

	case g.codingSeq.start <= delStart && g.codingSeq.end >= delEnd: // internal deletion
		g.codingSeq.end -= delLen

	case g.codingSeq.start >= delStart && g.codingSeq.start <= delEnd: // deletion on left
		g.codingSeq.start = delStart
		g.codingSeq.end -= delLen

	case g.codingSeq.end >= delStart && g.codingSeq.end <= delEnd: // deletion on right
		g.codingSeq.end = delStart
	}

	g.codingSeq.seq = g.cdnaSeq[g.codingSeq.start:g.codingSeq.end]

	// Fix 3' UTR
	switch {
	case delEnd < g.utrThree.start: // deletion before 3' UTR and does not overlap
		g.utrThree.start -= delLen
		g.utrThree.end -= delLen

	case delStart <= g.utrThree.start && delEnd >= g.utrThree.end: // completely deletes 3' UTR
		g.utrThree.start = delStart
		g.utrThree.end = delStart

	case delStart >= g.utrThree.start: // internal deletion
		g.utrThree.end = delEnd

	case delStart <= g.utrThree.start: // deletion on left
		g.utrThree.start = delStart
		g.utrThree.end -= delLen
	}

	g.utrThree.seq = g.cdnaSeq[g.utrThree.start:g.utrThree.end]

	g.cdnaSeq = g.cdnaSeq[:len(g.cdnaSeq)-delLen]
	return err
}

var (
	errInvalidPosition = errors.New("input position to deleteStable is invalid")
)

// deleteStable is identical to Delete, but conserves the
// underlying memory of the input destSeq slice.
func deleteStable(seq *[]dna.Base, delStart int, delEnd int) error {
	if delStart >= delEnd || delStart < 0 || delEnd > len(*seq) {
		return errInvalidPosition
	}
	if delStart == delEnd {
		return nil
	}
	*seq = (*seq)[:len(*seq)-(delEnd-delStart)]                // decrease the LEN of the slice to by num of deleted bases
	copy((*seq)[delStart:cap(*seq)], (*seq)[delEnd:cap(*seq)]) // shift bases in slice to remove deleted bases
	return nil
}

// deleteUpdateGenome updates the genomic sequence upon deletion
func deleteUpdateGenome(g *Gene, genomeIndexStartPos int, genomeIndexEndPos int) {
	deletionLen := genomeIndexEndPos - genomeIndexStartPos
	copy(g.genomeSeq[genomeIndexStartPos:], g.genomeSeq[genomeIndexEndPos:])
	g.genomeSeq = g.genomeSeq[:len(g.genomeSeq)-deletionLen]

	// Update startPos if necessary
	if genomeIndexStartPos == 0 {
		g.startPos += genomeIndexEndPos
	}
}

// deleteUpdateCds updates the coding sequence upon deletion
func deleteUpdateCds(g *Gene, genomeIndexStartPos int, genomeIndexEndPos int) (deletedCodingBases int, err error) {
	var cdsIndexesToDelete []int
	var codingDelStart, codingDelEnd int = -1, -1
	deletionLen := genomeIndexEndPos - genomeIndexStartPos
	for i := 0; i < len(g.cdsStarts); i++ {

		switch {
		// if cds is before deletion
		case genomeIndexStartPos > g.cdsEnds[i]:
			codingDelStart = int(g.featureArray[g.cdsEnds[i]]) + 1

		// If cds exon is fully deleted
		case genomeIndexStartPos <= g.cdsStarts[i] && genomeIndexEndPos > g.cdsEnds[i]:
			// mark cds for deletion
			cdsIndexesToDelete = append(cdsIndexesToDelete, i)
			deletedCodingBases += (g.cdsEnds[i] + 1) - g.cdsStarts[i]
			if codingDelStart == -1 {
				codingDelStart = int(g.featureArray[g.cdsStarts[i]])
			}
			codingDelEnd = int(g.featureArray[g.cdsEnds[i]]) + 1

		// If deletion falls within a cds
		case genomeIndexStartPos > g.cdsStarts[i] && genomeIndexStartPos < g.cdsEnds[i] &&
			genomeIndexEndPos > g.cdsStarts[i] && genomeIndexEndPos <= g.cdsEnds[i]:
			g.cdsEnds[i] -= deletionLen
			deletedCodingBases += deletionLen
			codingDelStart = int(g.featureArray[genomeIndexStartPos])
			codingDelEnd = int(g.featureArray[genomeIndexEndPos])

		// If cds exon is partially deleted on right
		case genomeIndexStartPos > g.cdsStarts[i] && genomeIndexStartPos <= g.cdsEnds[i]:
			deletedCodingBases += 1 + g.cdsEnds[i] - genomeIndexStartPos
			g.cdsEnds[i] = genomeIndexStartPos - 1
			codingDelStart = int(g.featureArray[genomeIndexStartPos])

		// If cds exon is partially deleted on left
		case genomeIndexEndPos > g.cdsStarts[i] && genomeIndexEndPos <= g.cdsEnds[i]:
			deletedCodingBases += genomeIndexEndPos - g.cdsStarts[i]
			g.cdsStarts[i] = genomeIndexEndPos - deletionLen
			g.cdsEnds[i] -= deletionLen
			codingDelEnd = int(g.featureArray[genomeIndexEndPos])

		// if cds is after deletion
		case genomeIndexEndPos <= g.cdsStarts[i]:
			if codingDelEnd == -1 {
				codingDelEnd = int(g.featureArray[g.cdsStarts[i]])
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
	if deletedCodingBases > 0 {
		err = safeDelete(g, &g.codingSeq.seq, codingDelStart, codingDelEnd, len(g.utrFive.seq))
	}

	return
}

// deleteUpdateFeatureArray updates the feature array upon deletion.
// CDS must be updated prior to calling this function.
func deleteUpdateFeatureArray(g *Gene, genomeIndexStartPos int, genomeIndexEndPos int, deletedCodingBases int) {
	copy(g.featureArray[genomeIndexStartPos:], g.featureArray[genomeIndexEndPos:])
	g.featureArray = g.featureArray[:len(g.featureArray)-(genomeIndexEndPos-(genomeIndexStartPos))]

	j := genomeIndexStartPos
	if genomeIndexStartPos < len(g.featureArray) { // check to ensure back portion of the gene was not deleted
		if g.featureArray[genomeIndexStartPos] >= 0 {
			for j = genomeIndexStartPos; g.featureArray[j] >= 0; j++ {
				g.featureArray[j] -= Feature(deletedCodingBases)
			}
		} else {
			for g.featureArray[j] < 0 {
				j++ // move to first cdsBase
				if j > len(g.featureArray)-1 {
					break
				}
			}
		}
		for _, val := range g.cdsStarts {
			if val >= j {
				for j = val; g.featureArray[j] >= 0; j++ {
					g.featureArray[j] -= Feature(deletedCodingBases)
				}
			}
		}
	}
}

// deleteUpdateUtr updates the 3' and 5' UTRs upon deletion.
// Feature array must be update prior to calling this function.
func deleteUpdateUtr(g *Gene, genomeIndexStartPos int, genomeIndexEndPos int, origCdsStart int, origCdsEnd int) (err error) {
	var utrFiveDelStart, utrFiveDelEnd int = -1, -1
	var utrThreeDelStart, utrThreeDelEnd int = -1, -1
	if genomeIndexStartPos < origCdsStart || genomeIndexEndPos > origCdsEnd {

		var utrFiveStartOffset, utrFiveEndOffset int = 0, 0
		var utrThreeStartOffset, utrThreeEndOffset int = 0, 0

		if genomeIndexStartPos < origCdsStart { // mutation overlaps 5' UTR, offset value needed
			for i := 0; g.featureArray[genomeIndexStartPos+i] < 0; i++ {
				if g.featureArray[genomeIndexStartPos+i] == UtrFive {
					utrFiveStartOffset++
					if genomeIndexStartPos+i > genomeIndexEndPos-1 {
						utrFiveEndOffset++
					}
				}
			}
		}

		if genomeIndexEndPos > origCdsEnd { // mutation overlaps 3' UTR, offset value needed
			for i := 0; g.featureArray[(genomeIndexEndPos-i)-1] < 0; i++ {
				if g.featureArray[(genomeIndexEndPos-i)-1] == UtrThree {
					utrThreeEndOffset++
					if (genomeIndexEndPos-i)-1 < genomeIndexStartPos {
						utrThreeStartOffset++
					}
				}
			}
		}

		utrFiveDelStart = len(g.utrFive.seq) - utrFiveStartOffset
		utrFiveDelEnd = len(g.utrFive.seq) - utrFiveEndOffset
		utrThreeDelStart = utrThreeStartOffset
		utrThreeDelEnd = utrThreeEndOffset
	}

	// Delete 5' UTR
	if utrFiveDelEnd != -1 {
		err = safeDelete(g, &g.utrFive.seq, utrFiveDelStart, utrFiveDelEnd, 0)
		if err != nil {
			return
		}
	}

	// Delete 3' UTR
	if utrThreeDelEnd != -1 {
		err = safeDelete(g, &g.utrThree.seq, utrThreeDelStart, utrThreeDelEnd, len(g.utrFive.seq)+len(g.codingSeq.seq))
		if err != nil {
			return
		}
	}

	return
}

// deleteEffectPrediciton performs the lionshare of the effect prediction computation upon deletion.
// Gene must be fully updated before calling this function.
func deleteEffectPrediction(g *Gene, deletedCodingBases int, answer EffectPrediction) (EffectPrediction, error) {
	var err error

	if deletedCodingBases == 0 { // Intronic or UTR (CDS unchanged)
		answer.Consequence = checkSplice(answer.CdnaDist)
		return answer, err

	} else { // Deletion affects CDS
		answer.AaPos = answer.CdnaPos / 3
		startFrame := answer.CdnaPos % 3                          // Where the start of the deletion is relative to codon boundaries
		numRemovedAA := (deletedCodingBases + startFrame + 2) / 3 // the +2 is to account for potential partial codon at del end
		removedAA := make([]dna.AminoAcid, numRemovedAA)
		copy(removedAA, g.protSeq[answer.AaPos:answer.AaPos+numRemovedAA])
		answer.AaRef = removedAA

		delFrame := deletedCodingBases % 3
		if delFrame != 0 { // Deletion causes frameshift
			answer.AaAlt = dna.TranslateSeqToTer(g.cdnaSeq[answer.CdnaPos+len(g.utrFive.seq)-delFrame:])

		} else if startFrame != 0 { // In-frame deletion does not fall on codon boundaries
			newCodonStart := answer.CdnaPos - startFrame
			answer.AaAlt = dna.TranslateSeq(g.codingSeq.seq[newCodonStart : newCodonStart+3])
		}
	}

	return answer, err
}

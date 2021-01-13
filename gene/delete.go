package gene

import (
	"errors"
	"github.com/vertgenlab/gonomics/dna"
)

var (
	ErrNoStopFound = errors.New("frameshift was unable to find a stop before the end of the gene")
)

//TODO EffectPrediction
// Deletion removes bases from the Gene, predicts the effect, and updates the Gene struct to reflect the change.
// genomeStartPos should be the base directly BEFORE the deleted bases. genomeEndPos should be the base directly AFTER the deleted bases.
// All positions should be given as base-zero genomic coordinates.
func Deletion(g *Gene, genomeStartPos int, genomeEndPos int) (EffectPrediction, error) {
	var answer EffectPrediction
	var err error

	genomeIndexStartPos, genomeIndexEndPos, err := deletionPreRunChecks(g, genomeStartPos, genomeEndPos)
	if err != nil {
		return EffectPrediction{}, err
	}

	// Save CDS start and stop for later
	origCdsStart := g.cdsStarts[0]
	origCdsEnd := g.cdsEnds[len(g.cdsEnds)-1]

	// Delete from Genome Sequence
	copy(g.genomeSeq[genomeIndexStartPos+1:], g.genomeSeq[genomeIndexEndPos:])
	g.genomeSeq = g.genomeSeq[:len(g.genomeSeq)-(genomeIndexEndPos-(genomeIndexStartPos+1))]

	// Update startPos if necessary
	if genomeIndexStartPos == -1 {
		g.startPos += genomeIndexEndPos - (genomeIndexStartPos + 1)
	}

	// Update CDS starts and ends and calculate how many cDNA bases have been deleted
	var cdsIndexesToDelete []int
	var deletedCodingBases int
	var codingDelStart, codingDelEnd int = -1, -1
	deletionLen := genomeIndexEndPos - (genomeIndexStartPos + 1)
	for i := 0; i < len(g.cdsStarts); i++ {

		switch {
		// if cds is before deletion
		case genomeIndexStartPos >= g.cdsEnds[i]:
			codingDelStart = int(g.featureArray[g.cdsEnds[i]])

		// If cds exon is fully deleted
		case genomeIndexStartPos < g.cdsStarts[i] && genomeIndexEndPos > g.cdsEnds[i]:
			// mark cds for deletion
			cdsIndexesToDelete = append(cdsIndexesToDelete, i)
			deletedCodingBases += (g.cdsEnds[i] + 1) - g.cdsStarts[i]
			if codingDelStart == -1 {
				codingDelStart = int(g.featureArray[g.cdsStarts[i]]) - 1
			}
			codingDelEnd = int(g.featureArray[g.cdsEnds[i]]) + 1

		// If deletion falls within a cds
		case genomeIndexStartPos >= g.cdsStarts[i] && genomeIndexStartPos < g.cdsEnds[i] &&
			genomeIndexEndPos > g.cdsStarts[i] && genomeIndexEndPos <= g.cdsEnds[i]:
			g.cdsEnds[i] -= deletionLen
			deletedCodingBases += deletionLen
			codingDelStart = int(g.featureArray[genomeIndexStartPos])
			codingDelEnd = int(g.featureArray[genomeIndexEndPos])

		// If cds exon is partially deleted on right
		case genomeIndexStartPos >= g.cdsStarts[i] && genomeIndexStartPos < g.cdsEnds[i]:
			deletedCodingBases += g.cdsEnds[i] - genomeIndexStartPos
			g.cdsEnds[i] = genomeIndexStartPos
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
		err = safeDelete(g, &g.codingSeq.seq, codingDelStart+1, codingDelEnd, len(g.utrFive.seq))
		if err != nil {
			return EffectPrediction{}, err
		}
	}

	// Determine UTR bases to delete
	var utrFiveDelStart, utrFiveDelEnd int = -1, -1
	var utrThreeDelStart, utrThreeDelEnd int = -1, -1
	if genomeIndexStartPos+1 < origCdsStart || genomeIndexEndPos > origCdsEnd {

		var utrFiveStartOffset, utrFiveEndOffset int = 0, 0
		var utrThreeStartOffset, utrThreeEndOffset int = 0, 0

		if genomeIndexStartPos+1 < origCdsStart { // mutation overlaps 5' UTR, offset value needed
			for i := 0; g.featureArray[genomeIndexStartPos+i+1] < 0; i++ {
				if g.featureArray[genomeIndexStartPos+i+1] == UtrFive {
					utrFiveStartOffset++
					if genomeIndexStartPos+i+1 > genomeIndexEndPos-1 {
						utrFiveEndOffset++
					}
				}
			}
		}

		if genomeIndexEndPos > origCdsEnd { // mutation overlaps 3' UTR, offset value needed
			for i := 0; g.featureArray[(genomeIndexEndPos-i)-1] < 0; i++ {
				if g.featureArray[(genomeIndexEndPos-i)-1] == UtrThree {
					utrThreeEndOffset++
					if (genomeIndexEndPos-i)-1 < genomeIndexStartPos+1 {
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
			return EffectPrediction{}, err
		}
	}

	// Delete 3' UTR
	if utrThreeDelEnd != -1 {
		err = safeDelete(g, &g.utrThree.seq, utrThreeDelStart, utrThreeDelEnd, len(g.utrFive.seq)+len(g.codingSeq.seq))
		if err != nil {
			return EffectPrediction{}, err
		}
	}

	// Update Feature Array
	copy(g.featureArray[genomeIndexStartPos+1:], g.featureArray[genomeIndexEndPos:])
	g.featureArray = g.featureArray[:len(g.featureArray)-(genomeIndexEndPos-(genomeIndexStartPos+1))]

	var j int = genomeIndexStartPos + 1
	if genomeIndexStartPos+1 < len(g.featureArray) { // check to ensure back portion of the gene was not deleted
		if g.featureArray[genomeIndexStartPos+1] >= 0 {
			for j = genomeIndexStartPos + 1; g.featureArray[j] >= 0; j++ {
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

	//TODO EffectPrediction
	g.protSeq = dna.TranslateSeq(g.codingSeq.seq) // TODO improve efficiency by updating the protein as changes are made
	return answer, err
}

// deletionPreRunChecks performs error checking and appends the diff log prior to making a deletion.
func deletionPreRunChecks(g *Gene, genomeStartPos int, genomeEndPos int) (int, int, error) {
	var log diff

	if genomeStartPos < -1 || genomeEndPos < 0 {
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

	var genomeIndexStartPos, genomeIndexEndPos int

	if g.posStrand {
		genomeIndexStartPos = genomeStartPos - g.startPos
		genomeIndexEndPos = genomeEndPos - g.startPos
	} else { // flipped for neg posStrand. start of deletion, and end of deletion must be flipped
		genomeIndexStartPos = g.startPos - genomeEndPos
		genomeIndexEndPos = g.startPos - genomeStartPos
	}

	if genomeIndexEndPos > len(g.genomeSeq) {
		genomeIndexEndPos = len(g.genomeSeq) // if deletion goes past the end of the gene, then trim the value
	}
	if genomeIndexStartPos >= len(g.genomeSeq)-1 {
		return -1, -1, ErrInputPosNotInGene
	}

	// Fill changeLog
	log.genomePos = genomeStartPos
	log.removed = make([]dna.Base, genomeIndexEndPos-(genomeIndexStartPos+1))
	copy(log.removed, g.genomeSeq[genomeIndexStartPos+1:genomeIndexEndPos])
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
	errInvalidPosition = errors.New("Input position to deleteStable is invalid")
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

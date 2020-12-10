package gene

import (
	"errors"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/numbers"
)

// PointMutation changes a single nucleotide to the desired base, predicts the effect,
// and updates the Gene struct to reflect the change.
// The position of the mutation should be given in base-zero genomic coordinates.
func PointMutation(g *Gene, genomePos int, alt dna.Base) (EffectPrediction, error) {
	var answer EffectPrediction
	var log diff
	var err error

	genomeIndexPos := numbers.AbsInt(genomePos - g.startPos) // abs needed to handle negative posStrand

	// Fill log before change
	log.genomePos = genomePos
	log.removed = make([]dna.Base, 1)
	log.removed[0] = g.genomeSeq[genomeIndexPos]
	log.added = make([]dna.Base, 1)
	log.added[0] = alt

	if !g.posStrand { // so that diff can be undone by another PointMutation call
		dna.ReverseComplement(log.removed)
	}

	if alt != dna.A && alt != dna.C && alt != dna.G && alt != dna.T {
		return answer, errors.New("alt base must be A, C, T, or G")
	}

	if genomePos < 0 {
		return answer, errors.New("genomePos must be positive")
	}
	if g.posStrand {
		if genomePos < g.startPos {
			return answer, errors.New("input genomePos is not in the gene")
		}
	} else {
		if genomePos > g.startPos {
			return answer, errors.New("input genomePos is not in the gene")
		}
		alt = dna.ComplementSingleBase(alt)
	}

	if genomeIndexPos > len(g.genomeSeq)-1 {
		return answer, errors.New("input genomePos is not in the gene")
	}

	// if no error then make the change and append log
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
		answer.CdnaPos, answer.CdnaDist, err = GenomicPosToCdna(g, genomePos)
		common.ExitIfError(err)
		answer.Consequence = checkSplice(answer.CdnaDist)
	}

	return answer, nil
}

//TODO EffectPrediction
// Insertion adds bases to the Gene, predicts the effect, and updates the Gene struct to reflect the change.
// The position should be the base-zero genomic coordinate of the base directly BEFORE the inserted bases.
func Insertion(g *Gene, genomePos int, alt []dna.Base) (EffectPrediction, error) {
	var answer EffectPrediction
	var log diff
	var genomeIndexPos int
	var err error

	// Fill changeLog
	log.genomePos = genomePos
	log.added = make([]dna.Base, len(alt))
	copy(log.added, alt)

	if !dna.IsSeqOfACTG(alt) {
		return answer, errors.New("all alt bases must be A, C, T, or G")
	}

	if genomePos < 0 {
		return answer, errors.New("genomePos must be positive")
	}
	if g.posStrand {
		if genomePos < g.startPos {
			return answer, errors.New("input genomePos is not in the gene")
		}
		genomeIndexPos = genomePos - g.startPos
	} else {
		if genomePos > g.startPos {
			return answer, errors.New("input genomePos is not in the gene")
		}
		genomeIndexPos = (g.startPos - genomePos) - 1 // -1 is so the genomeIndexPos is the base BEFORE the insertion
		dna.ReverseComplement(alt)
	}

	if genomeIndexPos > len(g.genomeSeq)-1 {
		return answer, errors.New("input genomePos is not in the gene")
	}

	// If no error, then append changeLog
	g.changeLog = append(g.changeLog, log)

	// Update genome sequence
	g.genomeSeq = dna.Insert(g.genomeSeq, int64(genomeIndexPos+1), alt)

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
		cdnaPos := int(g.featureArray[genomeIndexPos])

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

		//TODO: Effect Prediction WIP
		//answer.CdnaPos, answer.CdnaDist, err = GenomicPosToCdna(g, genomePos + 1)
		//frame := (cdnaPos + 1) % 3
		//var currCodon dna.Codon
		//
		//if frame != 0 { // if insertion disrupts a codon, note the codon being changed
		//	currCodon, err = CdnaPosToCodon(g, cdnaPos)
		//	answer.AaRef = []dna.AminoAcid{dna.TranslateCodon(&currCodon)}
		//}

		// Update cDNA
		g.cdnaSeq = dna.Insert(g.cdnaSeq, int64(cdnaPos+1), alt)

		//TODO: Effect Prediction WIP
		//answer.AaPos = cdnaPos / 3
		//
		//if len(alt) % 3 != 0 { // Causes Frameshift // TODO this is a bit jury rigged, but it works....
		//	answer.Consequence = Frameshift
		//	currCodon, err = CdnaPosToCodon(g, cdnaPos + 1)
		//	var stopFound bool = true
		//	for currCdnaPos := cdnaPos + 1 + 3; dna.TranslateCodon(&currCodon) != dna.Stop; currCdnaPos += 3 {
		//		answer.AaAlt = append(answer.AaAlt, dna.TranslateCodon(&currCodon))
		//		if currCdnaPos + (3 - frame) > len(g.cdnaSeq) { // if do dont hit a stop by the end of the cdnaSeq, then move to genomic DNA
		//			currGdnaPos := g.cdsEnds[len(g.cdsEnds)-1] + 1 - frame
		//			for ;dna.TranslateSeq(g.genomeSeq[currGdnaPos:currGdnaPos + 3])[0] != dna.Stop; currGdnaPos += 3 {
		//				answer.AaAlt = append(answer.AaAlt, dna.TranslateSeq(g.genomeSeq[currGdnaPos:currGdnaPos + 3])[0])
		//				if currGdnaPos + 3 > len(g.genomeSeq) {
		//					stopFound = false
		//					err = errors.New("frameshift was unable to find a stop before the end of the gene")
		//					break
		//				}
		//			}
		//			break
		//		} else {
		//			currCodon, err = CdnaPosToCodon(g, currCdnaPos)
		//		}
		//	}
		//	if stopFound {
		//		answer.AaAlt = append(answer.AaAlt, dna.Stop)
		//	}
		//} else { // In-Frame
		//	answer.Consequence = InFrameInsertion
		//	if frame != 0 { // Disrupts an existing codon
		//		answer.AaAlt = dna.TranslateSeq(g.cdnaSeq[(cdnaPos + 1) - frame : (cdnaPos + 1) + len(alt) + (3 - frame)])
		//		if answer.AaRef[0] == answer.AaAlt[0] {
		//			answer.AaRef = nil
		//			answer.AaAlt = answer.AaAlt[1:]
		//			answer.AaPos++
		//		}
		//	} else { // Does not disrupt existing codons AND is inframe
		//		answer.AaAlt = dna.TranslateSeq(alt)
		//	}
		//}

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
		//TODO: Effect Prediction WIP
		//var endCdnaOffset int
		//answer.CdnaPos, answer.CdnaDist, err = GenomicPosToCdna(g, genomePos + 1)
		//_, endCdnaOffset, err = GenomicPosToCdna(g, genomePos + 1 + (len(alt) - 1))
		//if numbers.AbsInt(endCdnaOffset) < numbers.AbsInt(answer.CdnaDist) {
		//	answer.Consequence = checkSplice(endCdnaOffset)
		//} else {
		//	answer.Consequence = checkSplice(answer.CdnaDist)
		//}
	}
	return answer, err
}

//TODO EffectPrediction
// Deletion removes bases from the Gene, predicts the effect, and updates the Gene struct to reflect the change.
// genomeStartPos should be the base directly BEFORE the deleted bases. genomeEndPos should be the base directly AFTER the deleted bases.
// All positions should be given as base-zero genomic coordinates.
func Deletion(g *Gene, genomeStartPos int, genomeEndPos int) (EffectPrediction, error) {
	var answer EffectPrediction
	var log diff
	var err error

	if genomeStartPos < 0 || genomeEndPos < 0 {
		return answer, errors.New("genome positions must be positive")
	}

	if genomeStartPos >= genomeEndPos {
		return answer, errors.New("genomeStartPos must be less than genomeEndPos")
	}

	if g.posStrand {
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
				cdnaDelEnd = int(g.featureArray[g.cdsStarts[i]])
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

	var j int = genomeIndexStartPos + 1
	if g.featureArray[genomeIndexStartPos] >= 0 {
		for j = genomeIndexStartPos + 1; g.featureArray[j] >= 0; j++ {
			g.featureArray[j] -= Feature(deletedCdnaBases)
		}
	} else {
		for g.featureArray[j] < 0 {
			j++ // move to first cdsBase
		}
	}
	for _, val := range g.cdsStarts {
		if val >= j {
			for j = val; g.featureArray[j] >= 0; j++ {
				g.featureArray[j] -= Feature(deletedCdnaBases)
			}
		}
	}

	//TODO EffectPrediction

	return answer, err
}

// Reset reverts all mutations done to a Gene.
func Reset(g *Gene) {
	var hasDel bool
	var err error
	for _, val := range g.changeLog {
		if len(val.added) == 0 && len(val.removed) > 0 {
			hasDel = true
			break
		}
	}

	if !hasDel { // Point mutations and insertions can be quickly reversed by performing the inverse function
		for i := len(g.changeLog) - 1; i >= 0; i-- {
			if len(g.changeLog[i].added) == 1 && len(g.changeLog[i].removed) == 1 {
				_, err = PointMutation(g, g.changeLog[i].genomePos, g.changeLog[i].removed[0])
				g.changeLog = g.changeLog[:len(g.changeLog)-2]
				if err != nil {
					hasDel = true
					break
				}
			} else if len(g.changeLog[i].added) >= 1 && len(g.changeLog[i].removed) == 0 {
				_, err = Deletion(g, g.changeLog[i].genomePos, g.changeLog[i].genomePos+len(g.changeLog[i].added)+1)
				g.changeLog = g.changeLog[:len(g.changeLog)-2]
				if err != nil {
					hasDel = true
					break
				}
			}
		}
	}

	if hasDel { // Deletions can remove important feature context which is non-trivial to store in the changelog, so insertions must restore from backup
		g.startPos = g.orig.startPos
		g.cdsStarts = g.cdsStarts[:len(g.orig.cdsStarts)]
		copy(g.cdsStarts, g.orig.cdsStarts)
		g.cdsEnds = g.cdsEnds[:len(g.orig.cdsEnds)]
		copy(g.cdsEnds, g.orig.cdsEnds)
		g.genomeSeq = g.genomeSeq[:len(g.orig.genomeSeq)]
		copy(g.genomeSeq, g.orig.genomeSeq)
		g.cdnaSeq = g.cdnaSeq[:len(g.orig.cdnaSeq)]
		copy(g.cdnaSeq, g.orig.cdnaSeq)
		g.featureArray = g.featureArray[:len(g.orig.featureArray)]
		copy(g.featureArray, g.orig.featureArray)
	}

	g.changeLog = nil
}

// checkSplice inputs a cDNA offset value and determines if a variant at this site may affect splicing.
func checkSplice(CdnaOffset int) MutationType {
	if numbers.AbsInt(CdnaOffset) <= 2 {
		return Splice
	}

	if numbers.AbsInt(CdnaOffset) <= 10 {
		return FarSplice
	}

	return Intronic
}

package gene

import (
	"errors"
	"log"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/numbers"
)

// Insertion adds bases to the Gene, predicts the effect, and updates the Gene struct to reflect the change.
// The position should be the base-zero genomic coordinate of the base directly BEFORE the inserted bases.
func Insertion(g *Gene, genomePos int, alt []dna.Base) (EffectPrediction, error) {
	var answer EffectPrediction
	answer.StopDist = -1
	var genomeIndexPos int
	var err error

	genomeIndexPos, err = insertionPreRunChecks(g, genomePos, alt)

	// Update genome sequence
	g.genomeSeq = dna.Insert(g.genomeSeq, genomeIndexPos+1, alt)

	// Update CDS starts and ends
	for idx, val := range g.cdsStarts {
		if val > genomeIndexPos {
			g.cdsStarts[idx] += len(alt)
			g.cdsEnds[idx] += len(alt)
		} else if g.cdsEnds[idx] > genomeIndexPos {
			g.cdsEnds[idx] += len(alt)
		}
	}

	if g.featureArray[genomeIndexPos] >= 0 && g.featureArray[genomeIndexPos+1] >= 0 { // Coding
		codingPos := int(g.featureArray[genomeIndexPos])

		// Update featureArray for coding
		fillVal := Feature(codingPos + 1)
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

		answer.CdnaPos, answer.CdnaDist, err = GenomicPosToCdna(g, genomePos)
		frame := (codingPos + 1) % 3
		var currCodon dna.Codon

		if frame != 0 { // if insertion disrupts a codon, note the codon being changed
			currCodon, err = CdnaPosToCodon(g, codingPos)
			if err != nil {
				log.Panic("error with cDNA pos")
			}
			answer.AaRef = []dna.AminoAcid{dna.TranslateCodon(currCodon)}
		}

		// Update cDNA
		err = insertAdjust(g, &g.cdnaSeq, codingPos+1+len(g.utrFive.seq), alt)
		if err != nil {
			return EffectPrediction{}, err
		}

		var newProtSeq []dna.AminoAcid

		answer.AaPos = codingPos / 3

		if len(alt)%3 != 0 { // Causes Frameshift
			answer.StopDist = 0
			answer.Consequence = Frameshift

			newProtSeq = frameshiftTranslate(g.codingSeq.seq[(codingPos+1)-frame:], g.utrThree.seq)

			if newProtSeq[len(newProtSeq)-1] == dna.Stop {
				answer.StopDist = len(newProtSeq) - 1
			} else {
				answer.StopDist = -2
			}

			var j int
			for j = 0; newProtSeq[j] == g.protSeq[answer.AaPos]; j++ { // trim matching amino acids from frameshift
				answer.AaPos++
				if answer.StopDist != -2 {
					answer.StopDist--
				}
				if answer.AaPos >= len(g.protSeq) {
					break
				}
			}

			answer.AaRef = []dna.AminoAcid{g.protSeq[answer.AaPos]}
			answer.AaAlt = []dna.AminoAcid{newProtSeq[j]}

			g.protSeq = newProtSeq //TODO append new prot seq to existing seq at point of insertion
		} else { // In-Frame
			newProtSeq = dna.TranslateSeqToTer(g.codingSeq.seq) //TODO this can be much more efficient
			answer.Consequence = InFrameInsertion
			if frame != 0 { // Disrupts an existing codon
				answer.AaAlt = dna.TranslateSeq(g.codingSeq.seq[(codingPos+1)-frame : (codingPos+1)+len(alt)+(3-frame)])
				if answer.AaRef[0] == answer.AaAlt[0] {
					answer.AaRef = nil
					answer.AaAlt = answer.AaAlt[1:]
					answer.AaPos++
				}
			} else { // Does not disrupt existing codons AND is inframe
				answer.AaAlt = dna.TranslateSeq(alt)
			}
		}
		g.protSeq = newProtSeq
	} else { // NonCoding
		// Update featureArray for noncoding
		// figure out how the inserted bases should be labeled in the feature array. Priority: 5'UTR=3'Utr>Intron
		fillVal := Feature(numbers.Min(int(g.featureArray[genomeIndexPos]), int(g.featureArray[genomeIndexPos+1])))

		switch fillVal {
		case UtrFive:
			utrFiveOffset := 0

			for i := 0; i <= genomeIndexPos; i++ {
				if g.featureArray[i] == UtrFive {
					utrFiveOffset++
				}
			}

			if utrFiveOffset > -1 { // Update 5' UTR
				err = insertAdjust(g, &g.cdnaSeq, g.utrFive.start+utrFiveOffset, alt)
				if err != nil {
					return EffectPrediction{}, err
				}
			}

		case UtrThree:
			utrThreeOffset := 0

			for i := genomeIndexPos; g.featureArray[i] < 0; i-- {
				if g.featureArray[i] == UtrThree {
					utrThreeOffset++
				}
				if i <= 0 { // case where entire gene is 3' UTR
					break
				}
			}

			if utrThreeOffset > -1 { // Update 3' UTR
				err = insertAdjust(g, &g.cdnaSeq, g.utrThree.start+utrThreeOffset, alt)
				if err != nil {
					return EffectPrediction{}, err
				}
			}
		}

		tmpSlice := make([]Feature, len(alt))
		g.featureArray = append(g.featureArray, tmpSlice...)
		copy(g.featureArray[genomeIndexPos+1+len(alt):], g.featureArray[genomeIndexPos+1:])
		for i := 0; i < len(alt); i++ {
			g.featureArray[genomeIndexPos+1+i] = fillVal
		}

		var endCdnaOffset int
		answer.CdnaPos, answer.CdnaDist, err = GenomicPosToCdna(g, genomePos+1)
		_, endCdnaOffset, err = GenomicPosToCdna(g, genomePos+1+(len(alt)-1))
		if err != nil {
			log.Panic("error with cDNA pos")
		}
		if numbers.AbsInt(endCdnaOffset) < numbers.AbsInt(answer.CdnaDist) {
			answer.Consequence = checkSplice(endCdnaOffset)
		} else {
			answer.Consequence = checkSplice(answer.CdnaDist)
		}
	}
	return answer, err
}

var (
	ErrExceedCapNewSliceCreated = errors.New("Capacity of input slice was exceeded, output slice was created at a new memory address")
)

// insertionPreRunChecks performs error checking and appends the diff log prior to making an insertion.
func insertionPreRunChecks(g *Gene, genomePos int, alt []dna.Base) (int, error) {
	var log diff
	var genomeIndexPos int

	// Fill changeLog
	log.genomePos = genomePos
	log.added = make([]dna.Base, len(alt))
	copy(log.added, alt)

	if !dna.IsSeqOfACGT(alt) {
		return -1, errors.New("nonstandard base")
	}

	if genomePos < 0 {
		return -1, ErrNegativeInputValue
	}
	if g.posStrand {
		if genomePos < g.startPos {
			return -1, ErrInputPosNotInGene
		}
		genomeIndexPos = genomePos - g.startPos
	} else {
		if genomePos > g.startPos {
			return -1, ErrInputPosNotInGene
		}
		genomeIndexPos = (g.startPos - genomePos) - 1 // -1 is so the genomeIndexPos is the base BEFORE the insertion
		dna.ReverseComplement(alt)
	}

	if genomeIndexPos > len(g.genomeSeq)-1 {
		return -1, ErrInputPosNotInGene
	}

	// If no error, then append changeLog
	g.changeLog = append(g.changeLog, log)

	return genomeIndexPos, nil
}

// insertAdjust is a wrapper for insertStable that adjusts the start and end positions of UTRs and CDS.
func insertAdjust(g *Gene, destSeq *[]dna.Base, insPos int, insSeq []dna.Base) error {
	var err error
	err = insertStable(destSeq, insPos, insSeq)
	if errors.Is(err, ErrExceedCapNewSliceCreated) {
		err = nil
	}

	if g.utrFive.start > insPos {
		g.utrFive.start += len(insSeq)
	}

	if g.utrFive.end > insPos {
		g.utrFive.end += len(insSeq)
	}

	if g.utrThree.start > insPos {
		g.utrThree.start += len(insSeq)
	}

	if g.utrThree.end > insPos {
		g.utrThree.end += len(insSeq)
	}

	if g.codingSeq.start > insPos {
		g.codingSeq.start += len(insSeq)
	}

	if g.codingSeq.end > insPos {
		g.codingSeq.end += len(insSeq)
	}

	g.utrFive.seq = g.cdnaSeq[g.utrFive.start:g.utrFive.end]
	g.utrThree.seq = g.cdnaSeq[g.utrThree.start:g.utrThree.end]
	g.codingSeq.seq = g.cdnaSeq[g.codingSeq.start:g.codingSeq.end]

	return nil
}

// insertStable is identical to Insert, but conserves the
// underlying memory of the input destSeq slice if possible.
func insertStable(destSeq *[]dna.Base, insPos int, insSeq []dna.Base) error {
	if insPos < 0 || insPos > len(*destSeq) {
		return errInvalidPosition
	}
	if len(*destSeq)+len(insSeq) > cap(*destSeq) {
		*destSeq = dna.Insert(*destSeq, insPos, insSeq)
		return ErrExceedCapNewSliceCreated
	}
	*destSeq = (*destSeq)[:len(*destSeq)+len(insSeq)]                                    // extend the LEN of the slice to include the inserted bases
	copy((*destSeq)[insPos+len(insSeq):cap(*destSeq)], (*destSeq)[insPos:cap(*destSeq)]) // shift bases in slice by len(insSeq) at the position of the insertion
	copy((*destSeq)[insPos:insPos+len(insSeq)], insSeq)                                  // add in the inserted bases
	return nil
}

// frameshiftTranslate translates the result of a frameshift going into the 3'UTR sequence if necessary.
// Returns a stop-terminated sequence of amino acids.
func frameshiftTranslate(shiftedCds []dna.Base, utrThree []dna.Base) []dna.AminoAcid {
	answer := make([]dna.AminoAcid, 0, (len(shiftedCds)+len(utrThree))/3)

	frameOffset := len(shiftedCds) % 3

	for i := 0; i < len(shiftedCds)-frameOffset; i += 3 {
		if i+3 > len(shiftedCds) {
			break
		}
		answer = append(answer, dna.TranslateSeq(shiftedCds[i : i+3])[0])
		if answer[len(answer)-1] == dna.Stop {
			return answer
		}
	}

	transitionCodon := append(shiftedCds[len(shiftedCds)-frameOffset:], utrThree[:3-frameOffset]...)
	answer = append(answer, dna.TranslateSeq(transitionCodon)[0])
	if answer[len(answer)-1] == dna.Stop {
		return answer
	}

	remainingUtrThree := utrThree[3-frameOffset:]

	for i := 0; i < len(remainingUtrThree) && i+3 <= len(remainingUtrThree); i += 3 {
		answer = append(answer, dna.TranslateSeq(remainingUtrThree[i : i+3])[0])
		if answer[len(answer)-1] == dna.Stop {
			return answer
		}
	}
	return answer
}

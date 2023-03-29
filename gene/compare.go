package gene

import (
	"errors"
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
)

func equal(alpha, beta *Gene) (bool, error) {
	if alpha.id != beta.id {
		return false, errors.New("id")
	}

	if alpha.startPos != beta.startPos {
		return false, errors.New("startPos")
	}

	if alpha.posStrand != beta.posStrand {
		return false, errors.New("posStrand")
	}

	if len(alpha.cdsStarts) != len(beta.cdsStarts) {
		return false, errors.New("cdsStarts")
	}

	for idx, val := range alpha.cdsStarts {
		if val != beta.cdsStarts[idx] {
			return false, errors.New("cdsStarts")
		}
	}

	if len(alpha.cdsEnds) != len(beta.cdsEnds) {
		return false, errors.New("cdsEnds")
	}

	for idx, val := range alpha.cdsEnds {
		if val != beta.cdsEnds[idx] {
			return false, errors.New("cdsEnds")
		}
	}

	if len(alpha.genomeSeq) != len(beta.genomeSeq) {
		return false, errors.New("genomeSeq")
	}

	for idx, val := range alpha.genomeSeq {
		if val != beta.genomeSeq[idx] {
			return false, errors.New("genomeSeq")
		}
	}

	if len(alpha.cdnaSeq) != len(beta.cdnaSeq) {
		return false, errors.New("cdnaSeq")
	}

	for idx, val := range alpha.cdnaSeq {
		if val != beta.cdnaSeq[idx] {
			return false, errors.New("cdnaSeq")
		}
	}

	if len(alpha.featureArray) != len(beta.featureArray) {
		return false, errors.New("featureArray")
	}

	for idx, val := range alpha.featureArray {
		if val != beta.featureArray[idx] {
			return false, errors.New("featureArray")
		}
	}

	if !equalSubSeq(alpha.utrFive, beta.utrFive) {
		return false, errors.New("utrFive")
	}

	if !equalSubSeq(alpha.utrThree, beta.utrThree) {
		return false, errors.New("utrThree")
	}

	if !equalSubSeq(alpha.codingSeq, beta.codingSeq) {
		return false, errors.New("codingSeq")
	}

	return true, nil
}

func equalSubSeq(alpha subSeq, beta subSeq) bool {
	if alpha.start != beta.start || alpha.end != beta.end {
		return false
	}
	if len(alpha.seq) != len(beta.seq) {
		return false
	}
	for idx, val := range alpha.seq {
		if val != beta.seq[idx] {
			return false
		}
	}
	return true
}

// for manual testing
func printEffPred(pred EffectPrediction) {
	var consequence string
	switch pred.Consequence {
	case Intronic:
		consequence = "Intronic"
	case Silent:
		consequence = "Silent"
	case Missense:
		consequence = "Missense"
	case Nonsense:
		consequence = "Nonsense"
	case Frameshift:
		consequence = "Frameshift"
	case Intergenic:
		consequence = "Intergenic"
	case Splice:
		consequence = "Splice"
	case FarSplice:
		consequence = "FarSplice"
	case DisruptStart:
		consequence = "DisruptStart"
	case DisruptStop:
		consequence = "DisruptStop"
	case InFrameInsertion:
		consequence = "InFrameInsertion"
	case InFrameDeletion:
		consequence = "InFrameDeletion"
	}
	fmt.Printf("Consequence: %s\n", consequence)
	fmt.Printf("cDNA Pos: %d%+d\n", pred.CdnaPos, pred.CdnaDist)
	fmt.Printf("AAPos: %d\n", pred.AaPos)
	fmt.Printf("AaRef: %s\n", dna.PeptideToString(pred.AaRef))
	fmt.Printf("AaAlt: %s\n", dna.PeptideToString(pred.AaAlt))
	fmt.Printf("StopDist: %d\n", pred.StopDist)
}

func equalPred(alpha, beta *EffectPrediction) (bool, error) {
	if alpha.Consequence != beta.Consequence {
		return false, errors.New("Consequence")
	}

	if alpha.CdnaPos != beta.CdnaPos {
		return false, errors.New("CdnaPos")
	}

	if alpha.AaPos != beta.AaPos {
		return false, errors.New("AaPos")
	}

	if alpha.StopDist != beta.StopDist {
		return false, errors.New("StopDist")
	}

	if alpha.CdnaDist != beta.CdnaDist {
		return false, errors.New("CdnaDist")
	}

	if len(alpha.AaRef) != len(beta.AaRef) {
		return false, errors.New("AaRef")
	}

	for idx, val := range alpha.AaRef {
		if val != beta.AaRef[idx] {
			return false, errors.New("cdsAaRefStarts")
		}
	}

	if len(alpha.AaAlt) != len(beta.AaAlt) {
		return false, errors.New("AaAlt")
	}

	for idx, val := range alpha.AaAlt {
		if val != beta.AaAlt[idx] {
			return false, errors.New("AaAlt")
		}
	}
	return true, nil
}

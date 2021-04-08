package variant

import (
	"github.com/vertgenlab/gonomics/dna"
)

// CodingChange describes a change of a protein sequence caused by
// a mutation affecting the coding sequence of the corresponding gene.
// CodingPos and ProteinPos refer to the zero-based start position of
// the change. The CodingPos always refers to the position changed in
// origin DNA sequence, however ProteinPos is not always equal to
// CodingPos / 3. If a change would result in a RemovedAa[0] == AddedAa[0]
// the first amino acid from each slice is removed and ProteinPos += 1.
//
// RemovedAa records all amino acids that are removed from the protein
// sequence. AddedAa records all amino acids that are added to the
// protein sequence. Substitution of one amino acid for another is
// reported as the original amino acid being removed, and the new amino
// acid being added. dna.Stop is included in RemovedAa and AddedAa.
type CodingChange struct {
	CodingPos  int
	ProteinPos int
	RemovedAa  []dna.AminoAcid
	AddedAa    []dna.AminoAcid
	Type       EffectType
}

// EffectType describes the overall effect of a CodingChange.
// Possible EffectTypes are described below.
type EffectType byte

const (
	Silent           EffectType = 6 // mutation changes codon, but still codes for the same amino acid
	Frameshift       EffectType = 5 // shifts reading frame
	Nonsense         EffectType = 4 // premature termination
	InFrameInsertion EffectType = 3 // insertion where insertedBases % 3 == 0
	InFrameDeletion  EffectType = 2 // deletion where deletedBases % 3 == 0
	Missense         EffectType = 1 // mutation changes coded amino acid
)

// Effector is satisfied by types that implement the Effect method, which returns
// a CodingChange describing how changes prescribed by the receiver type alter the
// resulting amino acid sequence corresponding to the input []dna.Base. The input
// sequence should be the CDS of a gene such that input[0:3] is the start codon.
// The last 3 dna.Base of the input is not required to be a stop codon. Trailing
// sequence after the stop codon (i.e. the 3' UTR) may be optionally added to the
// input sequence, however the deletion position must be within the coding sequence.
// This trailing sequence is used when a variant causes a frameshift and a new stop
// codon in not found prior to the original stop codon. In this case the translation
// will continue into the 3' UTR attempting to find a new stop codon.
//
// The offset integer input for the Effect method is added to the position in
// the receiver type and may be positive or negative. Often this offset will
// be the negative genomic position of the beginning of the start codon.
//
// Effector is satisfied by Substitution, Insertion, Deletion, and Delins.
// The effects of Structural variation are often too complex to be described by
// CodingChange and therefore should be handled separately.
type Effector interface {
	Effect(codingSeq []dna.Base, offsetStart int, offsetEnd int) (CodingChange, error)
}

// Effect for Substitution determines a coding change resulting from a substitution.
func (s Substitution) Effect(codingSeq []dna.Base, offsetStart int, offsetEnd int) (CodingChange, error) {
	var answer CodingChange
	offsetPos := s.Pos + offsetStart
	answer.CodingPos = offsetPos
	answer.ProteinPos = offsetPos / 3

	if offsetPos < 0 {
		return answer, ErrNegPos
	}

	if codingSeq[offsetPos] != s.Ref {
		return answer, ErrRefMatch
	}

	frame := offsetPos % 3
	codonStart := offsetPos - frame

	codon := dna.Codon{codingSeq[codonStart], codingSeq[codonStart+1], codingSeq[codonStart+2]}
	refAa := dna.GeneticCode[codon]

	codon[frame] = s.Alt // mutate codon
	altAa := dna.GeneticCode[codon]

	if refAa != altAa {
		answer.RemovedAa = []dna.AminoAcid{refAa}
		answer.AddedAa = []dna.AminoAcid{altAa}
	}

	switch {
	case altAa == refAa:
		answer.Type = Silent

	case altAa == dna.Stop:
		answer.Type = Nonsense

	default:
		answer.Type = Missense
	}

	return answer, nil
}

// Effect for Insertion determines a coding change resulting from an insertion.
func (i Insertion) Effect(codingSeq []dna.Base, offsetStart int, offsetEnd int) (CodingChange, error) {
	var answer CodingChange
	offsetPos := i.Pos + offsetStart
	answer.CodingPos = offsetPos
	answer.ProteinPos = offsetPos / 3
	var protOffset int

	if offsetPos < 0 {
		return answer, ErrNegPos
	}

	if offsetPos > len(codingSeq) {
		return answer, ErrInvalidPosition
	}

	frame := offsetPos % 3
	codonStart := offsetPos - frame

	switch {
	case len(i.Seq)%3 != 0: // frameshift
		answer.Type = Frameshift
		shiftedSeq := make([]dna.Base, len(i.Seq)+len(codingSeq[codonStart:]))
		copy(shiftedSeq[:frame], codingSeq[codonStart:offsetPos])  // add bases before insertion
		copy(shiftedSeq[frame:frame+len(i.Seq)], i.Seq)            // add inserted bases
		copy(shiftedSeq[frame+len(i.Seq):], codingSeq[offsetPos:]) // add bases after insertion
		answer.RemovedAa, answer.AddedAa, protOffset = aaChange(codingSeq[codonStart:], shiftedSeq)

	case frame != 0: // disrupts a codon
		answer.Type = InFrameInsertion
		insSeq := make([]dna.Base, len(i.Seq)+3)
		copy(insSeq[:frame], codingSeq[codonStart:offsetPos])              // add bases before insertion
		copy(insSeq[frame:frame+len(i.Seq)], i.Seq)                        // add inserted sequence
		copy(insSeq[frame+len(i.Seq):], codingSeq[offsetPos:codonStart+3]) // add bases after insertion
		answer.RemovedAa, answer.AddedAa, protOffset = aaChange(codingSeq[codonStart:codonStart+3], insSeq)

	default: // in frame & does not disrupt a codon
		answer.Type = InFrameInsertion
		answer.RemovedAa, answer.AddedAa, protOffset = aaChange(nil, i.Seq)
	}

	answer.ProteinPos += protOffset
	return answer, nil
}

// Effect for Deletion determines a coding change resulting from a deletion.
func (d Deletion) Effect(codingSeq []dna.Base, offsetStart int, offsetEnd int) (CodingChange, error) {
	var answer CodingChange
	offsetStartPos := d.Start + offsetStart
	offsetEndPos := d.End + offsetEnd
	answer.CodingPos = offsetStartPos
	answer.ProteinPos = offsetStartPos / 3
	var protOffset int

	if offsetStartPos < 0 {
		return answer, ErrNegPos
	}

	if offsetEndPos > len(codingSeq) {
		return answer, ErrInvalidPosition
	}

	delLen := offsetEndPos - offsetStartPos
	startFrame := offsetStartPos % 3
	endFrame := (offsetEndPos - 1) % 3
	codonStart := offsetStartPos - startFrame       // start of first codon affected (closed)
	codonEnd := ((offsetEndPos - 1) - endFrame) + 3 // end of last codon affected (open)

	switch {
	case delLen%3 != 0: // frameshift
		answer.Type = Frameshift
		shiftedSeq := make([]dna.Base, len(codingSeq[codonStart:])-delLen)
		copy(shiftedSeq[:startFrame], codingSeq[codonStart:offsetStartPos]) // add bases before deletion
		copy(shiftedSeq[startFrame:], codingSeq[offsetEndPos:])             // add bases after deletion
		answer.RemovedAa, answer.AddedAa, protOffset = aaChange(codingSeq[codonStart:], shiftedSeq)

	case offsetStartPos%3 != 0: // disrupts a codon
		answer.Type = InFrameDeletion
		newCodon := make([]dna.Base, 3)
		copy(newCodon[:startFrame], codingSeq[codonStart:offsetStartPos]) // add bases before deletion
		copy(newCodon[startFrame:], codingSeq[offsetEndPos:codonEnd])     // add bases after deletion
		answer.RemovedAa, answer.AddedAa, protOffset = aaChange(codingSeq[codonStart:codonEnd], newCodon)

	default: // in frame & does not disrupt a codon
		answer.Type = InFrameDeletion
		answer.RemovedAa, answer.AddedAa, protOffset = aaChange(codingSeq[codonStart:codonEnd], nil)
	}

	answer.ProteinPos += protOffset
	return answer, nil
}

// Effect for Delins determines a coding change resulting from a combined deletion and insertion.
func (di Delins) Effect(codingSeq []dna.Base, offsetStart int, offsetEnd int) (CodingChange, error) {
	var answer CodingChange
	offsetStartPos := di.Start + offsetStart
	offsetEndPos := di.End + offsetEnd
	answer.CodingPos = offsetStartPos
	answer.ProteinPos = offsetStartPos / 3
	var protOffset int

	if offsetStartPos < 0 {
		return answer, ErrNegPos
	}

	if offsetEndPos > len(codingSeq) {
		return answer, ErrInvalidPosition
	}

	delLen := offsetEndPos - offsetStartPos
	lenDiff := len(di.InsSeq) - delLen
	startFrame := offsetStartPos % 3
	endFrame := (offsetEndPos - 1) % 3
	codonStart := offsetStartPos - startFrame       // start of first codon affected (closed)
	codonEnd := ((offsetEndPos - 1) - endFrame) + 3 // end of last codon affected (open)
	switch {
	case lenDiff%3 != 0: // frameshift
		answer.Type = Frameshift
		shiftedSeq := make([]dna.Base, len(codingSeq[codonStart:])+lenDiff)
		copy(shiftedSeq[:startFrame], codingSeq[codonStart:offsetStartPos])    // add bases before delins
		copy(shiftedSeq[startFrame:startFrame+len(di.InsSeq)], di.InsSeq)      // add inserted sequence
		copy(shiftedSeq[startFrame+len(di.InsSeq):], codingSeq[offsetEndPos:]) // add bases after delins
		answer.RemovedAa, answer.AddedAa, protOffset = aaChange(codingSeq[codonStart:], shiftedSeq)

	case offsetStartPos%3 != 0: // disrupts a codon
		if lenDiff > 0 {
			answer.Type = InFrameInsertion
		} else {
			answer.Type = InFrameDeletion
		}
		insSeq := make([]dna.Base, len(codingSeq[codonStart:codonEnd])+lenDiff)
		copy(insSeq[:startFrame], codingSeq[codonStart:offsetStartPos])            // add bases before insertion
		copy(insSeq[startFrame:startFrame+len(di.InsSeq)], di.InsSeq)              // add inserted sequence
		copy(insSeq[startFrame+len(di.InsSeq):], codingSeq[offsetEndPos:codonEnd]) // add bases after insertion
		answer.RemovedAa, answer.AddedAa, protOffset = aaChange(codingSeq[codonStart:codonEnd], insSeq)

	default: // in frame & does not disrupt a codon
		if lenDiff > 0 {
			answer.Type = InFrameInsertion
		} else if lenDiff < 0 {
			answer.Type = InFrameDeletion
		} else {
			answer.Type = Missense
		}
		answer.RemovedAa, answer.AddedAa, protOffset = aaChange(codingSeq[codonStart:codonEnd], di.InsSeq)

		if len(answer.AddedAa) == 0 && len(answer.RemovedAa) == 0 {
			answer.Type = Silent
		}
	}

	answer.ProteinPos += protOffset
	return answer, nil
}

// aaChange determines the amino acid change of the input sequence. Input sequences must be in-frame
// beginning with the start of the first codon affected, and ending with the end of the last codon affected.
// protOffset should be added to the ProteinPos field.
func aaChange(ref []dna.Base, alt []dna.Base) (removed []dna.AminoAcid, added []dna.AminoAcid, protOffset int) {
	removed = dna.TranslateSeqToTer(ref)
	added = dna.TranslateSeqToTer(alt)

	// trim matching amino acids
	for len(removed) > 0 && len(added) > 0 && removed[0] == added[0] {
		removed = removed[1:]
		added = added[1:]
		protOffset++
	}

	return
}

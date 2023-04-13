package gene

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/numbers"
)

// Substitution changes a single nucleotide to the desired base, predicts the effect,
// and updates the Gene struct to reflect the change.
// The position of the mutation should be given in base-zero genomic coordinates.
func Substitution(g *Gene, genomePos int, alt dna.Base) (EffectPrediction, error) {
	var answer EffectPrediction
	answer.StopDist = -1
	var err error

	err = substitutionPreRunChecks(g, genomePos, &alt)
	if err != nil {
		return EffectPrediction{}, err
	}

	genomeIndexPos := numbers.AbsInt(genomePos - g.startPos) // abs needed to handle negative posStrand
	g.genomeSeq[genomeIndexPos] = alt

	cdnaIndexPos := int(g.featureArray[genomeIndexPos])

	if cdnaIndexPos >= 0 { // mutation is coding
		err = substitutionCoding(g, cdnaIndexPos, alt, &answer)
	} else { // mutation is noncoding
		err = substitutionNonCoding(g, genomePos, &answer)
	}

	if err != nil {
		return EffectPrediction{}, err
	}

	g.protSeq = dna.TranslateSeq(g.codingSeq.seq) // TODO improve efficiency by updating the protein as changes are made
	return answer, nil
}

// checkSplice inputs the distance to the nearest E-I boundary for a given INTRONIC variant
// and determines if a variant at this site may affect splicing.
func checkSplice(distToIntronExonBoundary int) MutationType {
	if numbers.AbsInt(distToIntronExonBoundary) <= 2 {
		return Splice
	}

	if numbers.AbsInt(distToIntronExonBoundary) <= 10 {
		return FarSplice
	}

	return Intronic
}

// substitutionPreRunChecks performs error checking and appends the diff log prior to making a point mutation.
func substitutionPreRunChecks(g *Gene, genomePos int, alt *dna.Base) error {
	var log diff

	genomeIndexPos := numbers.AbsInt(genomePos - g.startPos) // abs needed to handle negative posStrand

	// Fill log before change
	log.genomePos = genomePos
	log.removed = make([]dna.Base, 1)
	log.removed[0] = g.genomeSeq[genomeIndexPos]
	log.added = make([]dna.Base, 1)
	log.added[0] = *alt

	if !g.posStrand { // so that diff can be undone by another Substitution call
		dna.ReverseComplement(log.removed)
	}

	if !dna.IsSeqOfACGT([]dna.Base{*alt}) {

		return ErrNonACGTBase
	}

	if genomePos < 0 {
		return ErrNegativeInputValue
	}
	if g.posStrand {
		if genomePos < g.startPos {
			return ErrInputPosNotInGene
		}
	} else {
		if genomePos > g.startPos {
			return ErrInputPosNotInGene
		}
		*alt = dna.ComplementSingleBase(*alt)
	}

	if genomeIndexPos > len(g.genomeSeq)-1 {
		return ErrInputPosNotInGene
	}

	// if no error then  append log
	g.changeLog = append(g.changeLog, log)
	return nil
}

// substitutionNonCoding determines the effects of a non-coding point mutation.
func substitutionNonCoding(g *Gene, genomePos int, answer *EffectPrediction) error {
	var err error
	answer.CdnaPos, answer.CdnaDist, err = GenomicPosToCdna(g, genomePos)
	answer.Consequence = checkSplice(answer.CdnaDist)
	return err
}

// substitutionCoding determines the effects of a coding point mutation.
func substitutionCoding(g *Gene, cdnaIndexPos int, alt dna.Base, answer *EffectPrediction) error {
	answer.CdnaPos = cdnaIndexPos
	answer.AaPos = cdnaIndexPos / 3
	codon, err := CdnaPosToCodon(g, cdnaIndexPos)

	answer.AaRef = []dna.AminoAcid{dna.TranslateCodon(codon)}
	g.codingSeq.seq[cdnaIndexPos] = alt
	altCodon, altErr := CdnaPosToCodon(g, cdnaIndexPos)
	answer.AaAlt = []dna.AminoAcid{dna.TranslateCodon(altCodon)}

	if answer.AaRef[0] == answer.AaAlt[0] {
		answer.Consequence = Silent
	} else {
		if answer.AaAlt[0] == dna.Stop {
			answer.Consequence = Nonsense
			answer.StopDist = 0
		} else if answer.AaRef[0] == dna.Stop {
			answer.Consequence = DisruptStop
			//TODO calculate dist to new stop
		} else if answer.AaPos == 0 {
			answer.Consequence = DisruptStart
		} else {
			answer.Consequence = Missense
		}
	}

	if err != nil {
		return err
	} else {
		return altErr
	}
}

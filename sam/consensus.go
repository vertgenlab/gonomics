package sam

import (
	"fmt"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/numbers"
)

// consensusType represents the possible types of consensus variants in a Pile,
// which includes Bases, INDELs, or undefined (in the case where the Pile contains no data).
type consensusType byte

const (
	Base      consensusType = 0
	Insertion consensusType = 1
	Deletion  consensusType = 2
	Undefined consensusType = 3
)

// Consensus contains information on the most observed variant in a Pile struct.
type Consensus struct {
	Base      dna.Base
	Insertion []dna.Base
	Deletion  int
	Type      consensusType
}

// String for debug.
func (c Consensus) String() string {
	var typeString string
	switch c.Type {
	case 0:
		typeString = "Base"
	case 1:
		typeString = "Insertion"
	case 2:
		typeString = "Deletion"
	case 3:
		typeString = "Undefined"
	}
	return fmt.Sprintf("Base: %s. Insertion: %s. Deletion: %v. Type: %s.", dna.BaseToString(c.Base), dna.BasesToString(c.Insertion), c.Deletion, typeString)
}

// PileConsensus returns a Consensus struct from an input Pile struct
// representing information about the consensus variant observed in that Pile.
// An insertion is called as the consensus if the count number for the max insertion if
// maxCount > insertionThreshold * (total counts to bases).
func PileConsensus(p Pile, substitutionsOnly bool, insertionThreshold float64) Consensus {
	// first we check if consensus is a base
	var max int = p.CountF[dna.A] + p.CountR[dna.A]
	tiedConsensus := make([]Consensus, 1)
	tiedConsensus[0] = Consensus{Base: dna.A, Type: Base}

	max, tiedConsensus = getMaxBase(p, max, dna.C, tiedConsensus)
	max, tiedConsensus = getMaxBase(p, max, dna.G, tiedConsensus)
	max, tiedConsensus = getMaxBase(p, max, dna.T, tiedConsensus)
	if substitutionsOnly { //legacy feature
		//max, tiedConsensus = getMaxBase(p, max, dna.Gap, tiedConsensus) Dan I know that this is a legacy featuer but it's not technically substitutionsOnly, coordinates will change in the output
		if max < 1 {
			return Consensus{Type: Undefined}
		}
		return tiedConsensus[numbers.RandIntInRange(0, len(tiedConsensus))] //we don't need to check the length because this wil return tiedConsensus[0] even if length is 1.
	} else {
		max, tiedConsensus = getMaxDeletion(p, max, tiedConsensus)
		if max < 1 {
			return Consensus{Type: Undefined}
		}
		return getMaxInsertion(p, tiedConsensus, insertionThreshold)
	}
}

func getDeletionCounts(p Pile) int {
	var answer = 0
	var val int

	for _, val = range p.DelCountF {
		answer += val
	}
	for _, val = range p.DelCountR {
		answer += val
	}

	return answer
}

// helper function of PileConsensus, finds the max insertion in a Pile struct.
// InsThreshold is the proportion of Base reads for which an insertion needs to be observed to be called.
func getMaxInsertion(p Pile, tiedConsensus []Consensus, InsThreshold float64) Consensus {
	var count int
	var inMap bool
	var seenOnPosStrand = make(map[string]int, 0)
	var deletionSum = getDeletionCounts(p)
	var totalBaseCounts = p.CountF[0] + p.CountF[1] + p.CountF[2] + p.CountF[3] + p.CountR[0] + p.CountR[1] + p.CountR[2] + p.CountR[3] + deletionSum
	var insertionThreshold = int(InsThreshold * float64(totalBaseCounts))
	var maxInsertionScore, deletionScore = 0, 0

	for i := range p.InsCountF {
		seenOnPosStrand[i] = 1
		count = p.InsCountF[i]
		if _, inMap = p.InsCountR[i]; inMap {
			count += p.InsCountR[i]
		}
		switch tiedConsensus[0].Type {
		case Base:
			if count > insertionThreshold {
				tiedConsensus = tiedConsensus[:1]
				tiedConsensus[0].Type = Insertion
				tiedConsensus[0].Insertion = dna.StringToBases(i)
				maxInsertionScore = count
			} //note no ties for insertions. If the insertion is equal to the insertionThreshold, it is not called
		case Deletion:
			deletionScore = p.DelCountF[tiedConsensus[0].Deletion] + p.DelCountR[tiedConsensus[0].Deletion]
			if count > deletionScore {
				tiedConsensus = tiedConsensus[:1]
				tiedConsensus[0].Type = Insertion
				tiedConsensus[0].Insertion = dna.StringToBases(i)
				maxInsertionScore = count
			} //note no ties for insertions. If the insertion is equal to the insertionThreshold, it is not called
		case Insertion:
			if count > maxInsertionScore {
				tiedConsensus = tiedConsensus[:1]
				tiedConsensus[0].Insertion = dna.StringToBases(i)
				maxInsertionScore = count
			} else if count == maxInsertionScore {
				tiedConsensus = append(tiedConsensus, Consensus{Base: tiedConsensus[0].Base, Type: Insertion, Insertion: dna.StringToBases(i)})
			}
		default:
			return Consensus{Type: Undefined}
		}
	}

	//we have to check the p.InsCountR map for any insertions observed only on the minus strand.
	for i := range p.InsCountR {
		if _, inMap = seenOnPosStrand[i]; !inMap {
			count = p.InsCountR[i] //we don't have to sum since we are guaranteed not to have seen this insertion on the positive strand.
			switch tiedConsensus[0].Type {
			case Base:
				if count > insertionThreshold {
					tiedConsensus = tiedConsensus[:1]
					tiedConsensus[0].Type = Insertion
					tiedConsensus[0].Insertion = dna.StringToBases(i)
					maxInsertionScore = count
				}
			case Deletion:
				deletionScore = p.DelCountF[tiedConsensus[0].Deletion] + p.DelCountR[tiedConsensus[0].Deletion]
				if count > deletionScore {
					tiedConsensus = tiedConsensus[:1]
					tiedConsensus[0].Type = Insertion
					tiedConsensus[0].Insertion = dna.StringToBases(i)
					maxInsertionScore = count
				}
			case Insertion:
				if count > maxInsertionScore {
					tiedConsensus = tiedConsensus[:1]
					tiedConsensus[0].Insertion = dna.StringToBases(i)
					maxInsertionScore = count
				} else if count == maxInsertionScore {
					tiedConsensus = append(tiedConsensus, Consensus{Base: tiedConsensus[0].Base, Type: Insertion, Insertion: dna.StringToBases(i)})
				}
			default:
				return Consensus{Type: Undefined}
			}
		}
	}
	return tiedConsensus[numbers.RandIntInRange(0, len(tiedConsensus))]
}

// helper function of PileConsensus, finds the max deletion in a Pile struct.
func getMaxDeletion(p Pile, currMax int, tiedConsensus []Consensus) (int, []Consensus) {
	var count, i int
	var inMap bool
	var seenOnPosStrand = make(map[int]int, 0)
	for i = range p.DelCountF {
		seenOnPosStrand[i] = 1
		count = p.DelCountF[i]
		if _, inMap = p.DelCountR[i]; inMap {
			count += p.DelCountR[i]
		}
		if count > currMax {
			tiedConsensus = tiedConsensus[:1]
			tiedConsensus[0] = Consensus{Deletion: i, Type: Deletion}
			currMax = count
		}
		if count == currMax {
			tiedConsensus = append(tiedConsensus, Consensus{Deletion: i, Type: Deletion})
		}
	}

	for i = range p.DelCountR {
		if _, inMap = seenOnPosStrand[i]; !inMap {
			count = p.DelCountR[i]
			if count > currMax {
				tiedConsensus = tiedConsensus[:1]
				tiedConsensus[0] = Consensus{Deletion: i, Type: Deletion}
				currMax = count
			}
			if count == currMax {
				tiedConsensus = append(tiedConsensus, Consensus{Deletion: i, Type: Deletion})
			}
		}
	}
	return currMax, tiedConsensus
}

// getMaxBase is a helper function of PileConsensus, finds the consensus base in a Pile struct.
func getMaxBase(p Pile, currMax int, testBase dna.Base, tiedConsensus []Consensus) (int, []Consensus) {
	var count int = p.CountF[testBase] + p.CountR[testBase]
	if count > currMax {
		tiedConsensus = tiedConsensus[:1] // reset tied bases
		tiedConsensus[0] = Consensus{Base: testBase, Type: Base}
		return count, tiedConsensus
	}

	if count == currMax {
		tiedConsensus = append(tiedConsensus, Consensus{Base: testBase, Type: Base})
	}

	return currMax, tiedConsensus
}

package main

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/sam"
	"log"
	"math/rand"
)

func diploidInsertion(ans AnswerStruct, mlt MultiFaStruct, cacheStruct CacheStruct, p sam.Pile, refPos int, s Settings) (AnswerStruct, MultiFaStruct, CacheStruct, int) {
	var currRand = rand.Float64()
	var i int
	var currInsertion sam.DiploidInsertion
	var currInsertionSeqs [][]dna.Base

	// now we call diploid insertions and add the insertions to the output sequences
	currInsertion = sam.DiploidInsertionCallFromPile(p, cacheStruct.DiploidIndelPriorCache, cacheStruct.HomozygousIndelCache, cacheStruct.HeterozygousIndelCache, s.Epsilon)
	currInsertionSeqs = sam.DiploidInsertionToSeqs(currInsertion)
	refPos++

	switch currInsertion.Type {
	case sam.BBnoIns:
		// nothing to do, no insertion
	case sam.IaIa:
		for i = range currInsertionSeqs[0] {
			ans.AnswerA[ans.CurrFaIndex].Seq[ans.AnswerAPos] = currInsertionSeqs[0][i]
			ans.AnswerB[ans.CurrFaIndex].Seq[ans.AnswerBPos] = currInsertionSeqs[0][i]
			ans = advanceAnswerPos(ans)
		}
		for i = range currInsertionSeqs[0] {
			mlt = updateMultiFa(dna.Gap, currInsertionSeqs[0][i], currInsertionSeqs[0][i], mlt)
		}
	case sam.IaB:
		currRand = rand.Float64()
		for i = range currInsertionSeqs[0] {
			if currRand < 0.5 {
				ans = advanceAPos(ans)
				mlt = updateMultiFa(dna.Gap, currInsertionSeqs[0][i], dna.Gap, mlt)
			} else {
				ans = advanceBPos(ans)
				mlt = updateMultiFa(dna.Gap, dna.Gap, currInsertionSeqs[0][i], mlt)
			}
		}
	case sam.IaIb:
		currRand = rand.Float64()
		if currRand < 0.5 {
			for i = range currInsertionSeqs[0] {
				ans = advanceAPos(ans)
			}
			for i = range currInsertionSeqs[1] {
				ans = advanceBPos(ans)
			}
			for i = 0; i < numbers.Max(len(currInsertionSeqs[0]), len(currInsertionSeqs[1])); i++ { // for the length of the longer insertion
				mlt.CurrMultiFa[0].Seq[mlt.MultiFaPos] = dna.Gap
				if i < len(currInsertionSeqs[0]) {
					mlt.CurrMultiFa[1].Seq[mlt.MultiFaPos] = currInsertionSeqs[0][i]
				} else {
					mlt.CurrMultiFa[1].Seq[mlt.MultiFaPos] = dna.Gap
				}
				if i < len(currInsertionSeqs[1]) {
					mlt.CurrMultiFa[2].Seq[mlt.MultiFaPos] = currInsertionSeqs[1][i]
				} else {
					mlt.CurrMultiFa[2].Seq[mlt.MultiFaPos] = dna.Gap
				}
				mlt = advanceMultiFaPos(mlt)
			}
		} else {
			for i = range currInsertionSeqs[0] {
				ans = advanceBPos(ans)
			}
			for i = range currInsertionSeqs[1] {
				ans = advanceAPos(ans)
			}
			for i = 0; i < numbers.Max(len(currInsertionSeqs[0]), len(currInsertionSeqs[1])); i++ { // for the length of the longer insertion
				mlt.CurrMultiFa[0].Seq[mlt.MultiFaPos] = dna.Gap
				if i < len(currInsertionSeqs[0]) {
					mlt.CurrMultiFa[2].Seq[mlt.MultiFaPos] = currInsertionSeqs[0][i]
				} else {
					mlt.CurrMultiFa[2].Seq[mlt.MultiFaPos] = dna.Gap
				}
				if i < len(currInsertionSeqs[1]) {
					mlt.CurrMultiFa[1].Seq[mlt.MultiFaPos] = currInsertionSeqs[1][i]
				} else {
					mlt.CurrMultiFa[1].Seq[mlt.MultiFaPos] = dna.Gap
				}
				mlt = advanceMultiFaPos(mlt)
			}
		}
	}
	return ans, mlt, cacheStruct, refPos
}

func diploidDeletion(mlt MultiFaStruct, cacheStruct CacheStruct, p sam.Pile, refMap map[string][]dna.Base, refPos int, currChrom string, s Settings) (MultiFaStruct, CacheStruct, int, bool, int, int, int) {
	var currDeletion = sam.DiploidDeletionCallFromPile(p, cacheStruct.DiploidIndelPriorCache, cacheStruct.HomozygousIndelCache, cacheStruct.HeterozygousIndelCache, s.Epsilon)
	var i, currPloidy, haploidBases, positionsToSkip int
	var currRand float64
	var haploidStrand bool
	currPloidy = 2 //default diploid, unless changed to 1 later on.

	switch currDeletion.Type {
	case sam.BBNoDel:
		// no deletion, nothing to do
	case sam.DaDa:
		// we'll skip the next piles corresponding to the number of homozygous deleted bases
		positionsToSkip = currDeletion.Da
	case sam.DaB:
		currPloidy = 1
		haploidBases = currDeletion.Da
		currRand = rand.Float64()
		if currRand < 0.5 { //we randomly assign the haploid bases to a strand
			haploidStrand = true
		} else {
			haploidStrand = false
		}
	case sam.DaDb:
		currPloidy = 1
		// we first advance through the shared part of the deletion (the length of the shorter deletion)
		for i = 0; i < numbers.Min(currDeletion.Da, currDeletion.Db); i++ {
			mlt = updateMultiFa(refMap[currChrom][refPos], dna.Gap, dna.Gap, mlt)
			refPos++
		}
		// determine the haploid bases to be the difference in length between the two deletions
		haploidBases = numbers.Max(currDeletion.Da-currDeletion.Db, currDeletion.Db-currDeletion.Da)
		currRand = rand.Float64()
		if currRand < 0.5 { // we randomly assign the haploid bases to a strand
			haploidStrand = true
		} else {
			haploidStrand = false
		}
	default:
		log.Fatalf("Unrecognized deletion type: %v.\n", currDeletion.Type)
	}
	return mlt, cacheStruct, refPos, haploidStrand, currPloidy, haploidBases, positionsToSkip
}

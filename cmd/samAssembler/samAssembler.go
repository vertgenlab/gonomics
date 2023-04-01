// Command Group: "Data Conversion"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/sam"
	"log"
	"math/rand"
)

type Settings struct {
	SamFileName         string
	RefFile             string
	OutFileA            string
	OutFileB            string
	MultiFaDir          string
	tName               string
	qNameA              string
	qNameB              string
	Delta               float64
	Gamma               float64
	Epsilon             float64
	Kappa               float64
	LikelihoodCacheSize int
	SetSeed             int64
}

const bufferSize = 10_000_000

func samAssembler(s Settings) {
	rand.Seed(s.SetSeed)
	var i, currFaIndex, emptyRoomInBufferA, emptyRoomInBufferB, multiFaPos int
	var emptyRoomInMultiFaBuffer, answerAPos, answerBPos, refPos, positionsToSkip int
	var haploidBases int
	var haploidStrand bool //haploidStrand is true when the haploid bases are on the first strand (deletion on second strand)
	var currChrom string
	var currMultiFa []fasta.Fasta
	var newBufferRoom = make([]dna.Base, bufferSize)
	var firstTime = true
	var currPloidy = 2
	var currDiploidBases, currHaploidBases []dna.Base
	var currInsertion sam.DiploidInsertion
	var currInsertionSeqs [][]dna.Base
	var currDeletion sam.DiploidDeletion
	var currHaploidCall sam.HaploidCall
	var currRand float64

	// initialize caches for likelihood and priors
	var diploidBasePriorCache [][]float64 = sam.MakeDiploidBasePriorCache(s.Delta, s.Gamma)
	var diploidIndelPriorCache []float64 = sam.MakeDiploidIndelPriorCache(s.Kappa, s.Delta)
	var haploidBasePriorCache [][]float64 = sam.MakeHaploidBasePriorCache(s.Delta, s.Gamma)
	var haploidIndelPriorCache []float64 = sam.MakeHaploidIndelPriorCache(s.Delta, s.Kappa)
	var homozygousBaseCache [][]float64 = make([][]float64, s.LikelihoodCacheSize)
	for i = range homozygousBaseCache {
		homozygousBaseCache[i] = make([]float64, s.LikelihoodCacheSize)
	}
	var heterozygousBaseCache = make([][]float64, s.LikelihoodCacheSize)
	for i = range heterozygousBaseCache {
		heterozygousBaseCache[i] = make([]float64, s.LikelihoodCacheSize)
	}
	var homozygousIndelCache = make([][]float64, s.LikelihoodCacheSize)
	for i = range homozygousIndelCache {
		homozygousIndelCache[i] = make([]float64, s.LikelihoodCacheSize)
	}
	var heterozygousIndelCache = make([][]float64, s.LikelihoodCacheSize)
	for i = range heterozygousIndelCache {
		heterozygousIndelCache[i] = make([]float64, s.LikelihoodCacheSize)
	}

	// initialize reference genome
	ref := fasta.Read(s.RefFile)
	for i = range ref {
		dna.AllToUpper(ref[i].Seq)
	}
	refMap := fasta.ToMap(ref)

	// read pileups from sam/bam
	reads, header := sam.GoReadToChan(s.SamFileName)
	piles := sam.GoPileup(reads, header, false, nil, nil)

	//initialize output
	var answerA []fasta.Fasta = make([]fasta.Fasta, len(ref))
	var answerB []fasta.Fasta = make([]fasta.Fasta, len(ref))
	for i = range ref {
		answerA[i] = fasta.Fasta{Name: ref[i].Name, Seq: make([]dna.Base, bufferSize)}
		answerB[i] = fasta.Fasta{Name: ref[i].Name, Seq: make([]dna.Base, bufferSize)}
	}

	//now time for the main loop, we look through each pile
	for p := range piles {
		if positionsToSkip > 0 {
			if s.MultiFaDir != "" {
				currMultiFa[0].Seq[multiFaPos] = refMap[currChrom][refPos]
				currMultiFa[1].Seq[multiFaPos] = dna.Gap
				currMultiFa[2].Seq[multiFaPos] = dna.Gap
				multiFaPos, emptyRoomInMultiFaBuffer, currMultiFa = advanceMultiFaPos(multiFaPos, emptyRoomInMultiFaBuffer, currMultiFa, newBufferRoom)
			}
			positionsToSkip--
		}
		if firstTime {
			firstTime = false
			currChrom = header.Chroms[p.RefIdx].Name
			currFaIndex = getIndexForName(answerA, currChrom)
			answerAPos, answerBPos, refPos = 0, 0, 0
			emptyRoomInBufferA, emptyRoomInBufferB = bufferSize, bufferSize
			if s.MultiFaDir != "" {
				currMultiFa = []fasta.Fasta{{Name: s.tName, Seq: make([]dna.Base, bufferSize)},
					{Name: s.qNameA, Seq: make([]dna.Base, bufferSize)},
					{Name: s.qNameB, Seq: make([]dna.Base, bufferSize)}}
				emptyRoomInMultiFaBuffer = bufferSize
				multiFaPos = 0
			}
		}
		if currChrom != header.Chroms[p.RefIdx].Name { //if we've moved onto a new chromosome.
			for refPos < len(refMap[currChrom]) { //write out the rest of the current reference
				answerA[currFaIndex].Seq[answerAPos] = refMap[currChrom][refPos]
				answerB[currFaIndex].Seq[answerBPos] = refMap[currChrom][refPos]
				if s.MultiFaDir != "" {
					currMultiFa[0].Seq[multiFaPos] = refMap[currChrom][refPos]
					currMultiFa[1].Seq[multiFaPos] = refMap[currChrom][refPos]
					currMultiFa[2].Seq[multiFaPos] = refMap[currChrom][refPos]
					multiFaPos, emptyRoomInMultiFaBuffer, currMultiFa = advanceMultiFaPos(multiFaPos, emptyRoomInMultiFaBuffer, currMultiFa, newBufferRoom)
				}
				refPos++
				answerAPos, answerBPos, emptyRoomInBufferA, emptyRoomInBufferB, currFaIndex, answerA, answerB = advanceAnswerPos(answerAPos, answerBPos, emptyRoomInBufferA, emptyRoomInBufferB, currFaIndex, answerA, answerB, newBufferRoom)
			}
			//write out current chromosome
			answerA[currFaIndex].Seq = answerA[currFaIndex].Seq[:len(answerA[currFaIndex].Seq)-emptyRoomInBufferA] //clear out empty buffer positions
			answerB[currFaIndex].Seq = answerB[currFaIndex].Seq[:len(answerB[currFaIndex].Seq)-emptyRoomInBufferB]
			if s.MultiFaDir != "" {
				currMultiFa[0].Seq = currMultiFa[0].Seq[:len(currMultiFa[0].Seq)-emptyRoomInMultiFaBuffer]
				currMultiFa[1].Seq = currMultiFa[1].Seq[:len(currMultiFa[1].Seq)-emptyRoomInMultiFaBuffer]
				currMultiFa[2].Seq = currMultiFa[2].Seq[:len(currMultiFa[2].Seq)-emptyRoomInMultiFaBuffer]
				fasta.Write(fmt.Sprintf("%s/%s.fa", s.MultiFaDir, currChrom), currMultiFa)
				currMultiFa = []fasta.Fasta{{Name: s.tName, Seq: make([]dna.Base, bufferSize)},
					{Name: s.qNameA, Seq: make([]dna.Base, bufferSize)},
					{Name: s.qNameB, Seq: make([]dna.Base, bufferSize)}}
				emptyRoomInMultiFaBuffer = bufferSize
				multiFaPos = 0
			}
			//now we set up the new chromosome
			currChrom = header.Chroms[p.RefIdx].Name
			currFaIndex = getIndexForName(answerA, currChrom)
			emptyRoomInBufferA, emptyRoomInBufferB = bufferSize, bufferSize
			answerAPos, answerBPos, refPos = 0, 0, 0
		}

		//catch up to the current pile position, handles reference positions with no Pile coverage.
		for refPos < int(p.Pos-1) {
			answerA[currFaIndex].Seq[answerAPos] = refMap[currChrom][refPos]
			answerB[currFaIndex].Seq[answerBPos] = refMap[currChrom][refPos]
			if s.MultiFaDir != "" {
				currMultiFa[0].Seq[multiFaPos] = refMap[currChrom][refPos]
				currMultiFa[1].Seq[multiFaPos] = refMap[currChrom][refPos]
				currMultiFa[2].Seq[multiFaPos] = refMap[currChrom][refPos]
				multiFaPos, emptyRoomInMultiFaBuffer, currMultiFa = advanceMultiFaPos(multiFaPos, emptyRoomInMultiFaBuffer, currMultiFa, newBufferRoom)
			}
			answerAPos, answerBPos, emptyRoomInBufferA, emptyRoomInBufferB, currFaIndex, answerA, answerB = advanceAnswerPos(answerAPos, answerBPos, emptyRoomInBufferA, emptyRoomInBufferB, currFaIndex, answerA, answerB, newBufferRoom)
			refPos++
		}

		//now refPos should equal p.Pos - 1, because of our for loop before
		if refPos != int(p.Pos-1) {
			log.Fatalf("Something went wrong. RefPos is not equal to p.Pos -1.")
		}

		if currPloidy == 2 {
			//First we handle the base call for the current pile
			currDiploidBases = sam.DiploidBaseToBases(sam.DiploidBaseCallFromPile(p, refMap[currChrom][refPos], diploidBasePriorCache, homozygousBaseCache, heterozygousBaseCache, s.Epsilon))
			currRand = rand.Float64()
			if currRand < 0.5 {
				answerA[currFaIndex].Seq[answerAPos] = currDiploidBases[0]
				answerB[currFaIndex].Seq[answerBPos] = currDiploidBases[1]
				if s.MultiFaDir != "" {
					currMultiFa[0].Seq[multiFaPos] = refMap[currChrom][refPos]
					currMultiFa[1].Seq[multiFaPos] = currDiploidBases[0]
					currMultiFa[2].Seq[multiFaPos] = currDiploidBases[1]
					multiFaPos, emptyRoomInMultiFaBuffer, currMultiFa = advanceMultiFaPos(multiFaPos, emptyRoomInMultiFaBuffer, currMultiFa, newBufferRoom)
				}
			} else {
				answerA[currFaIndex].Seq[answerAPos] = currDiploidBases[1]
				answerB[currFaIndex].Seq[answerBPos] = currDiploidBases[0]
				if s.MultiFaDir != "" {
					currMultiFa[0].Seq[multiFaPos] = refMap[currChrom][refPos]
					currMultiFa[1].Seq[multiFaPos] = currDiploidBases[1]
					currMultiFa[2].Seq[multiFaPos] = currDiploidBases[0]
					multiFaPos, emptyRoomInMultiFaBuffer, currMultiFa = advanceMultiFaPos(multiFaPos, emptyRoomInMultiFaBuffer, currMultiFa, newBufferRoom)
				}
			}
			answerAPos, answerBPos, emptyRoomInBufferA, emptyRoomInBufferB, currFaIndex, answerA, answerB = advanceAnswerPos(answerAPos, answerBPos, emptyRoomInBufferA, emptyRoomInBufferB, currFaIndex, answerA, answerB, newBufferRoom)

			// now we call diploid insertions and add the insertions to the output sequences
			currInsertion = sam.DiploidInsertionCallFromPile(p, diploidIndelPriorCache, homozygousIndelCache, heterozygousIndelCache, s.Epsilon)
			currInsertionSeqs = sam.DiploidInsertionToSeqs(currInsertion)
			refPos++

			switch currInsertion.Type {
			case sam.BBnoIns:
				//nothing to do, no insertion
			case sam.IaIa:
				for i = range currInsertionSeqs[0] {
					answerA[currFaIndex].Seq[answerAPos] = currInsertionSeqs[0][i]
					answerB[currFaIndex].Seq[answerBPos] = currInsertionSeqs[0][i]
					answerAPos, answerBPos, emptyRoomInBufferA, emptyRoomInBufferB, currFaIndex, answerA, answerB = advanceAnswerPos(answerAPos, answerBPos, emptyRoomInBufferA, emptyRoomInBufferB, currFaIndex, answerA, answerB, newBufferRoom)
				}
				if s.MultiFaDir != "" {
					for i = range currInsertionSeqs[0] {
						currMultiFa[0].Seq[multiFaPos] = dna.Gap
						currMultiFa[1].Seq[multiFaPos] = currInsertionSeqs[0][i]
						currMultiFa[2].Seq[multiFaPos] = currInsertionSeqs[0][i]
						multiFaPos, emptyRoomInMultiFaBuffer, currMultiFa = advanceMultiFaPos(multiFaPos, emptyRoomInMultiFaBuffer, currMultiFa, newBufferRoom)
					}
				}
			case sam.IaB:
				currRand = rand.Float64()
				if currRand < 0.5 {
					for i = range currInsertionSeqs[0] {
						answerA[currFaIndex].Seq[answerAPos] = currInsertionSeqs[0][i]
						answerAPos++
						emptyRoomInBufferA--
						if emptyRoomInBufferA < 1 {
							answerA[currFaIndex].Seq = append(answerA[currFaIndex].Seq, newBufferRoom...)
							emptyRoomInBufferA += bufferSize
						}
					}
					if s.MultiFaDir != "" {
						for i = range currInsertionSeqs[0] {
							currMultiFa[0].Seq[multiFaPos] = dna.Gap
							currMultiFa[1].Seq[multiFaPos] = currInsertionSeqs[0][i]
							currMultiFa[2].Seq[multiFaPos] = dna.Gap
							multiFaPos, emptyRoomInMultiFaBuffer, currMultiFa = advanceMultiFaPos(multiFaPos, emptyRoomInMultiFaBuffer, currMultiFa, newBufferRoom)
						}
					}
				} else {
					for i = range currInsertionSeqs[0] {
						answerB[currFaIndex].Seq[answerBPos] = currInsertionSeqs[0][i]
						answerBPos++
						emptyRoomInBufferB--
						if emptyRoomInBufferB < 1 {
							answerB[currFaIndex].Seq = append(answerB[currFaIndex].Seq, newBufferRoom...)
							emptyRoomInBufferB += bufferSize
						}
						if s.MultiFaDir != "" {
							currMultiFa[0].Seq[multiFaPos] = dna.Gap
							currMultiFa[1].Seq[multiFaPos] = dna.Gap
							currMultiFa[2].Seq[multiFaPos] = currInsertionSeqs[0][i]
							multiFaPos, emptyRoomInMultiFaBuffer, currMultiFa = advanceMultiFaPos(multiFaPos, emptyRoomInMultiFaBuffer, currMultiFa, newBufferRoom)
						}
					}
				}
			case sam.IaIb:
				currRand = rand.Float64()
				if currRand < 0.5 {
					for i = range currInsertionSeqs[0] {
						answerA[currFaIndex].Seq[answerAPos] = currInsertionSeqs[0][i]
						answerAPos++
						emptyRoomInBufferA--
						if emptyRoomInBufferA < 1 {
							answerA[currFaIndex].Seq = append(answerA[currFaIndex].Seq, newBufferRoom...)
							emptyRoomInBufferA += bufferSize
						}
					}
					for i = range currInsertionSeqs[1] {
						answerB[currFaIndex].Seq[answerBPos] = currInsertionSeqs[1][i]
						answerBPos++
						emptyRoomInBufferB--
						if emptyRoomInBufferB < 1 {
							answerB[currFaIndex].Seq = append(answerB[currFaIndex].Seq, newBufferRoom...)
							emptyRoomInBufferB += bufferSize
						}
					}
					if s.MultiFaDir != "" {
						for i = 0; i < numbers.Max(len(currInsertionSeqs[0]), len(currInsertionSeqs[1])); i++ { //for the length of the longer insertion
							currMultiFa[0].Seq[multiFaPos] = dna.Gap
							if i < len(currInsertionSeqs[0]) {
								currMultiFa[1].Seq[multiFaPos] = currInsertionSeqs[0][i]
							} else {
								currMultiFa[1].Seq[multiFaPos] = dna.Gap
							}
							if i < len(currInsertionSeqs[1]) {
								currMultiFa[2].Seq[multiFaPos] = currInsertionSeqs[1][i]
							} else {
								currMultiFa[2].Seq[multiFaPos] = dna.Gap
							}
							multiFaPos, emptyRoomInMultiFaBuffer, currMultiFa = advanceMultiFaPos(multiFaPos, emptyRoomInMultiFaBuffer, currMultiFa, newBufferRoom)
						}
					}
				} else {
					for i = range currInsertionSeqs[0] {
						answerB[currFaIndex].Seq[answerBPos] = currInsertionSeqs[0][i]
						answerBPos++
						emptyRoomInBufferB--
						if emptyRoomInBufferB < 1 {
							answerB[currFaIndex].Seq = append(answerB[currFaIndex].Seq, newBufferRoom...)
							emptyRoomInBufferB += bufferSize
						}
					}
					for i = range currInsertionSeqs[1] {
						answerA[currFaIndex].Seq[answerAPos] = currInsertionSeqs[1][i]
						answerAPos++
						emptyRoomInBufferA--
						if emptyRoomInBufferA < 1 {
							answerA[currFaIndex].Seq = append(answerA[currFaIndex].Seq, newBufferRoom...)
							emptyRoomInBufferA += bufferSize
						}
					}
					if s.MultiFaDir != "" {
						for i = 0; i < numbers.Max(len(currInsertionSeqs[0]), len(currInsertionSeqs[1])); i++ { //for the length of the longer insertion
							currMultiFa[0].Seq[multiFaPos] = dna.Gap
							if i < len(currInsertionSeqs[0]) {
								currMultiFa[2].Seq[multiFaPos] = currInsertionSeqs[0][i]
							} else {
								currMultiFa[2].Seq[multiFaPos] = dna.Gap
							}
							if i < len(currInsertionSeqs[1]) {
								currMultiFa[1].Seq[multiFaPos] = currInsertionSeqs[1][i]
							} else {
								currMultiFa[1].Seq[multiFaPos] = dna.Gap
							}
							multiFaPos, emptyRoomInMultiFaBuffer, currMultiFa = advanceMultiFaPos(multiFaPos, emptyRoomInMultiFaBuffer, currMultiFa, newBufferRoom)
						}
					}
				}
			}

			//Now we handle diploid deletion calls
			currDeletion = sam.DiploidDeletionCallFromPile(p, diploidIndelPriorCache, homozygousIndelCache, heterozygousIndelCache, s.Epsilon)

			switch currDeletion.Type {
			case sam.BBNoDel:
				//no deletion, nothing to do
			case sam.DaDa:
				//we'll skip the next piles corresponding to the number of homozygous deleted bases
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
				//we first advance through the shared part of the deletion (the length of the shorter deletion)
				for i = 0; i < numbers.Min(currDeletion.Da, currDeletion.Db); i++ {
					if s.MultiFaDir != "" {
						currMultiFa[0].Seq[multiFaPos] = refMap[currChrom][refPos]
						currMultiFa[1].Seq[multiFaPos] = dna.Gap
						currMultiFa[2].Seq[multiFaPos] = dna.Gap
					}
					refPos++
				}
				// determine the haploid bases to be the difference in length between the two deletions
				haploidBases = numbers.Max(currDeletion.Da-currDeletion.Db, currDeletion.Db-currDeletion.Da)
				currRand = rand.Float64()
				if currRand < 0.5 { //we randomly assign the haploid bases to a strand
					haploidStrand = true
				} else {
					haploidStrand = false
				}
			default:
				log.Fatalf("Unrecognized deletion type: %v.\n", currDeletion.Type)
			}

		} else if currPloidy == 1 {
			currHaploidCall = sam.HaploidCallFromPile(p, refMap[currChrom][refPos], s.Epsilon, haploidBasePriorCache, haploidIndelPriorCache, homozygousBaseCache, heterozygousBaseCache, homozygousIndelCache)

			if haploidStrand {
				answerA[currFaIndex].Seq[answerAPos] = currHaploidCall.Base
				answerAPos++
				emptyRoomInBufferA--
				if emptyRoomInBufferA < 1 {
					answerA[currFaIndex].Seq = append(answerA[currFaIndex].Seq, newBufferRoom...)
					emptyRoomInBufferA += bufferSize
				}
				if s.MultiFaDir != "" {
					currMultiFa[0].Seq[multiFaPos] = refMap[currChrom][refPos]
					currMultiFa[1].Seq[multiFaPos] = currHaploidCall.Base
					currMultiFa[2].Seq[multiFaPos] = dna.Gap
					multiFaPos, emptyRoomInMultiFaBuffer, currMultiFa = advanceMultiFaPos(multiFaPos, emptyRoomInMultiFaBuffer, currMultiFa, newBufferRoom)
				}
				if currHaploidCall.Insertion != "" {
					currHaploidBases = dna.StringToBases(currHaploidCall.Insertion)
					for i = 0; i < len(currHaploidBases); i++ {
						answerA[currFaIndex].Seq[answerAPos] = currHaploidBases[i]
						answerAPos++
						emptyRoomInBufferA--
						if emptyRoomInBufferA < 1 {
							answerA[currFaIndex].Seq = append(answerA[currFaIndex].Seq, newBufferRoom...)
							emptyRoomInBufferA += bufferSize
						}
						if s.MultiFaDir != "" {
							currMultiFa[0].Seq[multiFaPos] = dna.Gap
							currMultiFa[1].Seq[multiFaPos] = currHaploidBases[i]
							currMultiFa[2].Seq[multiFaPos] = dna.Gap
							multiFaPos, emptyRoomInMultiFaBuffer, currMultiFa = advanceMultiFaPos(multiFaPos, emptyRoomInMultiFaBuffer, currMultiFa, newBufferRoom)
						}
					}
				}
				if currHaploidCall.Deletion != 0 {
					for i = 0; i < currHaploidCall.Deletion; i++ {
						if s.MultiFaDir != "" {
							currMultiFa[0].Seq[multiFaPos] = refMap[currChrom][refPos]
							currMultiFa[1].Seq[multiFaPos] = dna.Gap
							currMultiFa[2].Seq[multiFaPos] = dna.Gap
							multiFaPos, emptyRoomInMultiFaBuffer, currMultiFa = advanceMultiFaPos(multiFaPos, emptyRoomInMultiFaBuffer, currMultiFa, newBufferRoom)
						}
						refPos++
						if refPos >= len(refMap[currChrom]) {
							currPloidy = 2
							break
						}
						haploidBases--
						if haploidBases < 1 {
							currPloidy = 2
							break
						}
					}
				}
			} else {
				answerB[currFaIndex].Seq[answerBPos] = currHaploidCall.Base
				answerBPos++
				emptyRoomInBufferB--
				if emptyRoomInBufferB < 1 {
					answerB[currFaIndex].Seq = append(answerB[currFaIndex].Seq, newBufferRoom...)
					emptyRoomInBufferB += bufferSize
				}
				if s.MultiFaDir != "" {
					currMultiFa[0].Seq[multiFaPos] = refMap[currChrom][refPos]
					currMultiFa[1].Seq[multiFaPos] = dna.Gap
					currMultiFa[2].Seq[multiFaPos] = currHaploidCall.Base
					multiFaPos, emptyRoomInMultiFaBuffer, currMultiFa = advanceMultiFaPos(multiFaPos, emptyRoomInMultiFaBuffer, currMultiFa, newBufferRoom)
				}

				if currHaploidCall.Insertion != "" {
					currHaploidBases = dna.StringToBases(currHaploidCall.Insertion)
					for i = 0; i < len(currHaploidBases); i++ {
						answerB[currFaIndex].Seq[answerBPos] = currHaploidBases[i]
						answerBPos++
						emptyRoomInBufferB--
						if emptyRoomInBufferB < 1 {
							answerB[currFaIndex].Seq = append(answerB[currFaIndex].Seq, newBufferRoom...)
							emptyRoomInBufferB += bufferSize
						}
						if s.MultiFaDir != "" {
							currMultiFa[0].Seq[multiFaPos] = dna.Gap
							currMultiFa[1].Seq[multiFaPos] = dna.Gap
							currMultiFa[2].Seq[multiFaPos] = currHaploidBases[i]
							multiFaPos, emptyRoomInMultiFaBuffer, currMultiFa = advanceMultiFaPos(multiFaPos, emptyRoomInMultiFaBuffer, currMultiFa, newBufferRoom)
						}
					}
				}
				if currHaploidCall.Deletion != 0 {
					for i = 0; i < currHaploidCall.Deletion; i++ {
						if s.MultiFaDir != "" {
							currMultiFa[0].Seq[multiFaPos] = refMap[currChrom][refPos]
							currMultiFa[1].Seq[multiFaPos] = dna.Gap
							currMultiFa[2].Seq[multiFaPos] = dna.Gap
							multiFaPos, emptyRoomInMultiFaBuffer, currMultiFa = advanceMultiFaPos(multiFaPos, emptyRoomInMultiFaBuffer, currMultiFa, newBufferRoom)
						}
						refPos++
						if refPos >= len(refMap[currChrom]) {
							currPloidy = 2
							break
						}
						haploidBases--
						if haploidBases < 1 {
							currPloidy = 2
							break
						}
					}
				}
			}

			if haploidBases < 2 { //if we are on the last haploidBase, we re-enter diploid mode
				currPloidy = 2
			}
			refPos++
			haploidBases--
		} else {
			log.Fatalf("Error in samAssembly. Unrecognized ploidy: %v.\n", currPloidy)
		}
	}

	//once we're done with the piles we have to add the trailing ref bases and clear the buffer for the last chrom
	for refPos < len(refMap[currChrom]) {
		answerA[currFaIndex].Seq[answerAPos] = refMap[currChrom][refPos]
		answerB[currFaIndex].Seq[answerBPos] = refMap[currChrom][refPos]
		if s.MultiFaDir != "" {
			currMultiFa[0].Seq[multiFaPos] = refMap[currChrom][refPos]
			currMultiFa[1].Seq[multiFaPos] = refMap[currChrom][refPos]
			currMultiFa[2].Seq[multiFaPos] = refMap[currChrom][refPos]
			multiFaPos, emptyRoomInMultiFaBuffer, currMultiFa = advanceMultiFaPos(multiFaPos, emptyRoomInMultiFaBuffer, currMultiFa, newBufferRoom)
		}
		answerAPos, answerBPos, emptyRoomInBufferA, emptyRoomInBufferB, currFaIndex, answerA, answerB = advanceAnswerPos(answerAPos, answerBPos, emptyRoomInBufferA, emptyRoomInBufferB, currFaIndex, answerA, answerB, newBufferRoom)
		refPos++
	}

	//write last multiFa entry
	if s.MultiFaDir != "" {
		currMultiFa[0].Seq = currMultiFa[0].Seq[:len(currMultiFa[0].Seq)-emptyRoomInMultiFaBuffer]
		currMultiFa[1].Seq = currMultiFa[1].Seq[:len(currMultiFa[1].Seq)-emptyRoomInMultiFaBuffer]
		currMultiFa[2].Seq = currMultiFa[2].Seq[:len(currMultiFa[2].Seq)-emptyRoomInMultiFaBuffer]
		fasta.Write(fmt.Sprintf("%s/%s.fa", s.MultiFaDir, currChrom), currMultiFa)
	}

	//write answer
	answerA[currFaIndex].Seq = answerA[currFaIndex].Seq[:len(answerA[currFaIndex].Seq)-emptyRoomInBufferA]
	answerB[currFaIndex].Seq = answerB[currFaIndex].Seq[:len(answerB[currFaIndex].Seq)-emptyRoomInBufferB]
	fasta.Write(s.OutFileA, answerA)
	fasta.Write(s.OutFileB, answerB)
}

func advanceAnswerPos(answerAPos int, answerBPos int, emptyRoomInBufferA int, emptyRoomInBufferB int, currFaIndex int, answerA []fasta.Fasta, answerB []fasta.Fasta, newBufferRoom []dna.Base) (int, int, int, int, int, []fasta.Fasta, []fasta.Fasta) {
	answerAPos++
	answerBPos++
	emptyRoomInBufferA--
	emptyRoomInBufferB--
	if emptyRoomInBufferA < 1 {
		answerA[currFaIndex].Seq = append(answerA[currFaIndex].Seq, newBufferRoom...)
		emptyRoomInBufferA += bufferSize
	}
	if emptyRoomInBufferB < 1 {
		answerB[currFaIndex].Seq = append(answerB[currFaIndex].Seq, newBufferRoom...)
		emptyRoomInBufferB += bufferSize
	}
	return answerAPos, answerBPos, emptyRoomInBufferA, emptyRoomInBufferB, currFaIndex, answerA, answerB
}

func advanceMultiFaPos(multiFaPos int, emptyRoomInMultiFaBuffer int, currMultiFa []fasta.Fasta, newBufferRoom []dna.Base) (int, int, []fasta.Fasta) {
	multiFaPos++
	emptyRoomInMultiFaBuffer--
	if emptyRoomInMultiFaBuffer < 1 {
		currMultiFa[0].Seq = append(currMultiFa[0].Seq, newBufferRoom...)
		currMultiFa[1].Seq = append(currMultiFa[1].Seq, newBufferRoom...)
		currMultiFa[2].Seq = append(currMultiFa[2].Seq, newBufferRoom...)
		emptyRoomInMultiFaBuffer += bufferSize
	}
	return multiFaPos, emptyRoomInMultiFaBuffer, currMultiFa
}

func getIndexForName(f []fasta.Fasta, name string) int {
	for i := range f {
		if f[i].Name == name {
			return i
		}
	}
	log.Fatalf("Name: %s not found in fasta.", name)
	return -1
}

func usage() {
	fmt.Print(
		"samAssembler - Reference-based diploid assembly of aligned short reads.\n" +
			"Usage:\n" +
			"samAssembler individual.sam ref.fa outputA.fa outputB.fa\n" +
			"options:\n")
}

func main() {
	var expectedNumArgs int = 4
	var delta *float64 = flag.Float64("delta", 0.01, "Set the expected divergence frequency.")
	var gamma *float64 = flag.Float64("gamma", 3, "Set the expected transition bias.")
	var epsilon *float64 = flag.Float64("epsilon", 0.01, "Set the expected misclassification error rate.")
	var kappa *float64 = flag.Float64("kappa", 0.1, "Set the expected proportion of divergent sites that are INDELs.")
	var multiFaDir *string = flag.String("multiFaDir", "", "Output the reference and generated sequences as an aligned multiFa, each file by chrom.")
	var tName *string = flag.String("tName", "Target", "Set the tName in the optional multiFa output.")
	var qNameA *string = flag.String("qNameA", "QueryA", "Set the qName for the first generated chromosome in the optional multiFa output.")
	var qNameB *string = flag.String("qNameB", "QueryB", "Set the qName for the second generated chromosome in the optional multiFa output.")
	var likelihoodCacheSize *int = flag.Int("likelihoodCacheSize", 100, "Set the maximum dimension of the likelihood caches. Should be slightly larger than highest expected pile depth.")
	var setSeed *int64 = flag.Int64("setSeed", -1, "Use a specific seed for the RNG.")
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)
	refFile := flag.Arg(1)
	outFileA := flag.Arg(2)
	outFileB := flag.Arg(3)

	s := Settings{
		SamFileName:         inFile,
		RefFile:             refFile,
		OutFileA:            outFileA,
		OutFileB:            outFileB,
		MultiFaDir:          *multiFaDir,
		tName:               *tName,
		qNameA:              *qNameA,
		qNameB:              *qNameB,
		Delta:               *delta,
		Gamma:               *gamma,
		Epsilon:             *epsilon,
		Kappa:               *kappa,
		LikelihoodCacheSize: *likelihoodCacheSize,
		SetSeed:             *setSeed,
	}

	samAssembler(s)
}

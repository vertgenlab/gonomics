package motif

import (
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"log"
	"math"
)

type MatchCompSettings struct {
	MotifFile          string
	MotifType          string
	Records            []fasta.Fasta
	PropMatch          float64
	ChromName          string
	OutFile            string
	Pseudocounts       float64
	ResidualWindowSize int
	RefStart           int
	OutputAsProportion bool
	EnforceStrandMatch bool
	ResidualFilter     float64
	GcContent          float64
	MatrixFilter       bool
}

func MatchComp(s MatchCompSettings) {
	var err error
	var motifLen int
	var kmerHash, revKmerHash map[uint64]float64
	var consensusScore float64
	var couldScoreConsensus bool
	var revCompMotif PositionMatrix

	var motifsUnfiltered, motifs []PositionMatrix
	switch s.MotifType {
	case "Frequency":
		motifsUnfiltered = ReadJaspar(s.MotifFile, "Frequency")
		motifsUnfiltered = PfmSliceToPpmSlice(motifsUnfiltered, s.Pseudocounts)
		motifsUnfiltered = PpmSliceToPwmSlice(motifsUnfiltered, s.GcContent)
	case "Probability":
		motifsUnfiltered = ReadJaspar(s.MotifFile, "Probability")
		motifsUnfiltered = PpmSliceToPwmSlice(motifsUnfiltered, s.GcContent)
	case "Weight":
		motifsUnfiltered = ReadJaspar(s.MotifFile, "Weight")
	default:
		log.Fatalf("Error. Unexpected motif file format. Options are 'Frequency', 'Probability', and 'Weight'.")
	}

	out := fileio.EasyCreate(s.OutFile)

	// when using the matrixFilter option, filter motifs to only retain motifs with length <= 32
	if s.MatrixFilter {
		for i := range motifsUnfiltered {
			motifLen = len(motifsUnfiltered[i].Mat[0])
			if motifLen <= 32 {
				motifs = append(motifs, motifsUnfiltered[i])
			} else {
				fmt.Printf("Filtered out matrix with motif length greater than 32. Matrix ID: %v. Motif length: %v.\n", motifsUnfiltered[i].Id, motifLen)
			}
		}
	} else {
		motifs = make([]PositionMatrix, len(motifsUnfiltered))
		copy(motifs, motifsUnfiltered)
	}

	for i := range motifs {
		motifLen = len(motifs[i].Mat[0])
		if motifLen > 32 {
			log.Fatalf("Error: MatchComp cannot accommodate Position Matrices with a motif length greater than 32. Please filter out the matrix with this ID: %v.\n", motifs[i].Id)
		}
		var currSeq = ConsensusSequence(motifs[i], false)
		consensusScore, _, couldScoreConsensus = ScoreWindow(motifs[i], currSeq.Seq, 0)
		if !couldScoreConsensus {
			log.Fatalf("Error: problem in buildKmerHash. Could not score consensus sequence.")
		}

		altEndsConsidered := make(map[int]bool, 0)

		kmerHash = BuildKmerHash(motifs[i], s.PropMatch)
		scanRefSequenceComp(s.Records, kmerHash, motifs[i], s.ChromName, out, s.ResidualWindowSize, consensusScore, bed.Positive, s.RefStart, s.EnforceStrandMatch, s.OutputAsProportion, altEndsConsidered, s.ResidualFilter)

		revCompMotif = ReverseComplement(motifs[i])
		revKmerHash = BuildKmerHash(revCompMotif, s.PropMatch)
		scanRefSequenceComp(s.Records, revKmerHash, revCompMotif, s.ChromName, out, s.ResidualWindowSize, consensusScore, bed.Negative, s.RefStart, s.EnforceStrandMatch, s.OutputAsProportion, altEndsConsidered, s.ResidualFilter)

		//now we scan the alt sequences for any motifs lost in ref
		scanAltSequenceComp(s.Records, kmerHash, motifs[i], s.ChromName, out, s.ResidualWindowSize, consensusScore, bed.Positive, s.RefStart, s.EnforceStrandMatch, s.OutputAsProportion, altEndsConsidered, s.ResidualFilter)
		scanAltSequenceComp(s.Records, revKmerHash, revCompMotif, s.ChromName, out, s.ResidualWindowSize, consensusScore, bed.Negative, s.RefStart, s.EnforceStrandMatch, s.OutputAsProportion, altEndsConsidered, s.ResidualFilter)
	}

	err = out.Close()
	exception.PanicOnErr(err)
}

func scanRefSequenceComp(records []fasta.Fasta, kmerHash map[uint64]float64, pm PositionMatrix, chromName string, out *fileio.EasyWriter, residualWindowSize int, consensusScore float64, strand bed.Strand, refStart int, enforceStrandMatch bool, outputAsProportion bool, altEndsConsidered map[int]bool, residualFilter float64) {
	var needNewKey bool = true
	var currRefKey uint64
	var couldGetNewKey, inKmerHash, couldScoreSequence bool
	var currRefScore, currAltScore, currResidual, minResidual, minResidualAltScore float64
	var refPos, lastRefPos, lastAlnPos, currAltStart, currAltEnd int = 0, 0, 0, 0, 0 // original line: var refPos, lastRefPos, lastAlnPos, currAltStart, currAltEnd int = refStart, refStart, 0, 0, 0
	var bitMask uint64 = uint64(math.Pow(2, float64(2*len(pm.Mat[0]))) - 1)          // bitmask formula: B_n = 2^{2n} - 1
	var currBed bed.Bed
	revCompPm := ReverseComplement(pm)

	// scan through reference sequence first
	for alnPos := 0; alnPos < len(records[0].Seq); alnPos++ {
		refPos = fasta.AlnPosToRefPosCounter(records[0], alnPos, lastRefPos, lastAlnPos)
		//first we need a key for the alnPos. We can either make a new key or update our current key
		if needNewKey {
			currRefKey, alnPos, couldGetNewKey = getNewKey(records[0], alnPos, len(pm.Mat[0]))
			refPos = fasta.AlnPosToRefPosCounter(records[0], alnPos, lastRefPos, lastAlnPos)
			lastRefPos, lastAlnPos = refPos, alnPos
			if !couldGetNewKey {
				break //this means we've run out of windows on the reference sequence.
			}
			needNewKey = false
		} else {
			switch records[0].Seq[alnPos] {
			case dna.N:
				needNewKey = true
				continue
			case dna.Gap:
				continue
			case dna.A:
				currRefKey = currRefKey << 2
				currRefKey = currRefKey | 0
				currRefKey = currRefKey & bitMask
			case dna.C:
				currRefKey = currRefKey << 2
				currRefKey = currRefKey | 1
				currRefKey = currRefKey & bitMask
			case dna.G:
				currRefKey = currRefKey << 2
				currRefKey = currRefKey | 2
				currRefKey = currRefKey & bitMask
			case dna.T:
				currRefKey = currRefKey << 2
				currRefKey = currRefKey | 3
				currRefKey = currRefKey & bitMask
			default:
				log.Fatalf("Unrecognized base.")
			}
		}

		// now we check if the current key is a significant motif hit
		if currRefScore, inKmerHash = kmerHash[currRefKey]; inKmerHash {
			fmt.Printf("CurrRefKey: %v. CurrRefScore: %v.\n", currRefKey, currRefScore)
			minResidual = math.Inf(1)
			minResidualAltScore = 0
			fmt.Printf("Just set minResidual to inf, minResidualAltScore to 0\n")
			//TODO: fix this line so loop runs?
			// mirroring function line below
			//for currRefStart = numbers.Max(alnPos-len(pm.Mat[0])-residualWindowSize+1, 0); currRefStart <= numbers.Min(alnPos+residualWindowSize-len(pm.Mat[0])+1, len(records[0].Seq)); currRefStart++ {
			// yes len(pm.Mat[0]) correctly gets how many columns aka bases there are
			fmt.Printf("scanRefSequenceComp\n")
			fmt.Printf("for loop component checks. alnPos: %v, len(pm.Mat[0]): %v, residualWindowSize: %v, len(records[0].Seq): %v\n", alnPos, len(pm.Mat[0]), residualWindowSize, len(records[0].Seq))
			fmt.Printf("for loop endpoint checks. start: max(%v,0), end: min(%v, %v)\n", alnPos-len(pm.Mat[0])-residualWindowSize+1, alnPos+residualWindowSize-len(pm.Mat[0])+1, len(records[0].Seq))
			for currAltStart = numbers.Max(alnPos-len(pm.Mat[0])-residualWindowSize+1, 0); currAltStart <= numbers.Min(alnPos+residualWindowSize-len(pm.Mat[0])+1, len(records[0].Seq)); currAltStart++ {
				fmt.Printf("entered for loop\n")
				currAltScore, currAltEnd, couldScoreSequence = ScoreWindow(pm, records[1].Seq, currAltStart)
				//fmt.Printf("AlnPos: %v. currAltStart: %v. CurrAltEnd: %v. CurrAltScore: %v.\n", alnPos, currAltStart, currAltEnd, currAltScore)
				if !couldScoreSequence {
					fmt.Printf("!couldScoreSequence\n")
					break
				}
				currResidual = math.Abs(currRefScore - currAltScore)
				fmt.Printf("Just updated currResidual. currResidual: %v, currRefScore: %v, currAltScore: %v\n", currResidual, currRefScore, currAltScore)
				if currResidual < minResidual {
					minResidual = currResidual
					fmt.Printf("Just updated minResidual to = currResidual. minResidual: %v\n", minResidual)
					minResidualAltScore = currAltScore
					fmt.Printf("Just updated minResidualAltScore to = currAltScore. minResidualAltScore: %v\n", minResidualAltScore)
				}
				if !enforceStrandMatch {
					currAltScore, currAltEnd, couldScoreSequence = ScoreWindow(revCompPm, records[1].Seq, currAltStart)
					//fmt.Printf("AlnPos: %v. currAltStart: %v. CurrAltEnd: %v. CurrAltScore: %v.\n", alnPos, currAltStart, currAltEnd, currAltScore)
					if !couldScoreSequence {
						break
					}
					currResidual = math.Abs(currRefScore - currAltScore)
					if currResidual < minResidual {
						minResidual = currResidual
						fmt.Printf("Just updated minResidual to = currResidual. minResidual: %v\n", minResidual)
						minResidualAltScore = currAltScore
						fmt.Printf("Just updated minResidualAltScore to = currAltScore. minResidualAltScore: %v\n", minResidualAltScore)
					}
				}
				altEndsConsidered[currAltEnd] = true
			}
			minResidual = math.Abs(currRefScore - minResidualAltScore) // I think this line needs to exist after for loop, in case minResidual never gets updated in for loop, and also not using outputAsProportion?
			if outputAsProportion {
				currRefScore = currRefScore / consensusScore
				minResidualAltScore = minResidualAltScore / consensusScore
				minResidual = math.Abs(currRefScore - minResidualAltScore)
				fmt.Printf("Just updated minResidual and minResidualAltScore in outputAsProportion option\n")
			}
			if minResidual >= residualFilter {
				fmt.Printf("scanRefSequenceComp. About to write bed. ChromStart: %v, pm.Name: %v, currRefScore: %v, minResidualAltScore: %v, minResidual: %v\n", refPos-len(pm.Mat[0])+1, pm.Name, currRefScore, minResidualAltScore, minResidual)
				currBed = bed.Bed{
					Chrom:             chromName,
					ChromStart:        refStart + refPos - len(pm.Mat[0]) + 1, //original: refPos - len(pm.Mat[0]) + 1,
					ChromEnd:          refStart + refPos + 1,                  //original: refPos + 1,
					Name:              pm.Name,
					Score:             0,
					Strand:            strand,
					FieldsInitialized: 9,
					Annotation:        []string{fmt.Sprintf("%v", currRefScore), fmt.Sprintf("%v", minResidualAltScore), fmt.Sprintf("%v", minResidual)},
				}
				bed.WriteBed(out, currBed)
			}
		}
	}
}

func scanAltSequenceComp(records []fasta.Fasta, kmerHash map[uint64]float64, pm PositionMatrix, chromName string, out *fileio.EasyWriter, residualWindowSize int, consensusScore float64, strand bed.Strand, refStart int, enforceStrandMatch bool, outputAsProportion bool, altEndsConsidered map[int]bool, residualFilter float64) {
	var needNewKey = true
	var currAltKey uint64
	refPos := 0 //original line: var refPos = refStart
	var lastAlnPos, currRefStart int
	lastRefPos := 0 //original line: var lastRefPos = refStart
	var currAltScore, minResidual, minResidualRefScore, currRefScore, currResidual float64
	var couldGetNewKey, inKmerHash, foundInMap, couldScoreSequence bool
	var bitMask uint64 = uint64(math.Pow(2, float64(2*len(pm.Mat[0]))) - 1) // bitmask formula: B_n = 2^{2n} - 1
	var currBed bed.Bed
	revCompPm := ReverseComplement(pm)

	for alnPos := 0; alnPos < len(records[0].Seq); alnPos++ {
		refPos = fasta.AlnPosToRefPosCounter(records[0], alnPos, lastRefPos, lastAlnPos)
		if needNewKey {
			currAltKey, alnPos, couldGetNewKey = getNewKey(records[1], alnPos, len(pm.Mat[0]))
			refPos = fasta.AlnPosToRefPosCounter(records[0], alnPos, lastRefPos, lastAlnPos)
			lastRefPos, lastAlnPos = refPos, alnPos
			if !couldGetNewKey {
				break //this means we've run out of windows on the reference sequence.
			}
			needNewKey = false
		} else {
			switch records[1].Seq[alnPos] {
			case dna.N:
				needNewKey = true
				continue
			case dna.Gap:
				continue
			case dna.A:
				currAltKey = currAltKey << 2
				currAltKey = currAltKey | 0
				currAltKey = currAltKey & bitMask
			case dna.C:
				currAltKey = currAltKey << 2
				currAltKey = currAltKey | 1
				currAltKey = currAltKey & bitMask
			case dna.G:
				currAltKey = currAltKey << 2
				currAltKey = currAltKey | 2
				currAltKey = currAltKey & bitMask
			case dna.T:
				currAltKey = currAltKey << 2
				currAltKey = currAltKey | 3
				currAltKey = currAltKey & bitMask
			default:
				log.Fatalf("Unrecognized base.")
			}
		}

		// now we check if we have a significant key
		if currAltScore, inKmerHash = kmerHash[currAltKey]; inKmerHash {
			// now we have to check if we've already considered this hit when looking at the Ref sequence.
			if _, foundInMap = altEndsConsidered[alnPos]; !foundInMap {
				fmt.Printf("AlnPos: %v. RefPos: %v. currAltScore: %v.\n", alnPos, refPos, currAltScore/consensusScore)
				minResidual = math.Inf(1)
				minResidualRefScore = 0
				fmt.Printf("Just set minResidual to inf, and minResidualRefSCore to 0.\n")
				// TODO: maybe the line below is wrong because uses refStart, maybe should be diff variable. refStart should only be used when reporting final bed
				fmt.Printf("for loop component checks. alnPos: %v, len(pm.Mat[0]): %v, residualWindowSize: %v, len(records[0].Seq): %v\n", alnPos, len(pm.Mat[0]), residualWindowSize, len(records[0].Seq))
				fmt.Printf("for loop endpoint checks. start: max(%v,0), end: min(%v, %v)\n", alnPos-len(pm.Mat[0])-residualWindowSize+1, alnPos+residualWindowSize-len(pm.Mat[0])+1, len(records[0].Seq))
				for currRefStart = numbers.Max(alnPos-len(pm.Mat[0])-residualWindowSize+1, 0); currRefStart <= numbers.Min(alnPos+residualWindowSize-len(pm.Mat[0])+1, len(records[0].Seq)); currRefStart++ {
					// original line below
					//for currRefStart = numbers.Max(alnPos-len(pm.Mat[0])-residualWindowSize+1, refStart); currRefStart <= numbers.Min(alnPos+residualWindowSize-len(pm.Mat[0])+1, len(records[0].Seq)); currRefStart++ {
					fmt.Printf("entered for loop\n")
					currRefScore, _, couldScoreSequence = ScoreWindow(pm, records[0].Seq, currRefStart)
					if !couldScoreSequence {
						fmt.Printf("!couldScoreSequence\n")
						break
					}
					currResidual = math.Abs(currRefScore - currAltScore)
					fmt.Printf("Just updated currResidual\n")
					if currResidual < minResidual {
						minResidual = currResidual
						fmt.Printf("Just updated minResidual to = currResidual. minResidual: %v\n", minResidual)
						minResidualRefScore = currRefScore
						fmt.Printf("Just updated minResidualRefScore to = currRefScore. minResidualRefScore: %v\n", minResidualRefScore)
					}
					if !enforceStrandMatch {
						currRefScore, _, couldScoreSequence = ScoreWindow(revCompPm, records[0].Seq, currRefStart)
						//fmt.Printf("AlnPos: %v. currRefStart: %v. CurrRefEnd: %v. CurrRefScore: %v.\n", alnPos, currRefStart, currRefEnd, currRefScore)
						if !couldScoreSequence {
							break
						}
						currResidual = math.Abs(currRefScore - currAltScore)
						if currResidual < minResidual {
							minResidual = currResidual
							fmt.Printf("Just updated minResidual to = currResidual. minResidual: %v\n", minResidual)
							minResidualRefScore = currRefScore
							fmt.Printf("Just updated minResidualRefScore to = currRefScore. minResidualRefScore: %v\n", minResidualRefScore)
						}
					}
				}
				if outputAsProportion {
					currAltScore = currAltScore / consensusScore
					minResidualRefScore = minResidualRefScore / consensusScore
					minResidual = math.Abs(currAltScore - minResidualRefScore)
					fmt.Printf("Just updated minResidual, minResidualRefScore and currAltScore in outputAsProportion option\n")
				}
				if minResidual >= residualFilter {
					fmt.Printf("scanAltSequenceComp. About to write bed. ChromStart: %v, pm.Name: %v, minResidualRefScore: %v, currAltScore: %v, minResidual: %v\n", refPos-len(pm.Mat[0])+1, pm.Name, minResidualRefScore, currAltScore, minResidual)
					currBed = bed.Bed{
						Chrom:             chromName,
						ChromStart:        refStart + refPos - len(pm.Mat[0]) + 1, //original: refPos - len(pm.Mat[0]) + 1,
						ChromEnd:          refStart + refPos + 1,                  //original: refPos + 1,
						Name:              pm.Name,
						Score:             0,
						Strand:            strand,
						FieldsInitialized: 9,
						Annotation:        []string{fmt.Sprintf("%v", minResidualRefScore), fmt.Sprintf("%v", currAltScore), fmt.Sprintf("%v", minResidual)},
					}
					bed.WriteBed(out, currBed)
				}
			}
		}

	}
}

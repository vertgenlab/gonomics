// Command Group: "Data Conversion"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/sam"
	"log"
	"math/rand"
)

// Settings specifies all user options and arguments.
type Settings struct {
	SamFileName         string
	RefFile             string
	OutFileA            string
	OutFileB            string
	MultiFaDir          string
	qNameA              string
	qNameB              string
	Delta               float64
	Gamma               float64
	Epsilon             float64
	Kappa               float64
	LikelihoodCacheSize int
	SetSeed             int64
	Verbose             int
	FlatPrior           bool
}

// AnswerStruct contains all variables related to the two output assembled sequences, here labeled A and B.
// This includes the []fasta.Fasta sequences, their current positions, the remaining buffer room for the
// sequence slice, the index for the current chromosome string, and an empty slice of []dna.Base to resupply
// empty buffer room.
type AnswerStruct struct {
	AnswerA            []fasta.Fasta
	AnswerB            []fasta.Fasta
	AnswerAPos         int
	AnswerBPos         int
	EmptyRoomInBufferA int
	EmptyRoomInBufferB int
	CurrFaIndex        int
	NewBufferRoom      []dna.Base
}

// MultiFaStruct houses all variables related to the optional multiFa output, including
// the multiFa struct itself (CurrMultiFa), the current position, remaining room
// in the sequence buffer, and an empty slice of []dna.Base to resupply empty buffer room.
type MultiFaStruct struct {
	CurrMultiFa              []fasta.Fasta
	MultiFaPos               int
	EmptyRoomInMultiFaBuffer int
	NewBufferRoom            []dna.Base
}

const bufferSize = 10_000_000

func samAssembler(s Settings) {
	rand.Seed(s.SetSeed)
	var i, refPos, positionsToSkip, haploidBases int
	var haploidStrand bool // haploidStrand is true when the haploid bases are on the first strand (deletion on second strand)
	var currChrom string
	var newBufferRoom = make([]dna.Base, bufferSize)
	var firstTime = true
	var currPloidy = 2
	var diploidBaseCall sam.DiploidBase
	var currDiploidBases, currHaploidBases []dna.Base
	var currHaploidCall sam.HaploidCall
	var currRand float64

	// initialize caches for likelihood and priors
	cacheStruct := cacheSetup(s)

	// initialize reference genome
	ref := fasta.Read(s.RefFile)
	for i = range ref {
		dna.AllToUpper(ref[i].Seq)
	}
	refMap := fasta.ToMap(ref)

	// read pileups from sam/bam
	reads, header := sam.GoReadToChan(s.SamFileName)
	piles := sam.GoPileup(reads, header, false, nil, nil)

	// initialize output
	ans := AnswerStruct{
		AnswerA:       make([]fasta.Fasta, len(ref)),
		AnswerB:       make([]fasta.Fasta, len(ref)),
		NewBufferRoom: newBufferRoom,
	}
	for i = range ref {
		ans.AnswerA[i] = fasta.Fasta{Name: ref[i].Name, Seq: make([]dna.Base, bufferSize)}
		ans.AnswerB[i] = fasta.Fasta{Name: ref[i].Name, Seq: make([]dna.Base, bufferSize)}
	}
	var mlt MultiFaStruct

	// now time for the main loop, we look through each pile
	for p := range piles {
		if s.Verbose > 0 && !firstTime && p.Pos%1_000_000 == 0 {
			fmt.Printf("Current Chrom: %s. Current Position: %v.\n", currChrom, p.Pos)
		}
		if positionsToSkip > 0 {
			mlt = updateMultiFa(refMap[currChrom][refPos], dna.Gap, dna.Gap, mlt)
			refPos++
			positionsToSkip--
			continue
		}
		if firstTime {
			firstTime = false
			currChrom = header.Chroms[p.RefIdx].Name
			ans.CurrFaIndex = getIndexForName(ans.AnswerA, currChrom)
			ans.AnswerAPos, ans.AnswerBPos, refPos = 0, 0, 0
			ans.EmptyRoomInBufferA, ans.EmptyRoomInBufferB = bufferSize, bufferSize
			mlt = MultiFaStruct{
				CurrMultiFa: []fasta.Fasta{{Name: currChrom, Seq: make([]dna.Base, bufferSize)},
					{Name: s.qNameA, Seq: make([]dna.Base, bufferSize)},
					{Name: s.qNameB, Seq: make([]dna.Base, bufferSize)}},
				EmptyRoomInMultiFaBuffer: bufferSize,
				NewBufferRoom:            newBufferRoom,
				MultiFaPos:               0,
			}
		}
		if currChrom != header.Chroms[p.RefIdx].Name { // if we've moved onto a new chromosome.
			for refPos < len(refMap[currChrom]) { // write out the rest of the current reference
				ans.AnswerA[ans.CurrFaIndex].Seq[ans.AnswerAPos] = refMap[currChrom][refPos]
				ans.AnswerB[ans.CurrFaIndex].Seq[ans.AnswerBPos] = refMap[currChrom][refPos]
				ans = advanceAnswerPos(ans)
				mlt = updateMultiFa(refMap[currChrom][refPos], refMap[currChrom][refPos], refMap[currChrom][refPos], mlt)
				refPos++
			}
			// write out current chromosome multiFa
			ans = clearAnswerBuffer(ans)
			mlt = clearMultiFaBuffer(mlt)
			if s.MultiFaDir != "" {
				fasta.Write(fmt.Sprintf("%s/%s.fa", s.MultiFaDir, currChrom), mlt.CurrMultiFa)
			}
			// now we set up the new chromosome
			currChrom = header.Chroms[p.RefIdx].Name
			mlt.CurrMultiFa = []fasta.Fasta{{Name: currChrom, Seq: make([]dna.Base, bufferSize)},
				{Name: s.qNameA, Seq: make([]dna.Base, bufferSize)},
				{Name: s.qNameB, Seq: make([]dna.Base, bufferSize)}}
			mlt.EmptyRoomInMultiFaBuffer = bufferSize
			mlt.MultiFaPos = 0
			ans.CurrFaIndex = getIndexForName(ans.AnswerA, currChrom)
			ans.EmptyRoomInBufferA, ans.EmptyRoomInBufferB = bufferSize, bufferSize
			ans.AnswerAPos, ans.AnswerBPos, refPos = 0, 0, 0
		}

		// catch up to the current pile position, handles reference positions with no Pile coverage.
		for refPos < int(p.Pos-1) {
			ans.AnswerA[ans.CurrFaIndex].Seq[ans.AnswerAPos] = refMap[currChrom][refPos]
			ans.AnswerB[ans.CurrFaIndex].Seq[ans.AnswerBPos] = refMap[currChrom][refPos]
			mlt = updateMultiFa(refMap[currChrom][refPos], refMap[currChrom][refPos], refMap[currChrom][refPos], mlt)
			ans = advanceAnswerPos(ans)
			refPos++
		}

		// now refPos should equal p.Pos - 1, because of our for loop before
		if refPos != int(p.Pos-1) {
			log.Fatalf("Something went wrong. RefPos is not equal to p.Pos -1.")
		}

		if currPloidy == 2 {
			// First we handle the base call for the current pile
			diploidBaseCall = sam.DiploidBaseCallFromPile(p, refMap[currChrom][refPos], cacheStruct.DiploidBasePriorCache, cacheStruct.HomozygousBaseCache, cacheStruct.HeterozygousBaseCache, s.Epsilon)
			currDiploidBases = sam.DiploidBaseToBases(diploidBaseCall)
			currRand = rand.Float64()
			if currRand < 0.5 {
				ans.AnswerA[ans.CurrFaIndex].Seq[ans.AnswerAPos] = currDiploidBases[0]
				ans.AnswerB[ans.CurrFaIndex].Seq[ans.AnswerBPos] = currDiploidBases[1]
				mlt = updateMultiFa(refMap[currChrom][refPos], currDiploidBases[0], currDiploidBases[1], mlt)
			} else {
				ans.AnswerA[ans.CurrFaIndex].Seq[ans.AnswerAPos] = currDiploidBases[1]
				ans.AnswerB[ans.CurrFaIndex].Seq[ans.AnswerBPos] = currDiploidBases[0]
				mlt = updateMultiFa(refMap[currChrom][refPos], currDiploidBases[1], currDiploidBases[0], mlt)
			}
			ans = advanceAnswerPos(ans)

			// Now we handle diploid insertion calls in a helper function
			ans, mlt, cacheStruct, refPos = diploidInsertion(ans, mlt, cacheStruct, p, refPos, s)

			// Now we handle diploid deletion calls
			mlt, cacheStruct, refPos, haploidStrand, currPloidy, haploidBases, positionsToSkip = diploidDeletion(mlt, cacheStruct, p, refMap, refPos, currChrom, s)

		} else if currPloidy == 1 {
			currHaploidCall = sam.HaploidCallFromPile(p, refMap[currChrom][refPos], s.Epsilon, cacheStruct.HaploidBasePriorCache, cacheStruct.HaploidIndelPriorCache, cacheStruct.HomozygousBaseCache, cacheStruct.HeterozygousBaseCache, cacheStruct.HomozygousIndelCache)

			if haploidStrand {
				ans = advanceAPos(ans)
				mlt = updateMultiFa(refMap[currChrom][refPos], currHaploidCall.Base, dna.Gap, mlt)
				if currHaploidCall.Insertion != "" {
					currHaploidBases = dna.StringToBases(currHaploidCall.Insertion)
					for i = 0; i < len(currHaploidBases); i++ {
						ans = advanceAPos(ans)
						mlt = updateMultiFa(dna.Gap, currHaploidBases[i], dna.Gap, mlt)
					}
				}
				if currHaploidCall.Deletion != 0 {
					for i = 0; i < currHaploidCall.Deletion; i++ {
						mlt = updateMultiFa(refMap[currChrom][refPos], dna.Gap, dna.Gap, mlt)
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
				ans = advanceBPos(ans)
				mlt = updateMultiFa(refMap[currChrom][refPos], dna.Gap, currHaploidCall.Base, mlt)

				if currHaploidCall.Insertion != "" {
					currHaploidBases = dna.StringToBases(currHaploidCall.Insertion)
					for i = 0; i < len(currHaploidBases); i++ {
						ans = advanceBPos(ans)
						mlt = updateMultiFa(dna.Gap, dna.Gap, currHaploidBases[i], mlt)
					}
				}
				if currHaploidCall.Deletion != 0 {
					for i = 0; i < currHaploidCall.Deletion; i++ {
						mlt = updateMultiFa(refMap[currChrom][refPos], dna.Gap, dna.Gap, mlt)
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

			if haploidBases < 2 { // if we are on the last haploidBase, we re-enter diploid mode
				currPloidy = 2
			}
			refPos++
			haploidBases--
		} else {
			log.Fatalf("Error: Unrecognized ploidy: %v.\n", currPloidy)
		}
	}

	// once we're done with the piles we have to add the trailing ref bases and clear the buffer for the last chrom
	for refPos < len(refMap[currChrom]) {
		ans.AnswerA[ans.CurrFaIndex].Seq[ans.AnswerAPos] = refMap[currChrom][refPos]
		ans.AnswerB[ans.CurrFaIndex].Seq[ans.AnswerBPos] = refMap[currChrom][refPos]
		mlt = updateMultiFa(refMap[currChrom][refPos], refMap[currChrom][refPos], refMap[currChrom][refPos], mlt)
		mlt = advanceMultiFaPos(mlt)
		ans = advanceAnswerPos(ans)
		refPos++
	}

	// write last multiFa entry
	if s.MultiFaDir != "" {
		mlt = clearMultiFaBuffer(mlt)
		fasta.Write(fmt.Sprintf("%s/%s.fa", s.MultiFaDir, currChrom), mlt.CurrMultiFa)
	}

	// write answer
	ans = clearAnswerBuffer(ans)
	fasta.Write(s.OutFileA, ans.AnswerA)
	fasta.Write(s.OutFileB, ans.AnswerB)
}

// advanceAnswerPos advances the position of both haplotypes in the AnswerStruct.
func advanceAnswerPos(ans AnswerStruct) AnswerStruct {
	ans = advanceAPos(ans)
	return advanceBPos(ans)
}

// advanceAPos advances the position of haplotype A in the AnswerStruct, allocating additional buffer room if required.
func advanceAPos(ans AnswerStruct) AnswerStruct {
	ans.AnswerAPos++
	ans.EmptyRoomInBufferA--
	if ans.EmptyRoomInBufferA < 1 {
		ans.AnswerA[ans.CurrFaIndex].Seq = append(ans.AnswerA[ans.CurrFaIndex].Seq, ans.NewBufferRoom...)
		ans.EmptyRoomInBufferA += bufferSize
	}
	return ans
}

// advanceBPos advances the position of haplotype B in the AnswerStruct, allocating additional buffer room if required.
func advanceBPos(ans AnswerStruct) AnswerStruct {
	ans.AnswerBPos++
	ans.EmptyRoomInBufferB--
	if ans.EmptyRoomInBufferB < 1 {
		ans.AnswerB[ans.CurrFaIndex].Seq = append(ans.AnswerB[ans.CurrFaIndex].Seq, ans.NewBufferRoom...)
		ans.EmptyRoomInBufferB += bufferSize
	}
	return ans
}

// updateMultiFa assigns bases to the three sequences in the MultiFaStruct, and advances the multiFa position.
func updateMultiFa(zero dna.Base, first dna.Base, second dna.Base, mlt MultiFaStruct) MultiFaStruct {
	mlt.CurrMultiFa[0].Seq[mlt.MultiFaPos] = zero
	mlt.CurrMultiFa[1].Seq[mlt.MultiFaPos] = first
	mlt.CurrMultiFa[2].Seq[mlt.MultiFaPos] = second
	return advanceMultiFaPos(mlt)
}

// advanceMultiFaPos marches the position of the MultiFaStruct forward by one.
// We therefore adjust the remaining room in the buffer and allocate additional buffer room if required.
func advanceMultiFaPos(mlt MultiFaStruct) MultiFaStruct {
	mlt.MultiFaPos++
	mlt.EmptyRoomInMultiFaBuffer--
	if mlt.EmptyRoomInMultiFaBuffer < 1 {
		mlt.CurrMultiFa[0].Seq = append(mlt.CurrMultiFa[0].Seq, mlt.NewBufferRoom...)
		mlt.CurrMultiFa[1].Seq = append(mlt.CurrMultiFa[1].Seq, mlt.NewBufferRoom...)
		mlt.CurrMultiFa[2].Seq = append(mlt.CurrMultiFa[2].Seq, mlt.NewBufferRoom...)
		mlt.EmptyRoomInMultiFaBuffer += bufferSize
	}
	return mlt
}

// getIndexForName finds the index of a []fasta.Fasta with a specified input name.
func getIndexForName(f []fasta.Fasta, name string) int {
	for i := range f {
		if f[i].Name == name {
			return i
		}
	}
	log.Fatalf("Name: %s not found in fasta.", name)
	return -1
}

// clearMultiFaBuffer removes trailing bases from the end of the fasta seqeunces, based on the remaining room
// in the multiFa buffer.
func clearMultiFaBuffer(mlt MultiFaStruct) MultiFaStruct {
	mlt.CurrMultiFa[0].Seq = mlt.CurrMultiFa[0].Seq[:len(mlt.CurrMultiFa[0].Seq)-mlt.EmptyRoomInMultiFaBuffer]
	mlt.CurrMultiFa[1].Seq = mlt.CurrMultiFa[1].Seq[:len(mlt.CurrMultiFa[1].Seq)-mlt.EmptyRoomInMultiFaBuffer]
	mlt.CurrMultiFa[2].Seq = mlt.CurrMultiFa[2].Seq[:len(mlt.CurrMultiFa[2].Seq)-mlt.EmptyRoomInMultiFaBuffer]
	return mlt
}

// clearAnswerBuffer removes trailing bases at the end of the fasta sequence, based on the remaining room in each buffer.
func clearAnswerBuffer(ans AnswerStruct) AnswerStruct {
	ans.AnswerA[ans.CurrFaIndex].Seq = ans.AnswerA[ans.CurrFaIndex].Seq[:len(ans.AnswerA[ans.CurrFaIndex].Seq)-ans.EmptyRoomInBufferA]
	ans.AnswerB[ans.CurrFaIndex].Seq = ans.AnswerB[ans.CurrFaIndex].Seq[:len(ans.AnswerB[ans.CurrFaIndex].Seq)-ans.EmptyRoomInBufferB]
	return ans
}

func usage() {
	fmt.Print(
		"samAssembler - Reference-based diploid assembly of aligned short reads." +
			"Can be used in two modes: 'build' generates diploid assemblies." +
			"'score' validates assembly accuracy with a five way alignment including the known divergent sequences" +
			"In 'score', the user specifies a list of filenames, corresponding to each multiFa to score together (usually one file per chromosome)." +
			"For score, type may be chosen from the following: smallBaseMatrix, baseMatrixByRefBase\n" +
			"Usage:\n" +
			"samAssembler build individual.sam ref.fa outputA.fa outputB.fa\n" +
			"OR\n" +
			"samAssembler score scoreType inFileList outFile.txt\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var delta *float64 = flag.Float64("delta", 0.01, "Set the expected divergence frequency.")
	var gamma *float64 = flag.Float64("gamma", 3, "Set the expected transition bias.")
	var epsilon *float64 = flag.Float64("epsilon", 0.01, "Set the expected misclassification error rate.")
	var kappa *float64 = flag.Float64("kappa", 0.1, "Set the expected proportion of divergent sites that are INDELs.")
	var multiFaDir *string = flag.String("multiFaDir", "", "Output the reference and generated sequences as an aligned multiFa, each file by chrom.")
	var qNameA *string = flag.String("qNameA", "QueryA", "Set the qName for the first generated chromosome in the optional multiFa output.")
	var qNameB *string = flag.String("qNameB", "QueryB", "Set the qName for the second generated chromosome in the optional multiFa output.")
	var likelihoodCacheSize *int = flag.Int("likelihoodCacheSize", 100, "Set the maximum dimension of the likelihood caches. Should be slightly larger than highest expected pile depth.")
	var setSeed *int64 = flag.Int64("setSeed", -1, "Use a specific seed for the RNG.")
	var verbose *int = flag.Int("verbose", 0, "Set to 1 to enable additional debug prints.")
	var flatPrior *bool = flag.Bool("flatPrior", false, "Use a flat prior instead of the default informative prior distribution.")
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) < 1 {
		flag.Usage()
		log.Fatalf("Expecting at least one argument, received 0.")
	}

	var mode = flag.Arg(0)

	switch mode {
	case "build":
		var expectedNumArgs int = 5
		if len(flag.Args()) != expectedNumArgs {
			flag.Usage()
			log.Fatalf("Error: expecting %d arguments, but got %d\n",
				expectedNumArgs, len(flag.Args()))
		}

		inFile := flag.Arg(1)
		refFile := flag.Arg(2)
		outFileA := flag.Arg(3)
		outFileB := flag.Arg(4)

		s := Settings{
			SamFileName:         inFile,
			RefFile:             refFile,
			OutFileA:            outFileA,
			OutFileB:            outFileB,
			MultiFaDir:          *multiFaDir,
			qNameA:              *qNameA,
			qNameB:              *qNameB,
			Delta:               *delta,
			Gamma:               *gamma,
			Epsilon:             *epsilon,
			Kappa:               *kappa,
			LikelihoodCacheSize: *likelihoodCacheSize,
			SetSeed:             *setSeed,
			Verbose:             *verbose,
			FlatPrior:           *flatPrior,
		}

		samAssembler(s)
	case "score":
		var expectedNumArgs int = 4
		if len(flag.Args()) != expectedNumArgs {
			flag.Usage()
			log.Fatalf("Error: expecting %d arguments, but got %d\n",
				expectedNumArgs, len(flag.Args()))
		}
		scoreType := flag.Arg(1)
		inDir := flag.Arg(2)
		outFile := flag.Arg(3)
		samAssemblerScore(scoreType, inDir, outFile)
	default:
		log.Fatalf("Unknown mode. samAssembler can be run with the first argument as 'build' or 'score'.")
	}
}

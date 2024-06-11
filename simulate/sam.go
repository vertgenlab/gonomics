package simulate

import (
	"fmt"
	"log"
	"math/rand"
	"strings"

	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/sam"
)

// IlluminaPairedSam generates pairs of sam reads randomly distributed across the input DNA sequence.
// The inputs are the name of the input DNA sequence, the sequence itself, the number of read pairs to generate, the length of each read,
// the average fragment size, the standard deviation of fragment sizes, the error rate where a base in
// the read will not match the input DNA, a numbers.binomialAlias that is used to speed up calculations, and
// output file handles for sam, bam, and a bool denoting if bam (or sam) should be the output.
// Whichever handle (sam or bam) is not being used can be nil.
func IlluminaPairedSam(refName string, ref []dna.Base, numPairs, readLen, avgFragmentSize int, avgFragmentStdDev float64, flatErrorRate float64, ancientErrorRate float64, flatBinomialAlias numbers.BinomialAlias, ancientBinomialAlias numbers.BinomialAlias, geometricParam float64, out *fileio.EasyWriter, bw *sam.BamWriter, bamOutput bool, deaminationDistributionSlice []int, rng *rand.Rand) {
	var fragmentSize, midpoint, startFor, endRev, newCapacity int
	var currFor, currRev sam.Sam
	var fragment []dna.Base = make([]dna.Base, 0, avgFragmentSize+int(5*avgFragmentStdDev)) //set initial fragment slice capacity to 5 deviations above average fragment size
	var newFragment []dna.Base                                                              //this will be used if we run out of capacity in 'fragment'
	if avgFragmentSize < readLen {
		log.Fatalf("Error: average fragment size %v is less than read length %v.\n", avgFragmentSize, readLen)
	}

	for i := 0; i < numPairs; i++ {
		fragmentSize = numbers.Max(readLen, int(numbers.SampleInverseNormal(float64(avgFragmentSize), avgFragmentStdDev)))
		midpoint = numbers.RandIntInRange(0, len(ref))
		startFor = numbers.Max(midpoint-(fragmentSize/2), 0)
		endRev = numbers.Min(midpoint+(fragmentSize/2), len(ref))

		if fragmentSize < readLen {
			readLen = fragmentSize
		}

		//here we check if fragment has sufficient capacity for the next fragment we drew. If not, we extend capacity.
		if len(fragment)+fragmentSize > cap(fragment) {
			newCapacity = len(fragment) + fragmentSize
			newFragment = make([]dna.Base, fragmentSize, newCapacity)
			fragment = newFragment
		} else {
			fragment = fragment[:fragmentSize]
		}
		copy(fragment, ref[startFor:endRev])

		if endRev > len(ref) || startFor < 0 {
			currFor, currRev = generateUnmapped(readLen)
		} else {
			if ancientErrorRate > 0 {
				ancientDamage(fragment, ancientBinomialAlias, geometricParam, deaminationDistributionSlice)
			}
			currFor, currRev = generateSamReadNoFlag(fmt.Sprintf("%s_Read:%d", refName, i), refName, fragment, readLen, startFor, flatErrorRate, flatBinomialAlias, rng)
		}

		if currFor.Cigar == nil && currRev.Cigar == nil {
			i -= 1 // retry
			continue
		}
		addPairedFlags(&currFor, &currRev)
		if currFor.Cigar != nil && currRev.Cigar != nil {
			currFor.RNext = "="
			currRev.RNext = "="
		} else {
			currFor.RNext = currRev.RName
			currRev.RNext = currFor.RName
		}

		currFor.PNext = currRev.Pos
		currRev.PNext = currFor.Pos
		if bamOutput {
			sam.WriteToBamFileHandle(bw, currFor, 0)
			sam.WriteToBamFileHandle(bw, currRev, 0)
		} else {
			sam.WriteToFileHandle(out, currFor)
			sam.WriteToFileHandle(out, currRev)
		}
	}
}

// generateUnmapped creates a paired end read with a random sequence that is unmapped.
func generateUnmapped(readLen int) (sam.Sam, sam.Sam) {
	var forRead, revRead sam.Sam
	forRead.RName = "*"
	revRead.RName = "*"
	forRead.Seq = make([]dna.Base, readLen)
	revRead.Seq = make([]dna.Base, readLen)
	for i := range forRead.Seq {
		forRead.Seq[i] = dna.Base(numbers.RandIntInRange(0, 4))
		revRead.Seq[i] = dna.Base(numbers.RandIntInRange(0, 4))
	}
	return forRead, revRead
}

// generateSamReadNoFlag generates a sam record for the input position.
// the second return is a deaminationDistributionSlice, recording the positions of observed cytosine deamination events.
// Soft clips sequence that is off template and does not generate Flag, RNext, or PNext.
func generateSamReadNoFlag(readName string, refName string, fragment []dna.Base, readLength int, fragmentStart int, flatErrorRate float64, flatAlias numbers.BinomialAlias, src rand.Source) (sam.Sam, sam.Sam) {
	var currForSam = sam.Sam{
		QName: readName,
		Seq:   make([]dna.Base, readLength),
	}
	var currRevSam = sam.Sam{
		QName: readName,
		Seq:   make([]dna.Base, readLength),
	}

	// generate quality scores for forward and reverse read
	var bldrFor strings.Builder
	for range currForSam.Seq {
		bldrFor.WriteRune(rune(numbers.RandIntInRange(30, 40) + 33)) // high quality seq + ascii offset
	}
	currForSam.Qual = bldrFor.String()
	var bldrRev strings.Builder
	for range currRevSam.Seq {
		bldrRev.WriteRune(rune(numbers.RandIntInRange(30, 40) + 33)) // high quality seq + ascii offset
	}
	currRevSam.Qual = bldrRev.String()

	currForSam.MapQ = uint8(numbers.RandIntInRangeSrc(30, 40, src))
	currRevSam.MapQ = uint8(numbers.RandIntInRangeSrc(30, 40, src))
	currForSam.RName, currRevSam.RName = refName, refName

	copy(currForSam.Seq, fragment[0:readLength])
	copy(currRevSam.Seq, fragment[len(fragment)-readLength:])

	// now we will simulate sequencing/PCR error with a flat error rate
	if flatErrorRate > 0 {
		currForSam = sequencingError(currForSam, flatAlias, src)
		currRevSam = sequencingError(currRevSam, flatAlias, src)
	}

	// generate other values
	currForSam.Pos = uint32(fragmentStart) + 1
	currRevSam.Pos = uint32(fragmentStart+len(fragment)-readLength) + 1
	currForSam.TLen = int32(readLength)
	currRevSam.TLen = int32(readLength)

	// assemble cigar
	currForSam.Cigar = append(currForSam.Cigar, cigar.Cigar{RunLength: readLength, Op: 'M'})
	currRevSam.Cigar = append(currRevSam.Cigar, cigar.Cigar{RunLength: readLength, Op: 'M'})

	return currForSam, currRevSam
}

// addPairedFlags adds the flag for a pair of sam records.
func addPairedFlags(f, r *sam.Sam) {
	var fIsRevComp bool = rand.Float64() > 0.5
	if fIsRevComp {
		*f, *r = *r, *f // so that the reads always point towards one another
	}
	f.Flag += 1 + 64
	r.Flag += 1 + 128
	switch {
	case f.Cigar != nil && r.Cigar != nil: // both mapped
		f.Flag += 2
		r.Flag += 2
		if fIsRevComp {
			f.Flag += 16
			r.Flag += 32
		} else {
			f.Flag += 32
			r.Flag += 16
		}

	case f.Cigar == nil && r.Cigar == nil: // both unmapped
		f.Flag += 4 + 8
		r.Flag += 4 + 8

	case f.Cigar != nil && r.Cigar == nil: // f mapped r unmapped
		f.Flag += 8
		r.Flag += 4
		if fIsRevComp {
			f.Flag += 16
			r.Flag += 32
		}

	case f.Cigar == nil && r.Cigar != nil: // f unmapped r mapped
		f.Flag += 4
		r.Flag += 8
		if !fIsRevComp {
			f.Flag += 32
			r.Flag += 16
		}
	}
}

// sequencingError takes in a sam.Sam record and a BinomialAlias to
// edit bases randomly across the read to simulate PCR/sequencing error. Errors are made with no
// positional dependence and with a flat error spectrum.
func sequencingError(currSam sam.Sam, alias numbers.BinomialAlias, src rand.Source) sam.Sam {
	numFlatErrors := numbers.RandBinomial(alias)         // sample a binomial distribution to get the number of sequencing errors
	mutatedPositions := make(map[int]int, numFlatErrors) // store positions we've mutated so that we can sample without replacement
	var foundInMap bool
	var currRandInt int
	var currError int = 0

	seed := rand.New(src)

	for currError < numFlatErrors {
		currRandInt = numbers.RandIntInRange(0, len(currSam.Seq)) // sample a base on the read
		if _, foundInMap = mutatedPositions[currRandInt]; !foundInMap {
			mutatedPositions[currRandInt] = 1
			currSam.Seq[currRandInt] = changeBase(currSam.Seq[currRandInt], seed)
			currError++
		}
	}
	return currSam
}

// ancientDamage introduces cytosine deamination events into a DNA fragment ([]dna.Base).
// ancientAlias is a numbers.BinomialAlias which allows rapid generation of binomial-distributed random variates.
// cytosine deamination events are distributed geometrically from fragment ends, based on an input geometric parameter.
// deaminationDistributionSlice records the locations of cytosine deamination events for debugging.
func ancientDamage(currFrag []dna.Base, ancientAlias numbers.BinomialAlias, geometricParam float64, deaminationDistributionSlice []int) {
	var currRandPos int
	var distanceToEnd int
	var foundInMap bool
	var whichSideRand float64
	numAncientErrorAttempts := numbers.RandBinomial(ancientAlias)
	damagedPositions := make(map[int]int, numAncientErrorAttempts)
	var currError int = 0
	for currError < numAncientErrorAttempts {
		distanceToEnd = numbers.RandGeometric(geometricParam)
		for distanceToEnd >= len(currFrag) {
			distanceToEnd = numbers.RandGeometric(geometricParam)
		}
		whichSideRand = rand.Float64() // we call a random number here to decide if the deamination attempt occurs on the left or right of the fragment
		if whichSideRand < 0.5 {
			currRandPos = len(currFrag) - distanceToEnd - 1 // distance from end, instead of distance from start
		} else {
			currRandPos = distanceToEnd // distance from start of read
		}
		if _, foundInMap = damagedPositions[currRandPos]; !foundInMap {
			damagedPositions[currRandPos] = 1
			switch currFrag[currRandPos] {
			case dna.A:
				// nothing to do
			case dna.C:
				currFrag[currRandPos] = dna.T
				if distanceToEnd < len(deaminationDistributionSlice) {
					deaminationDistributionSlice[distanceToEnd]++
				}
			case dna.G:
				currFrag[currRandPos] = dna.A
				if distanceToEnd < len(deaminationDistributionSlice) {
					deaminationDistributionSlice[distanceToEnd]++
				}
			case dna.T:
				// nothing to do
			default:
				log.Fatalf("Error: Unrecognized base: %v.\n", currFrag[currRandPos])
			}
			currError++
		}
	}
}

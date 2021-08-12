// Command Group: "Statistics & Population Genetics"

package main

import (
	"flag"
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

type Settings struct {
	InFile                     string
	ChromName                  string
	VelOut                     string
	AccelOut                   string
	InitialVelOut              string
	SearchSpaceBed             string
	SearchSpaceProportion      float64
	WindowSize                 int
	UseSnpDistance             bool
	Verbose                    bool
	Epsilon                    float64
	AllowNegative              bool
	ZeroDistanceWeightConstant float64
}

type BranchCache struct {
	ChromStart int
	ChromEnd   int
	b1         float64
	b3         float64
}

type Distances struct {
	D01 float64
	D02 float64
	D03 float64
	D12 float64
	D13 float64
	D23 float64
}

type BranchLengths struct {
	B1 float64
	B2 float64
	B3 float64
	B4 float64
	B5 float64
}

type SubTree struct {
	Dab float64
	Dac float64
	Dbc float64
	va  float64
	vb  float64
	vc  float64
}

func multiFaAcceleration(s Settings) {
	records := fasta.Read(s.InFile)
	var searchSpace []bed.Bed
	var referenceLength, i, j, threshold int
	var bitArray []byte

	//a couple of safety checks
	if len(records) != 4 {
		log.Fatalf("multiFaAcceleration accepts a multiFa file with 4 records, found %v.", len(records))
	}
	if len(records[1].Seq) != len(records[0].Seq) || len(records[2].Seq) != len(records[0].Seq) || len(records[3].Seq) != len(records[0].Seq) {
		log.Fatalf("Error. All records must be of the same sequence length.")
	}

	referenceLength = fasta.AlnPosToRefPos(records[0], len(records[0].Seq)-1)

	if s.SearchSpaceBed != "" {
		searchSpace = bed.Read(s.SearchSpaceBed)
		bitArray = make([]byte, referenceLength)
		for i = range searchSpace {
			if searchSpace[i].Chrom == s.ChromName {
				for j = searchSpace[i].ChromStart; j < searchSpace[i].ChromEnd; j++ {
					bitArray[j] = 1
				}
			}
		}
		threshold = int(s.SearchSpaceProportion * float64(s.WindowSize))
	}

	var currDistances Distances
	var distanceCache = make(map[Distances]BranchLengths)
	var referenceCounter int = 0
	var reachedEnd bool = false
	var b1, b3 float64
	var currCount int
	var pass, containedInMap bool

	//variables for normalization
	var velSum, initialSum float64 = 0, 0
	var branchCacheSlice = make([]BranchCache, 0)
	for alignmentCounter := 0; reachedEnd == false && referenceCounter < referenceLength-s.WindowSize; alignmentCounter++ {
		if s.Verbose && alignmentCounter%1000000 == 0 {
			log.Printf("alignmentCounter: %v\n", alignmentCounter)
		}
		currCount, pass = thresholdCheckPasses(s, currCount, threshold, bitArray, referenceCounter)
		if records[0].Seq[alignmentCounter] != dna.Gap {
			if pass {
				if s.UseSnpDistance {
					reachedEnd = fourWaySnpDistances(records, alignmentCounter, s, &currDistances)
				} else {
					reachedEnd = fourWayMutationDistances(records, alignmentCounter, s, &currDistances)
				}

				if _, containedInMap = distanceCache[currDistances]; !containedInMap { //if this tree has not been seen before, calculate branch lengths
					distanceCache[currDistances] = alternatingLeastSquares(currDistances, s)
				}

				b1 = distanceCache[currDistances].B1
				b3 = distanceCache[currDistances].B3

				if !reachedEnd {
					velSum += b1
					initialSum += b3

					branchCacheSlice = append(branchCacheSlice, BranchCache{referenceCounter, referenceCounter + s.WindowSize, b1, b3})
				}
			}
			referenceCounter++
		}
	}

	var averageVel, averageInitialVel, b1Normal, b3Normal float64
	averageVel = velSum / float64(len(branchCacheSlice))
	averageInitialVel = initialSum / float64(len(branchCacheSlice))

	var err error
	velBed := fileio.EasyCreate(s.VelOut)
	accelBed := fileio.EasyCreate(s.AccelOut)
	initialVelBed := fileio.EasyCreate(s.InitialVelOut)

	for i = range branchCacheSlice {
		b1Normal = branchCacheSlice[i].b1 / averageVel
		b3Normal = branchCacheSlice[i].b3 / averageInitialVel
		bed.WriteBed(velBed, bed.Bed{Chrom: s.ChromName, ChromStart: branchCacheSlice[i].ChromStart, ChromEnd: branchCacheSlice[i].ChromEnd, Name: fmt.Sprintf("%e", b1Normal), FieldsInitialized: 4})
		bed.WriteBed(initialVelBed, bed.Bed{Chrom: s.ChromName, ChromStart: branchCacheSlice[i].ChromStart, ChromEnd: branchCacheSlice[i].ChromEnd, Name: fmt.Sprintf("%e", b3Normal), FieldsInitialized: 4})
		bed.WriteBed(accelBed, bed.Bed{Chrom: s.ChromName, ChromStart: branchCacheSlice[i].ChromStart, ChromEnd: branchCacheSlice[i].ChromEnd, Name: fmt.Sprintf("%e", b1Normal-b3Normal), FieldsInitialized: 4})
	}

	err = velBed.Close()
	exception.PanicOnErr(err)
	err = accelBed.Close()
	exception.PanicOnErr(err)
	err = initialVelBed.Close()
	exception.PanicOnErr(err)
}

func alternatingLeastSquares(d Distances, s Settings) BranchLengths {
	var answer = BranchLengths{1, 1, 1, 1, 1}
	var Q float64 = calculateQ(d, answer, s)
	var nextQ float64
	var currDiff float64 = s.Epsilon + 1 //set currDiff to something larger than epsilon so that we make it into the loop the first time.
	var sub SubTree
	var maxIteration, i = 1000, 0

	for currDiff > s.Epsilon && i < maxIteration {
		pruneLeft(d, answer, &sub, s)
		answer.B1, answer.B2, answer.B3 = optimizeSubtree(&sub, s)
		pruneRight(d, answer, &sub, s)
		answer.B4, answer.B5, answer.B3 = optimizeSubtree(&sub, s)
		nextQ = calculateQ(d, answer, s)
		currDiff = math.Abs(nextQ - Q)
		Q = nextQ
		i++
	}
	if i >= maxIteration {
		log.Fatalf("Failed to converge.")
	}
	return answer
}

func optimizeSubtree(sub *SubTree, s Settings) (float64, float64, float64) {
	sub.va = (sub.Dab + sub.Dac - sub.Dbc) / 2.0
	sub.vb = (sub.Dab + sub.Dbc - sub.Dac) / 2.0
	sub.vc = (sub.Dac + sub.Dbc - sub.Dac) / 2.0

	if s.AllowNegative {
		return sub.va, sub.vb, sub.vc
	}
	if sub.va < 0 && sub.vb < 0 && sub.vc < 0 {
		if s.Verbose {
			log.Printf("WARNING: All branches are negative.")
		}
		sub.va, sub.vb, sub.vc = 0, 0, 0
	} else if sub.va < 0 && sub.vb < 0 {
		sub.va = 0
		sub.vb = 0
		sub.vc = nonNegativeApproximation(sub.Dac, sub.Dbc, sub.va, sub.vb, s)
	} else if sub.va < 0 && sub.vc < 0 {
		sub.va = 0
		sub.vc = 0
		sub.vb = nonNegativeApproximation(sub.Dbc, sub.Dab, sub.vc, sub.va, s)
	} else if sub.vb < 0 && sub.vc < 0 {
		sub.vb = 0
		sub.vc = 0
		sub.va = nonNegativeApproximation(sub.Dab, sub.Dac, sub.vb, sub.vc, s)
	} else if sub.va < 0 {
		sub.va = 0
		sub.vb = nonNegativeApproximation(sub.Dab, sub.Dbc, sub.va, sub.vc, s)
		sub.vc = nonNegativeApproximation(sub.Dac, sub.Dbc, sub.va, sub.vb, s)
	} else if sub.vb < 0 {
		sub.vb = 0
		sub.va = nonNegativeApproximation(sub.Dab, sub.Dac, sub.vb, sub.vc, s)
		sub.vc = nonNegativeApproximation(sub.Dab, sub.Dbc, sub.va, sub.vc, s)
	} else if sub.vc < 0 {
		sub.vc = 0
		sub.va = nonNegativeApproximation(sub.Dab, sub.Dac, sub.vb, sub.vc, s)
		sub.vb = nonNegativeApproximation(sub.Dab, sub.Dbc, sub.va, sub.vc, s)
	}
	return sub.va, sub.vb, sub.vc
}

func nonNegativeApproximation(d1 float64, d2 float64, v1 float64, v2 float64, s Settings) float64 {
	if d1 == 0 {
		if d2 == 0 {
			return numbers.MaxFloat64(0, (s.ZeroDistanceWeightConstant*(d1-v1)+s.ZeroDistanceWeightConstant*(d2-v2))/(2*s.ZeroDistanceWeightConstant))
		} else {
			return numbers.MaxFloat64(0, s.ZeroDistanceWeightConstant*(d1-v1)+(1.0/math.Pow(d2, 2)*(d2-v2))) / (s.ZeroDistanceWeightConstant + (1.0 / math.Pow(d2, 2)))
		}
	} else if d2 == 0 {
		return numbers.MaxFloat64(0, (1.0/(math.Pow(d1, 2))*(d1-v1)+s.ZeroDistanceWeightConstant*(d2-v2))/((1.0/math.Pow(d1, 2))+s.ZeroDistanceWeightConstant))
	}
	return numbers.MaxFloat64(0, (1.0/(math.Pow(d1, 2))*(d1-v1)+(1.0/math.Pow(d2, 2)*(d2-v2)))/((1.0/math.Pow(d1, 2))+(1.0/math.Pow(d2, 2))))
}

func pruneLeft(d Distances, b BranchLengths, sub *SubTree, s Settings) {
	sub.Dab = d.D01
	if d.D03 == 0 {
		if d.D02 == 0 {
			sub.Dac = (s.ZeroDistanceWeightConstant*(d.D02-b.B4) + s.ZeroDistanceWeightConstant*(d.D03-b.B5)) / (2 * s.ZeroDistanceWeightConstant)
		} else {
			sub.Dac = ((1.0/math.Pow(d.D02, 2))*(d.D02-b.B4) + (s.ZeroDistanceWeightConstant * (d.D03 - b.B5))) / (s.ZeroDistanceWeightConstant + (1.0 / math.Pow(d.D02, 2)))
		}
	} else if d.D02 == 0 {
		sub.Dac = (s.ZeroDistanceWeightConstant*(d.D02-b.B4) + ((1.0 / math.Pow(d.D03, 2)) * (d.D03 - b.B5))) / ((1.0 / math.Pow(d.D03, 2)) + s.ZeroDistanceWeightConstant)
	} else {
		sub.Dac = ((1.0/math.Pow(d.D02, 2))*(d.D02-b.B4) + (1.0/math.Pow(d.D03, 2))*(d.D03-b.B5)) / ((1.0 / math.Pow(d.D03, 2)) + (1.0 / math.Pow(d.D02, 2)))
	}

	if d.D13 == 0 {
		if d.D12 == 0 {
			sub.Dbc = (s.ZeroDistanceWeightConstant*(d.D12-b.B4) + s.ZeroDistanceWeightConstant*(d.D13-b.B5)) / (2 * s.ZeroDistanceWeightConstant)
		} else {
			sub.Dbc = (s.ZeroDistanceWeightConstant*(d.D13-b.B5) + (1.0/math.Pow(d.D12, 2))*(d.D12-b.B4)) / ((1.0 / math.Pow(d.D12, 2)) + s.ZeroDistanceWeightConstant)
		}
	} else if d.D12 == 0 {
		sub.Dbc = (s.ZeroDistanceWeightConstant*(d.D12-b.B4) + (1.0/math.Pow(d.D13, 2))*(d.D13-b.B5)) / ((1.0 / math.Pow(d.D13, 2)) + s.ZeroDistanceWeightConstant)
	} else {
		sub.Dbc = ((1.0/math.Pow(d.D12, 2))*(d.D12-b.B4) + (1.0/math.Pow(d.D13, 2))*(d.D13-b.B5)) / ((1.0 / math.Pow(d.D13, 2)) + (1.0 / math.Pow(d.D12, 2)))
	}
}

func pruneRight(d Distances, b BranchLengths, sub *SubTree, s Settings) {
	sub.Dac = d.D23

	if d.D02 == 0 {
		if d.D12 == 0 {
			sub.Dac = (s.ZeroDistanceWeightConstant*(d.D02-b.B1) + s.ZeroDistanceWeightConstant*(d.D12-b.B2)) / (2.0 * s.ZeroDistanceWeightConstant)
		} else {
			sub.Dac = ((1.0/math.Pow(d.D12, 2))*(d.D12-b.B2) + s.ZeroDistanceWeightConstant*(d.D02-b.B1)) / ((1.0 / math.Pow(d.D12, 2)) + s.ZeroDistanceWeightConstant)
		}
	} else if d.D12 == 0 {
		sub.Dac = ((1.0/math.Pow(d.D02, 2))*(d.D02-b.B1) + s.ZeroDistanceWeightConstant*(d.D12-b.B2)) / ((1.0 / math.Pow(d.D02, 2)) + s.ZeroDistanceWeightConstant)
	} else {
		sub.Dab = ((1.0/math.Pow(d.D02, 2))*(d.D02-b.B1) + (1.0/math.Pow(d.D12, 2))*(d.D12-b.B2)) / ((1.0 / math.Pow(d.D02, 2)) + (1.0 / math.Pow(d.D12, 2)))
	}

	if d.D03 == 0 {
		if d.D13 == 0 {
			sub.Dbc = (s.ZeroDistanceWeightConstant*(d.D03-b.B1) + s.ZeroDistanceWeightConstant*(d.D13-b.B2)) / (2.0 * s.ZeroDistanceWeightConstant)
		} else {
			sub.Dbc = ((1.0/math.Pow(d.D13, 2))*(d.D13-b.B2) + s.ZeroDistanceWeightConstant*(d.D03-b.B1)) / ((1.0 / math.Pow(d.D13, 2)) + s.ZeroDistanceWeightConstant)
		}
	} else if d.D13 == 0 {
		sub.Dbc = ((1.0/math.Pow(d.D03, 2))*(d.D03-b.B1) + s.ZeroDistanceWeightConstant*(d.D13-b.B2)) / ((1.0 / math.Pow(d.D03, 2)) + s.ZeroDistanceWeightConstant)
	} else {
		sub.Dbc = ((1.0/math.Pow(d.D03, 2))*(d.D03-b.B1) + (1.0/math.Pow(d.D13, 2))*(d.D13-b.B2)) / ((1.0 / math.Pow(d.D03, 2)) + (1.0 / math.Pow(d.D13, 2)))
	}
}

func calculateQ(d Distances, b BranchLengths, s Settings) float64 {
	var sum float64 = 0
	if d.D01 != 0 { //avoid divide by zero error
		sum += math.Pow(d.D01-b.B1-b.B2, 2) / math.Pow(d.D01, 2)
	} else {
		sum += math.Pow(d.D01-b.B1-b.B2, 2) * s.ZeroDistanceWeightConstant
	}
	if d.D02 != 0 {
		sum += math.Pow(d.D02-b.B1-b.B3-b.B4, 2) / math.Pow(d.D02, 2)
	} else {
		sum += math.Pow(d.D02-b.B1-b.B3-b.B4, 2) * s.ZeroDistanceWeightConstant
	}
	if d.D03 != 0 {
		sum += math.Pow(d.D03-b.B1-b.B3-b.B5, 2) / math.Pow(d.D03, 2)
	} else {
		sum += math.Pow(d.D03-b.B1-b.B3-b.B5, 2) * s.ZeroDistanceWeightConstant
	}
	if d.D12 != 0 {
		sum += math.Pow(d.D12-b.B2-b.B3-b.B4, 2) / math.Pow(d.D12, 2)
	} else {
		sum += math.Pow(d.D12-b.B2-b.B3-b.B4, 2) * s.ZeroDistanceWeightConstant
	}
	if d.D13 != 0 {
		sum += math.Pow(d.D13-b.B2-b.B3-b.B5, 2) / math.Pow(d.D13, 2)
	} else {
		sum += math.Pow(d.D13-b.B2-b.B3-b.B5, 2) * s.ZeroDistanceWeightConstant
	}
	if d.D23 != 0 {
		sum += math.Pow(d.D23-b.B4-b.B5, 2) / math.Pow(d.D23, 2)
	} else {
		sum += math.Pow(d.D23-b.B4-b.B5, 2) * s.ZeroDistanceWeightConstant
	}
	return sum
}

//bitArray is on reference coordinates, not alignment coordinates, so the window is simply equal to windowSize.
func thresholdCheckPasses(s Settings, currCount int, threshold int, bitArray []byte, referenceCounter int) (int, bool) {
	if s.SearchSpaceBed == "" { //no search space file, no need to look further
		return 0, true
	}
	if referenceCounter == 0 {
		currCount = 0
		for i := 0; i < s.WindowSize; i++ {
			currCount += int(bitArray[i])
		}
	} else {
		if bitArray[referenceCounter-1] == 1 {
			currCount--
		}
		if bitArray[referenceCounter+s.WindowSize-1] > 0 {
			currCount++
		}
	}
	return currCount, currCount >= threshold
}

func fourWayMutationDistances(records []fasta.Fasta, alignmentCounter int, s Settings, d *Distances) bool {
	//first we clear the values in d.
	d.D01, d.D02, d.D03, d.D12, d.D13, d.D23 = 0, 0, 0, 0, 0, 0
	var d01tmp int
	var reachedEnd bool
	var alnEnd int
	d01tmp, reachedEnd, alnEnd = fasta.PairwiseMutationDistanceReferenceWindow(records[0], records[1], alignmentCounter, s.WindowSize)
	d.D01 = float64(d01tmp)
	d.D02 = float64(fasta.PairwiseMutationDistanceInRange(records[0], records[2], alignmentCounter, alnEnd))
	d.D03 = float64(fasta.PairwiseMutationDistanceInRange(records[0], records[3], alignmentCounter, alnEnd))
	d.D12 = float64(fasta.PairwiseMutationDistanceInRange(records[1], records[2], alignmentCounter, alnEnd))
	d.D13 = float64(fasta.PairwiseMutationDistanceInRange(records[1], records[3], alignmentCounter, alnEnd))
	d.D23 = float64(fasta.PairwiseMutationDistanceInRange(records[2], records[3], alignmentCounter, alnEnd))
	return reachedEnd
}

func fourWaySnpDistances(records []fasta.Fasta, alignmentCounter int, s Settings, d *Distances) bool {
	//first we clear the values in d.
	d.D01, d.D02, d.D03, d.D12, d.D13, d.D23 = 0, 0, 0, 0, 0, 0
	var baseCount, i int = 0, 0
	var reachedEnd bool = false

	if len(records) != 4 {
		log.Fatalf("multiFaAcceleration must take in a four-way multiple alignment.")
	}
	for i = alignmentCounter; baseCount < s.WindowSize && i < len(records[0].Seq); i++ {
		if records[0].Seq[i] != dna.Gap {
			baseCount++
		}
		if isUngappedColumn(records, i) {
			if records[0].Seq[i] != records[1].Seq[i] {
				d.D01++
			}
			if records[0].Seq[i] != records[2].Seq[i] {
				d.D02++
			}
			if records[0].Seq[i] != records[3].Seq[i] {
				d.D03++
			}
			if records[1].Seq[i] != records[2].Seq[i] {
				d.D12++
			}
			if records[1].Seq[i] != records[3].Seq[i] {
				d.D13++
			}
			if records[2].Seq[i] != records[3].Seq[i] {
				d.D23++
			}
		}
	}
	if baseCount != s.WindowSize {
		reachedEnd = true
	}
	return reachedEnd
}

func isUngappedColumn(records []fasta.Fasta, index int) bool {
	for i := range records {
		if !isUngappedBase(records[i].Seq[index]) {
			return false
		}
	}
	return true
}

func isUngappedBase(b dna.Base) bool {
	if b == dna.A || b == dna.T || b == dna.C || b == dna.G {
		return true
	}
	if b == dna.LowerA || b == dna.LowerC || b == dna.LowerG || b == dna.LowerT {
		return true
	}
	return false
}

func usage() {
	fmt.Print(
		"multiFaAcceleration - Performs velocity and acceleration on a four way multiple alignment in multiFa format." +
			"A four way multiple alignment must contain four species (index 0 to 3) in the topology that aln[0] is the most derived and species 1 to 3 are successive outgroups.\n" +
			"Three bed files are returned. The first produces the velocity score, the second returns the acceleration score, and the third returns the initial velocity score for each window of the genome for aln[0].\n" +
			"Usage:\n" +
			" multiFaAcceleration chromName in.fa velocity.bed acceleration.bed initialVelocity.bed\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 5
	var searchSpaceBed *string = flag.String("searchSpaceBed", "", "Limits the generation of data to windows contained within these regions.")
	var searchSpaceProportion *float64 = flag.Float64("searchSpaceProportion", 0.5, "Proportion of window that must overlap search space in order to be evaluated.")
	var windowSize *int = flag.Int("windowSize", 500, "Set the size of the sliding window.")
	var useSnpDistance *bool = flag.Bool("useSnpDistance", false, "Calculate pairwise distances with SNPs instead of the default mutation distance, which counts INDELs.")
	var verbose *bool = flag.Bool("verbose", false, "Enables debug prints.")
	var epsilon *float64 = flag.Float64("epsilon", 1e-8, "Set the error threshold for alternating least squares branch length calculation.")
	var allowNegative *bool = flag.Bool("allowNegative", false, "Allow the algorithm to evaluate negative branch lengths. This program will constrain the optimal solution to non-negative branch lengths by default.")
	var zeroDistanceWeightConstant *float64 = flag.Float64("zeroDistanceWeightConstant", 1000, "Set the relative error weight applied to pairs of species with a pairwise distance of zero.")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	chromName := flag.Arg(0)
	inFile := flag.Arg(1)
	velOut := flag.Arg(2)
	accelOut := flag.Arg(3)
	initialVOut := flag.Arg(4)

	s := Settings{
		InFile:                     inFile,
		ChromName:                  chromName,
		VelOut:                     velOut,
		AccelOut:                   accelOut,
		InitialVelOut:              initialVOut,
		SearchSpaceBed:             *searchSpaceBed,
		SearchSpaceProportion:      *searchSpaceProportion,
		UseSnpDistance:             *useSnpDistance,
		WindowSize:                 *windowSize,
		Verbose:                    *verbose,
		Epsilon:                    *epsilon,
		AllowNegative:              *allowNegative,
		ZeroDistanceWeightConstant: *zeroDistanceWeightConstant,
	}

	multiFaAcceleration(s)
}

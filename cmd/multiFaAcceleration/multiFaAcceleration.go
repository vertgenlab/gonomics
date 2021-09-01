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

//Settings contains all program arguments and options, and is passed throughout helper functions to aid readability.
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

//Once we have branch lengths for each valid window, we will need to normalize the values relative to each other. Thus, we store the branch lengths in this intermediate cache before writing to file.
type BranchCache struct {
	ChromStart int
	ChromEnd   int
	B1         float64
	B3         float64
}

//A set of all observed pairwise distances between the four species. While these numbers must be integers, we store them as float64 to avoid casting in other functions as a way to improve readability.
type Distances struct {
	D01 float64
	D02 float64
	D03 float64
	D12 float64
	D13 float64
	D23 float64
}

//The set of branch lengths corresponding to a particular distance matrix.
type BranchLengths struct {
	B1 float64
	B2 float64
	B3 float64
	B4 float64
	B5 float64
}

//SubTree is a tree with three leaves and three branches, joined at the single internal node.
//Dij represent the observed pairwise distances between two species. vi represents the length of the branch between a species i and the internal node at the current stage of optimization.
type SubTree struct {
	Dab float64
	Dac float64
	Dbc float64
	Va  float64
	Vb  float64
	Vc  float64
}

func multiFaAcceleration(s Settings) {
	records := fasta.Read(s.InFile)
	var searchSpace []bed.Bed
	var referenceLength, i, j, threshold int
	var bitArray []bool

	//a couple of safety checks
	if len(records) != 4 {
		log.Fatalf("multiFaAcceleration accepts a multiFa file with 4 records, found %v.", len(records))
	}
	if len(records[1].Seq) != len(records[0].Seq) || len(records[2].Seq) != len(records[0].Seq) || len(records[3].Seq) != len(records[0].Seq) {
		log.Fatalf("Error. All records must be of the same sequence length.")
	}

	referenceLength = fasta.AlnPosToRefPos(records[0], len(records[0].Seq)-1)

	if s.SearchSpaceBed != "" {
		//if we want to limit our search to a particular set of regions, we populate the bitArray with 1s instead of 0s at positions overlapping the input bed entries.
		searchSpace = bed.Read(s.SearchSpaceBed)
		bitArray = make([]bool, referenceLength)
		for i = range searchSpace {
			if searchSpace[i].Chrom == s.ChromName {
				for j = searchSpace[i].ChromStart; j < searchSpace[i].ChromEnd; j++ {
					bitArray[j] = true
				}
			}
		}
		threshold = int(s.SearchSpaceProportion * float64(s.WindowSize)) //the minimum number of bases at which a window must overlap the search space in order to be considered a valid window.
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
		currCount, pass = thresholdCheckPasses(s, currCount, threshold, bitArray, referenceCounter) //first we see if we have a valid window.
		if records[0].Seq[alignmentCounter] != dna.Gap {                                            //and if we are at a reference position.
			if pass {
				if s.UseSnpDistance {
					reachedEnd = fourWaySnpDistances(records, alignmentCounter, s, &currDistances)
				} else {
					reachedEnd = fourWayMutationDistances(records, alignmentCounter, s, &currDistances)
				}

				if _, containedInMap = distanceCache[currDistances]; !containedInMap { //if this tree has not been seen before, calculate branch lengths
					distanceCache[currDistances] = alternatingLeastSquares(currDistances, s)
				}
				//now our distances should be in the cache, and we can do a simple lookup.
				b1 = distanceCache[currDistances].B1
				b3 = distanceCache[currDistances].B3

				if !reachedEnd { //if we haven't run out of chromosome, we have another valid window to appear in our output.
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

	//with our normalization parameters calculated, we can normalize the branch lengths for each window and write to file.
	for i = range branchCacheSlice {
		b1Normal = branchCacheSlice[i].B1 / averageVel
		b3Normal = branchCacheSlice[i].B3 / averageInitialVel
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

//this helper function calculates the optimal branch lengths for a given set of distances. See readme for a detailed description of this algorithm.
func alternatingLeastSquares(d Distances, s Settings) BranchLengths {
	var answer = BranchLengths{1, 1, 1, 1, 1}
	var Q float64 = calculateQ(d, answer, s)
	var nextQ float64
	var currDiff float64 = s.Epsilon + 1 //set currDiff to something larger than epsilon so that we make it into the loop the first time.
	var sub SubTree
	var maxIteration, i = 1000, 0
	var oldAnswer = BranchLengths{1, 1, 1, 1, 1}

	for currDiff > s.Epsilon && i < maxIteration {
		oldAnswer = answer
		pruneLeft(d, answer, &sub, s.ZeroDistanceWeightConstant)
		answer.B1, answer.B2, answer.B3 = optimizeSubtree(&sub, s)
		pruneRight(d, answer, &sub, s.ZeroDistanceWeightConstant)
		answer.B4, answer.B5, answer.B3 = optimizeSubtree(&sub, s)
		nextQ = calculateQ(d, answer, s)
		//DEBUG: log.Printf("nextQ: %e. currDiff: %e. Here were the branch lengths: %f. %f. %f. %f. %f.", nextQ, currDiff, answer.B1, answer.B2, answer.B3, answer.B4, answer.B5)
		currDiff = math.Abs(Q - nextQ)
		if nextQ > Q { //nextQ is higher than Q, which means we got "worse"
			answer = oldAnswer //we will exit the loop next time, so we want the old answer, which has the lower of the two terminal Q estimates.
			currDiff = 0
		}
		Q = nextQ
		i++
	}
	if i >= maxIteration {
		log.Fatalf("Failed to converge on a tree with these distances. D01: %f, D02: %f, D03: %f, D12: %f, D13: %f, D23: %f.", d.D01, d.D02, d.D03, d.D12, d.D13, d.D23)
	}
	return answer
}

//a helper function of alternatingLeastSquares. Calculates the optimal branch lengths for the three branches in a subtree.
func optimizeSubtree(sub *SubTree, s Settings) (float64, float64, float64) {
	sub.Va = (sub.Dab + sub.Dac - sub.Dbc) / 2.0
	sub.Vb = (sub.Dab + sub.Dbc - sub.Dac) / 2.0
	sub.Vc = (sub.Dac + sub.Dbc - sub.Dac) / 2.0

	if s.AllowNegative {
		return sub.Va, sub.Vb, sub.Vc
	}
	if sub.Va < 0 && sub.Vb < 0 && sub.Vc < 0 {
		if s.Verbose {
			log.Printf("WARNING: All branches are negative.") //TODO: Should this error out?
		}
		sub.Va, sub.Vb, sub.Vc = 0, 0, 0
	} else if sub.Va < 0 && sub.Vb < 0 {
		sub.Va = 0
		sub.Vb = 0
		sub.Vc = nonNegativeApproximation(sub.Dac, sub.Dbc, sub.Va, sub.Vb, s.ZeroDistanceWeightConstant)
	} else if sub.Va < 0 && sub.Vc < 0 {
		sub.Va = 0
		sub.Vc = 0
		sub.Vb = nonNegativeApproximation(sub.Dbc, sub.Dab, sub.Vc, sub.Va, s.ZeroDistanceWeightConstant)
	} else if sub.Vb < 0 && sub.Vc < 0 {
		sub.Vb = 0
		sub.Vc = 0
		sub.Va = nonNegativeApproximation(sub.Dab, sub.Dac, sub.Vb, sub.Vc, s.ZeroDistanceWeightConstant)
	} else if sub.Va < 0 {
		sub.Va = 0
		sub.Vb = nonNegativeApproximation(sub.Dab, sub.Dbc, sub.Va, sub.Vc, s.ZeroDistanceWeightConstant)
		sub.Vc = nonNegativeApproximation(sub.Dac, sub.Dbc, sub.Va, sub.Vb, s.ZeroDistanceWeightConstant)
	} else if sub.Vb < 0 {
		sub.Vb = 0
		sub.Va = nonNegativeApproximation(sub.Dab, sub.Dac, sub.Vb, sub.Vc, s.ZeroDistanceWeightConstant)
		sub.Vc = nonNegativeApproximation(sub.Dab, sub.Dbc, sub.Va, sub.Vc, s.ZeroDistanceWeightConstant)
	} else if sub.Vc < 0 {
		sub.Vc = 0
		sub.Va = nonNegativeApproximation(sub.Dab, sub.Dac, sub.Vb, sub.Vc, s.ZeroDistanceWeightConstant)
		sub.Vb = nonNegativeApproximation(sub.Dab, sub.Dbc, sub.Va, sub.Vc, s.ZeroDistanceWeightConstant)
	}
	return sub.Va, sub.Vb, sub.Vc
}

//If we constrain branch lengths to be nonNegative, we apply this correction when the minimum Q is achieved at negative branch lengths for a subtree.
func nonNegativeApproximation(d1 float64, d2 float64, v1 float64, v2 float64, ZeroDistanceWeightConstant float64) float64 {
	if d1 == 0 {
		if d2 == 0 {
			return numbers.MaxFloat64(0, (ZeroDistanceWeightConstant*(d1-v1)+ZeroDistanceWeightConstant*(d2-v2))/(2*ZeroDistanceWeightConstant))
		} else {
			return numbers.MaxFloat64(0, ZeroDistanceWeightConstant*(d1-v1)+(1.0/math.Pow(d2, 2)*(d2-v2))) / (ZeroDistanceWeightConstant + (1.0 / math.Pow(d2, 2)))
		}
	} else if d2 == 0 {
		return numbers.MaxFloat64(0, (1.0/(math.Pow(d1, 2))*(d1-v1)+ZeroDistanceWeightConstant*(d2-v2))/((1.0/math.Pow(d1, 2))+ZeroDistanceWeightConstant))
	}
	return numbers.MaxFloat64(0, (1.0/(math.Pow(d1, 2))*(d1-v1)+(1.0/math.Pow(d2, 2)*(d2-v2)))/((1.0/math.Pow(d1, 2))+(1.0/math.Pow(d2, 2))))
}

//Reduce the four species tree to the subtree containing species 0, 1, and the ancestor of 2/3.
func pruneLeft(d Distances, b BranchLengths, sub *SubTree, ZeroDistanceWeightConstant float64) {
	sub.Dab = d.D01
	if d.D03 == 0 {
		if d.D02 == 0 {
			sub.Dac = (ZeroDistanceWeightConstant*(d.D02-b.B4) + ZeroDistanceWeightConstant*(d.D03-b.B5)) / (2 * ZeroDistanceWeightConstant)
		} else {
			sub.Dac = ((1.0/math.Pow(d.D02, 2))*(d.D02-b.B4) + (ZeroDistanceWeightConstant * (d.D03 - b.B5))) / (ZeroDistanceWeightConstant + (1.0 / math.Pow(d.D02, 2)))
		}
	} else if d.D02 == 0 {
		sub.Dac = (ZeroDistanceWeightConstant*(d.D02-b.B4) + ((1.0 / math.Pow(d.D03, 2)) * (d.D03 - b.B5))) / ((1.0 / math.Pow(d.D03, 2)) + ZeroDistanceWeightConstant)
	} else {
		sub.Dac = ((1.0/math.Pow(d.D02, 2))*(d.D02-b.B4) + (1.0/math.Pow(d.D03, 2))*(d.D03-b.B5)) / ((1.0 / math.Pow(d.D03, 2)) + (1.0 / math.Pow(d.D02, 2)))
	}

	if d.D13 == 0 {
		if d.D12 == 0 {
			sub.Dbc = (ZeroDistanceWeightConstant*(d.D12-b.B4) + ZeroDistanceWeightConstant*(d.D13-b.B5)) / (2 * ZeroDistanceWeightConstant)
		} else {
			sub.Dbc = (ZeroDistanceWeightConstant*(d.D13-b.B5) + (1.0/math.Pow(d.D12, 2))*(d.D12-b.B4)) / ((1.0 / math.Pow(d.D12, 2)) + ZeroDistanceWeightConstant)
		}
	} else if d.D12 == 0 {
		sub.Dbc = (ZeroDistanceWeightConstant*(d.D12-b.B4) + (1.0/math.Pow(d.D13, 2))*(d.D13-b.B5)) / ((1.0 / math.Pow(d.D13, 2)) + ZeroDistanceWeightConstant)
	} else {
		sub.Dbc = ((1.0/math.Pow(d.D12, 2))*(d.D12-b.B4) + (1.0/math.Pow(d.D13, 2))*(d.D13-b.B5)) / ((1.0 / math.Pow(d.D13, 2)) + (1.0 / math.Pow(d.D12, 2)))
	}
}

//Reduce the four species tree to the subtree containing 2, 3, and the ancestor of 0/1.
func pruneRight(d Distances, b BranchLengths, sub *SubTree, ZeroDistanceWeightConstant float64) {
	sub.Dac = d.D23

	if d.D02 == 0 {
		if d.D12 == 0 {
			sub.Dac = (ZeroDistanceWeightConstant*(d.D02-b.B1) + ZeroDistanceWeightConstant*(d.D12-b.B2)) / (2.0 * ZeroDistanceWeightConstant)
		} else {
			sub.Dac = ((1.0/math.Pow(d.D12, 2))*(d.D12-b.B2) + ZeroDistanceWeightConstant*(d.D02-b.B1)) / ((1.0 / math.Pow(d.D12, 2)) + ZeroDistanceWeightConstant)
		}
	} else if d.D12 == 0 {
		sub.Dac = ((1.0/math.Pow(d.D02, 2))*(d.D02-b.B1) + ZeroDistanceWeightConstant*(d.D12-b.B2)) / ((1.0 / math.Pow(d.D02, 2)) + ZeroDistanceWeightConstant)
	} else {
		sub.Dab = ((1.0/math.Pow(d.D02, 2))*(d.D02-b.B1) + (1.0/math.Pow(d.D12, 2))*(d.D12-b.B2)) / ((1.0 / math.Pow(d.D02, 2)) + (1.0 / math.Pow(d.D12, 2)))
	}

	if d.D03 == 0 {
		if d.D13 == 0 {
			sub.Dbc = (ZeroDistanceWeightConstant*(d.D03-b.B1) + ZeroDistanceWeightConstant*(d.D13-b.B2)) / (2.0 * ZeroDistanceWeightConstant)
		} else {
			sub.Dbc = ((1.0/math.Pow(d.D13, 2))*(d.D13-b.B2) + ZeroDistanceWeightConstant*(d.D03-b.B1)) / ((1.0 / math.Pow(d.D13, 2)) + ZeroDistanceWeightConstant)
		}
	} else if d.D13 == 0 {
		sub.Dbc = ((1.0/math.Pow(d.D03, 2))*(d.D03-b.B1) + ZeroDistanceWeightConstant*(d.D13-b.B2)) / ((1.0 / math.Pow(d.D03, 2)) + ZeroDistanceWeightConstant)
	} else {
		sub.Dbc = ((1.0/math.Pow(d.D03, 2))*(d.D03-b.B1) + (1.0/math.Pow(d.D13, 2))*(d.D13-b.B2)) / ((1.0 / math.Pow(d.D03, 2)) + (1.0 / math.Pow(d.D13, 2)))
	}
}

//For a set of distances and corresponding branch lengths, determine the value of Q, the Fitch-Margoliash least squares error.
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
func thresholdCheckPasses(s Settings, currCount int, threshold int, bitArray []bool, referenceCounter int) (int, bool) {
	if s.SearchSpaceBed == "" { //no search space file, no need to look further
		return 0, true
	}
	if referenceCounter == 0 {
		currCount = 0
		for i := 0; i < s.WindowSize; i++ {
			if bitArray[i] {
				currCount++
			}
		}
	} else {
		if bitArray[referenceCounter-1] {
			currCount--
		}
		if bitArray[referenceCounter+s.WindowSize-1] {
			currCount++
		}
	}
	return currCount, currCount >= threshold
}

//Generate distances from mutation distances, which includes SNPs and INDELs, where each INDEL counts as one mutation regardless of length.
func fourWayMutationDistances(records []fasta.Fasta, alignmentCounter int, s Settings, D *Distances) bool {
	//first we clear the values in D.
	D.D01, D.D02, D.D03, D.D12, D.D13, D.D23 = 0, 0, 0, 0, 0, 0
	var D01tmp int
	var reachedEnd bool
	var alnEnd int
	D01tmp, reachedEnd, alnEnd = fasta.PairwiseMutationDistanceReferenceWindow(records[0], records[1], alignmentCounter, s.WindowSize)
	D.D01 = float64(D01tmp)
	D.D02 = float64(fasta.PairwiseMutationDistanceInRange(records[0], records[2], alignmentCounter, alnEnd))
	D.D03 = float64(fasta.PairwiseMutationDistanceInRange(records[0], records[3], alignmentCounter, alnEnd))
	D.D12 = float64(fasta.PairwiseMutationDistanceInRange(records[1], records[2], alignmentCounter, alnEnd))
	D.D13 = float64(fasta.PairwiseMutationDistanceInRange(records[1], records[3], alignmentCounter, alnEnd))
	D.D23 = float64(fasta.PairwiseMutationDistanceInRange(records[2], records[3], alignmentCounter, alnEnd))
	return reachedEnd
}

//Generate distances from SNP distance, which includes only SNPs.
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

//a helper function of fourWaySnpDistances, determines if an alignment column is comprised of bases (not gaps) for each species.
func isUngappedColumn(records []fasta.Fasta, index int) bool {
	for i := range records {
		if !isUngappedBase(records[i].Seq[index]) {
			return false
		}
	}
	return true
}

//a helper function of isUngappedColumn. True if a dna.Base is a base, not an N, gap, or dot.
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

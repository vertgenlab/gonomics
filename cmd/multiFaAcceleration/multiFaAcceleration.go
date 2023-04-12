// Command Group: "Statistics & Population Genetics"

package main

import (
	"flag"
	"fmt"
	"log"

	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/phylo"
)

// Settings contains all program arguments and options, and is passed throughout helper functions to aid readability.
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
	RawVelOut                  string
	RawInitialOut              string
	CavalliSforzaQ             bool
}

// Once we have branch lengths for each valid window, we will need to normalize the values relative to each other. Thus, we store the branch lengths in this intermediate cache before writing to file.
type BranchCache struct {
	ChromStart int
	ChromEnd   int
	BhumHca    float64
	BhcaHga    float64
}

func multiFaAcceleration(s Settings) {
	records := fasta.Read(s.InFile)
	var referenceLength, i, threshold int
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
		bitArray = phylo.MakeBitArrayFromSearchSpaceBed(s.SearchSpaceBed, referenceLength, s.ChromName)
		threshold = int(s.SearchSpaceProportion * float64(s.WindowSize)) //the minimum number of bases at which a window must overlap the search space in order to be considered a valid window.
	}

	var currDistances phylo.AccelDistancesAndWeights
	var distanceCache = make(map[phylo.AccelDistancesAndWeights]phylo.AccelBranchLengths)
	var referenceCounter int = 0
	var reachedEnd bool = false
	var BhumHca, BhcaHga float64
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
					reachedEnd = phylo.AccelFourWaySnpDistancesAndWeights(records, alignmentCounter, s.WindowSize, &currDistances, s.ZeroDistanceWeightConstant, s.CavalliSforzaQ)
				} else {
					reachedEnd = phylo.AccelFourWayMutationDistancesAndWeights(records, alignmentCounter, s.WindowSize, &currDistances, s.ZeroDistanceWeightConstant, s.CavalliSforzaQ)
				}

				if _, containedInMap = distanceCache[currDistances]; !containedInMap { //if this tree has not been seen before, calculate branch lengths
					distanceCache[currDistances] = phylo.BranchLengthsAlternatingLeastSquares(currDistances, s.AllowNegative, s.Verbose, s.ZeroDistanceWeightConstant, s.Epsilon, s.CavalliSforzaQ)
				}
				if s.Verbose {
					fmt.Printf("RefPos: %v\tDistances: %v\t Branches: %v.\n", referenceCounter, currDistances, distanceCache[currDistances])
				}
				//now our distances should be in the cache, and we can do a simple lookup.
				BhumHca = distanceCache[currDistances].BhumHca
				BhcaHga = distanceCache[currDistances].BhcaHga

				if !reachedEnd { //if we haven't run out of chromosome, we have another valid window to appear in our output.
					velSum += BhumHca
					initialSum += BhcaHga
					branchCacheSlice = append(branchCacheSlice, BranchCache{referenceCounter, referenceCounter + s.WindowSize, BhumHca, BhcaHga})
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
	var RawVelOutBed, RawInitialOutBed *fileio.EasyWriter
	if s.RawVelOut != "" {
		RawVelOutBed = fileio.EasyCreate(s.RawVelOut)
	}
	if s.RawInitialOut != "" {
		RawInitialOutBed = fileio.EasyCreate(s.RawInitialOut)
	}

	//with our normalization parameters calculated, we can normalize the branch lengths for each window and write to file.
	for i = range branchCacheSlice {
		b1Normal = branchCacheSlice[i].BhumHca / averageVel
		b3Normal = branchCacheSlice[i].BhcaHga / averageInitialVel
		bed.WriteBed(velBed, bed.Bed{Chrom: s.ChromName, ChromStart: branchCacheSlice[i].ChromStart, ChromEnd: branchCacheSlice[i].ChromEnd, Name: fmt.Sprintf("%.8g", b1Normal), FieldsInitialized: 4})
		bed.WriteBed(initialVelBed, bed.Bed{Chrom: s.ChromName, ChromStart: branchCacheSlice[i].ChromStart, ChromEnd: branchCacheSlice[i].ChromEnd, Name: fmt.Sprintf("%.8g", b3Normal), FieldsInitialized: 4})
		bed.WriteBed(accelBed, bed.Bed{Chrom: s.ChromName, ChromStart: branchCacheSlice[i].ChromStart, ChromEnd: branchCacheSlice[i].ChromEnd, Name: fmt.Sprintf("%.8g", b1Normal-b3Normal), FieldsInitialized: 4})
		if s.RawVelOut != "" {
			bed.WriteBed(RawVelOutBed, bed.Bed{Chrom: s.ChromName, ChromStart: branchCacheSlice[i].ChromStart, ChromEnd: branchCacheSlice[i].ChromEnd, Name: fmt.Sprintf("%.8g", branchCacheSlice[i].BhumHca), FieldsInitialized: 4})
		}
		if s.RawInitialOut != "" {
			bed.WriteBed(RawInitialOutBed, bed.Bed{Chrom: s.ChromName, ChromStart: branchCacheSlice[i].ChromStart, ChromEnd: branchCacheSlice[i].ChromEnd, Name: fmt.Sprintf("%.8g", branchCacheSlice[i].BhcaHga), FieldsInitialized: 4})
		}
	}

	err = velBed.Close()
	exception.PanicOnErr(err)
	err = accelBed.Close()
	exception.PanicOnErr(err)
	err = initialVelBed.Close()
	exception.PanicOnErr(err)
	if s.RawVelOut != "" {
		err = RawVelOutBed.Close()
		exception.PanicOnErr(err)
	}
	if s.RawInitialOut != "" {
		err = RawInitialOutBed.Close()
		exception.PanicOnErr(err)
	}
}

// this helper function returns true if a sufficient number of positions in the candidate window overlap the search space, as measured by non-zero values in the bitArray.
// bitArray is on reference coordinates, not alignment coordinates, so the window is simply equal to windowSize.
// the first return corresponds to the number of positions that overlap the search space.
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
	var rawVelBranchLength *string = flag.String("rawVelBranchLengths", "", "Set an output file name to return the raw branch length for the branch associated with the velocity score.")
	var rawInitialVelBranchLength *string = flag.String("rawInitialVelBranchLengths", "", "Set an output file name to return the raw branch length for the branch associated with the initial velocity score.")
	var cavalliSforzaQ *bool = flag.Bool("CavalliSforzaEdwardsQ", false, "Usees the CavalliSforza-Edwards error term Q instead of the Fitch-Margoliash error term.")

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
		RawVelOut:                  *rawVelBranchLength,
		RawInitialOut:              *rawInitialVelBranchLength,
		CavalliSforzaQ:             *cavalliSforzaQ,
	}

	multiFaAcceleration(s)
}

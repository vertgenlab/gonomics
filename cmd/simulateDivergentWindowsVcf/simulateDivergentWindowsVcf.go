package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/simulate"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
	"math/rand"
	"os"
	"sort"
	"strings"
)

func init() {
	dbgVar := os.Getenv("GODEBUG")

	if dbgVar == "" {
		err := os.Setenv("GODEBUG", "randautoseed=0")
		exception.PanicOnErr(err)
		return
	}

	variables := strings.Split(dbgVar, ",")
	var found bool
	for i := range variables {
		if strings.HasPrefix(variables[i], "randautoseed=") {
			variables[i] = "randautoseed=0"
			found = true
		}
	}

	if !found {
		variables = append(variables, "randautoseed=0")
	}

	err := os.Setenv("GODEBUG", strings.Join(variables, ","))
	exception.PanicOnErr(err)
	rand.Seed(0)
}

type Window struct {
	NumDivergent int
	Variants     []vcf.Vcf
}

func simulateDivergentWindowsVcf(s Settings) {
	if s.NumWindowSites > s.NumTotalSites {
		log.Fatalf("The number of total simulated Vcf variants must be greater than the desired number of variants per window.")
	}
	if s.UpperPercentile > 1 || s.UpperPercentile < 0 {
		log.Fatalf("UpperPercentile argument must be between one and zero.")
	}
	if s.LowerPercentile > 1 || s.LowerPercentile < 0 {
		log.Fatalf("LowerPercentile argument must be between one and zero.")
	}

	rand.Seed(s.SetSeed)
	var err error
	var TotalSites []vcf.Vcf = make([]vcf.Vcf, s.NumTotalSites)
	var windows []Window = make([]Window, s.NumWindows)

	for i := 0; i < s.NumTotalSites; i++ {
		TotalSites[i] = simulate.SingleVcf(s.Alpha, s.NumAlleles, s.BoundAlpha, s.BoundBeta, s.BoundMultiplier, i+1)
	}

	for i := 0; i < s.NumWindows; i++ {
		windows[i].Variants = make([]vcf.Vcf, s.NumWindowSites)
		//Shuffle the vcf records, our subset will be composed to the first entries in the shuffled order.
		rand.Seed(s.SetSeed * int64(i))
		rand.Shuffle(len(TotalSites), func(i, j int) { TotalSites[i], TotalSites[j] = TotalSites[j], TotalSites[i] })
		copy(windows[i].Variants, TotalSites[:s.NumWindowSites]) //keep only as many results as specified
		windows[i].NumDivergent = countDivergent(windows[i].Variants)
	}

	//now sort the windows by the number of divergent sites from low to high
	sort.Slice(windows, func(i, j int) bool { return windows[i].NumDivergent < windows[j].NumDivergent })
	lowerOut := fileio.EasyCreate(s.LowerOut)

	for i := 0; i < int(s.LowerPercentile*float64(s.NumWindows)); i++ {
		for j := range windows[i].Variants {
			vcf.WriteVcf(lowerOut, windows[i].Variants[j])
		}
	}
	err = lowerOut.Close()
	exception.PanicOnErr(err)

	upperOut := fileio.EasyCreate(s.UpperOut)
	for i := int(s.UpperPercentile * float64(s.NumWindows)); i < len(windows); i++ {
		for j := range windows[i].Variants {
			vcf.WriteVcf(upperOut, windows[i].Variants[j])
		}
	}
	err = upperOut.Close()
	exception.PanicOnErr(err)
}

// countDivergent returns the number of variants in a slice of Vcf structs that are in the divergent state.
func countDivergent(v []vcf.Vcf) int {
	var answer int = 0
	for i := range v {
		if vcf.IsAltAncestor(v[i]) {
			answer++
		}
	}
	return answer
}

func usage() {
	fmt.Print(
		"simulateDivergentWindowsVcf - Simulates Vcf data, partitioned by divergence-based ascertainment.\n" +
			"Usage:\n" +
			" simulateDivergentWindowsVcf [options] upperDiverge.vcf lowerDiverge.vcf\n" +
			"options:\n")
	flag.PrintDefaults()
}

type Settings struct {
	UpperOut        string
	LowerOut        string
	Alpha           float64
	NumAlleles      int
	NumTotalSites   int
	NumWindowSites  int
	NumWindows      int
	SetSeed         int64
	BoundAlpha      float64
	BoundBeta       float64
	BoundMultiplier float64
	UpperPercentile float64
	LowerPercentile float64
}

func main() {
	var expectedNumArgs int = 2
	var numTotalSites *int = flag.Int("numTotalSites", 10000, "Specifies the total number of simulated variants used to pick windows. Must be larger than NumWindowSites.")
	var numWindowSites *int = flag.Int("numWindowSites", 100, "Specifies the number of segregating sites per window.")
	var numWindows *int = flag.Int("numWindows", 1000, "Specifies the number of windows.")
	var setSeed *int64 = flag.Int64("setSeed", -1, "Use a specific seed for the RNG.")
	var alpha *float64 = flag.Float64("alpha", 0.01, "Specifies the selection parameter alpha for drawing individual gVCF alleles from a stationarity distribution.")
	var numAlleles *int = flag.Int("numAlleles", 100, "Specifies the number of alleles for gVCF samples.")
	var boundAlpha *float64 = flag.Float64("boundAlpha", 0.001, "Set the alpha parameter for the bounding function.")
	var boundBeta *float64 = flag.Float64("boundBeta", 0.001, "Set the beta parameter for the bounding function.")
	var boundMultiplier *float64 = flag.Float64("boundMultiplier", 10000, "Set the multiplier for the bounding function.")
	var upperPercentile *float64 = flag.Float64("upperPercentile", 0.99, "Set the percentile for variants in the upper divergence set.")
	var lowerPercentile *float64 = flag.Float64("lowerPercentile", 0.01, "Set the percentile for variants in the lower divergence set.")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()
	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}
	upperOut := flag.Arg(0)
	lowerOut := flag.Arg(1)

	s := Settings{
		UpperOut:        upperOut,
		LowerOut:        lowerOut,
		Alpha:           *alpha,
		NumAlleles:      *numAlleles,
		NumTotalSites:   *numTotalSites,
		NumWindowSites:  *numWindowSites,
		NumWindows:      *numWindows,
		SetSeed:         *setSeed,
		BoundAlpha:      *boundAlpha,
		BoundBeta:       *boundBeta,
		BoundMultiplier: *boundMultiplier,
		UpperPercentile: *upperPercentile,
		LowerPercentile: *lowerPercentile,
	}

	simulateDivergentWindowsVcf(s)
}

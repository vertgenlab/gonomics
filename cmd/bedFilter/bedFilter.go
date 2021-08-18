// Command Group: "BED Tools"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"log"
	"math"
	"math/rand"
)

func bedFilter(s Settings) {
	common.RngSeed(s.RandSeed, s.SetSeed)
	var length int
	var pass bool = false
	var r float64
	bedChan := bed.GoReadToChan(s.InFile)
	out := fileio.EasyCreate(s.OutFile)

	for curr := range bedChan {
		pass = true
		length = curr.ChromEnd - curr.ChromStart
		if curr.FieldsInitialized > 4 {
			if curr.Score < s.MinScore {
				pass = false
			}
			if curr.Score > s.MaxScore {
				pass = false
			}
		}
		if length < s.MinLength {
			pass = false
		}
		if length > s.MaxLength {
			pass = false
		}
		if curr.ChromStart < s.MinStart {
			pass = false
		}
		if curr.ChromStart > s.MaxStart {
			pass = false
		}
		if curr.ChromEnd < s.MinEnd {
			pass = false
		}
		if curr.ChromEnd > s.MaxEnd {
			pass = false
		}
		if s.MinNameFloat != -1*math.MaxFloat64 {
			if common.StringToFloat64(curr.Name) < s.MinNameFloat {
				pass = false
			}
		}
		if s.MaxNameFloat != math.MaxFloat64 {
			if common.StringToFloat64(curr.Name) > s.MaxNameFloat {
				pass = false
			}
		}
		if s.Chrom != "" {
			if curr.Chrom != s.Chrom {
				pass = false
			}
		}
		if s.SubSet < 1.0 {
			r = rand.Float64()
			if r > s.SubSet {
				pass = false
			}
		}
		if pass {
			bed.WriteBed(out, curr)
		}
	}
	var err error
	err = out.Close()
	exception.PanicOnErr(err)
}

func usage() {
	fmt.Print(
		"bedFilter\n" +
			"Usage:\n" +
			"bedFilter input.bed output.bed\n" +
			"options:\n")
	flag.PrintDefaults()
}

type Settings struct {
	InFile    string
	OutFile   string
	MinScore  int
	MaxScore  int
	MinLength int
	MaxLength int
	MinStart  int
	MaxStart  int
	MinEnd    int
	MaxEnd    int
	MinNameFloat float64
	MaxNameFloat float64
	Chrom     string
	SubSet    float64
	RandSeed  bool
	SetSeed   int64
}

func main() {
	var expectedNumArgs int = 2
	var minScore *int = flag.Int("minScore", -1*numbers.MaxInt, "Specifies the minimum score in the fourth field.")
	var maxScore *int = flag.Int("maxScore", numbers.MaxInt, "Specifies the maximum score in the fourth field.")
	var minLength *int = flag.Int("minLength", 0, "Specifies the minimum length of the region.")
	var maxLength *int = flag.Int("maxLength", numbers.MaxInt, "Specifies the maximum length of the region.")
	var minStart *int = flag.Int("minStart", 0, "Specifies the minimum starting position of the region.")
	var maxStart *int = flag.Int("maxStart", numbers.MaxInt, "Specifies the maximum starting position of the region.")
	var minEnd *int = flag.Int("minEnd", 0, "Specifies the minimum ending position of the region.")
	var maxEnd *int = flag.Int("maxEnd", numbers.MaxInt, "Specifies the maximum ending position of the region.")
	var minNameFloat *float64 = flag.Float64("minNameFloat", -1*math.MaxFloat64, "Specifies the minimum floating point number value for bed entries where floating point numbers are stored in the name field.")
	var maxNameFloat *float64 = flag.Float64("maxNameFloat", math.MaxFloat64, "Specifies the maximum floating point number value for bed entries where floating point numbers are stored in the name field.")
	var chrom *string = flag.String("chrom", "", "Specifies the chromosome name.")
	var subSet *float64 = flag.Float64("subSet", 1.0, "Proportion of entries to retain in output, range from 0 to 1.")
	var randSeed *bool = flag.Bool("randSeed", false, "Uses a random seed for the RNG.")
	var setSeed *int64 = flag.Int64("setSeed", -1, "Use a specific seed for the RNG.")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()
	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}
	infile := flag.Arg(0)
	outfile := flag.Arg(1)
	s := Settings{
		InFile:    infile,
		OutFile:   outfile,
		MinScore:  *minScore,
		MaxScore:  *maxScore,
		MinLength: *minLength,
		MaxLength: *maxLength,
		MinStart:  *minStart,
		MaxStart:  *maxStart,
		MinEnd:    *minEnd,
		MaxEnd:    *maxEnd,
		MinNameFloat: *minNameFloat,
		MaxNameFloat: *maxNameFloat,
		Chrom:     *chrom,
		SubSet:    *subSet,
		RandSeed:  *randSeed,
		SetSeed:   *setSeed,
	}
	bedFilter(s)
}

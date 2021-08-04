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
)

type Settings struct {
	InFile string
	OutFile string
	MinScore int
	MaxScore int
	MinLength int
	MaxLength int
	MinStart int
	MaxStart int
	MinEnd int
	MaxEnd int
	MinNameFloat float64
	MaxNameFloat float64
	Chrom string
}

func bedFilter(s Settings) {
	var pass bool = false
	var length int
	out := fileio.EasyCreate(s.OutFile)

	bedChan := bed.GoReadToChan(s.InFile)
	for i := range bedChan {
		pass = true
		length = i.ChromEnd - i.ChromStart
		if i.Score < s.MinScore {
			pass = false
		}
		if i.Score > s.MaxScore {
			pass = false
		}
		if length < s.MinLength {
			pass = false
		}
		if length > s.MaxLength {
			pass = false
		}
		if i.ChromStart < s.MinStart {
			pass = false
		}
		if i.ChromStart > s.MaxStart {
			pass = false
		}
		if i.ChromEnd < s.MinEnd {
			pass = false
		}
		if i.ChromEnd > s.MaxEnd {
			pass = false
		}
		if s.MinNameFloat != -1 * math.MaxFloat64 {
			if common.StringToFloat64(i.Name) < s.MinNameFloat {
				pass = false
			}
		}
		if s.MaxNameFloat != math.MaxFloat64 {
			if common.StringToFloat64(i.Name) > s.MaxNameFloat {
				pass = false
			}
		}
		if s.Chrom != "" {
			if s.Chrom != i.Chrom {
				pass = false
			}
		}
		if pass {
			bed.WriteBed(out, i)
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

func main() {
	var expectedNumArgs int = 2
	var minScore *int = flag.Int("minScore", -1 * numbers.MaxInt, "Specifies the minimum score in the fourth field.")
	var maxScore *int = flag.Int("maxScore", numbers.MaxInt, "Specifies the maximum score in the fourth field.")
	var minLength *int = flag.Int("minLength", 0, "Specifies the minimum length of the region.")
	var maxLength *int = flag.Int("maxLength", numbers.MaxInt, "Specifies the maximum length of the region.")
	var minStart *int = flag.Int("minStart", 0, "Specifies the minimum starting position of the region.")
	var maxStart *int = flag.Int("maxStart", numbers.MaxInt, "Specifies the maximum starting position of the region.")
	var minEnd *int = flag.Int("minEnd", 0, "Specifies the minimum ending position of the region.")
	var maxEnd *int = flag.Int("maxEnd", numbers.MaxInt, "Specifies the maximum ending position of the region.")
	var minNameFloat *float64 = flag.Float64("minNameFloat", -1 * math.MaxFloat64, "Specifies the minimum floating point number value for bed entries where floating point numbers are stored in the name field.")
	var maxNameFloat *float64 = flag.Float64("maxNameFloat", math.MaxFloat64, "Specifies the maximum floating point number value for bed entries where floating point numbers are stored in the name field.")
	var chrom *string = flag.String("chrom", "", "Specifies the chromosome name.")

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
		InFile: infile,
		OutFile: outfile,
		MinScore: *minScore,
		MinLength: *minLength,
		MaxScore: *maxScore,
		MaxLength: *maxLength,
		MinStart: *minStart,
		MaxStart: *maxStart,
		MinEnd: *minEnd,
		MaxEnd: *maxEnd,
		MinNameFloat: *minNameFloat,
		MaxNameFloat: *maxNameFloat,
		Chrom: *chrom,
	}
	bedFilter(s)
}

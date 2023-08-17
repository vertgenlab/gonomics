// Command Group: "BED Tools"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed/bedpe"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"log"
	"math"
)

type Settings struct {
	InFile         string
	OutFile        string
	MinScore       int
	MaxScore       int
	MinDistance    int
	MaxDistance    int
	MinStart       int
	MaxStart       int
	MinEnd         int
	MaxEnd         int
	Chrom          string
	OnlyInterChrom bool
	OnlyIntraChrom bool
	//SubSet         float64
}

// uses start positions of bedpe
func bedpeFilter(s Settings) {
	var distance int
	var pass = false
	if s.OnlyIntraChrom && s.OnlyInterChrom {
		log.Fatal("Error: Cannot have both Only Intra Chrom contacts and Only Inter Chrom contacts set to true.")
	}
	bedpeChan := bedpe.GoReadToChan(s.InFile)
	out := fileio.EasyCreate(s.OutFile)

	for curr := range bedpeChan {
		pass = true
		distance = int(math.Abs(float64(curr.A.ChromStart - curr.B.ChromStart)))
		if curr.A.Score < s.MinScore {
			pass = false
		}
		if curr.A.Score > s.MaxScore {
			pass = false
		}
		if distance < s.MinDistance {
			pass = false
		}
		if distance > s.MaxDistance {
			pass = false
		}
		if curr.A.ChromStart < s.MinStart || curr.B.ChromStart < s.MinStart {
			pass = false
		}
		if curr.A.ChromStart > s.MaxStart || curr.B.ChromStart > s.MaxStart {
			pass = false
		}
		if curr.A.ChromEnd < s.MinEnd || curr.B.ChromStart < s.MinEnd {
			pass = false
		}
		if curr.A.ChromEnd > s.MaxEnd || curr.B.ChromStart > s.MaxEnd {
			pass = false
		}
		if s.OnlyIntraChrom && curr.A.Chrom != curr.B.Chrom {
			pass = false
		}
		if s.OnlyInterChrom && curr.A.Chrom == curr.B.Chrom {
			pass = false
		}
		if curr.A.Chrom != s.Chrom && curr.B.Chrom != s.Chrom {
			pass = false
		}
		if pass {
			bedpe.WriteToFileHandle(out, curr)
		}
	}
	var err error
	err = out.Close()
	exception.PanicOnErr(err)
}

func usage() {
	fmt.Print(
		"bedpeFilter\n" +
			"Usage:\n" +
			"bedpeFilter input.bedpe output.bedpe\n" +
			"options:\n")
	flag.PrintDefaults()
}
func main() {
	var expectedNumArgs int = 2
	var minScore *int = flag.Int("minScore", -1*numbers.MaxInt, "Specifies the minimum score in the fourth field.")
	var maxScore *int = flag.Int("maxScore", numbers.MaxInt, "Specifies the maximum score in the fourth field.")
	var minDistance *int = flag.Int("minDistance", 0, "Specifies the minimum distance between the feet of the bedpe.")
	var maxDistance *int = flag.Int("maxDistance", numbers.MaxInt, "Specifies the maximum distance between the feet of the bedpe.")
	var minStart *int = flag.Int("minStart", 0, "Specifies the minimum starting position of either region such that if either region satisfies the requirement the whole contact will be in the outFile.")
	var maxStart *int = flag.Int("maxStart", numbers.MaxInt, "Specifies the maximum starting position of either region such that if either region satisfies the requirement the whole contact will be in the outFile.")
	var minEnd *int = flag.Int("minEnd", 0, "Specifies the minimum ending position of either region such that if either region satisfies the requirement the whole contact will be in the outFile.")
	var maxEnd *int = flag.Int("maxEnd", numbers.MaxInt, "Specifies the maximum ending position of either region such that if either region satisfies the requirement the whole contact will be in the outFile.")
	var chrom *string = flag.String("chrom", "", "Specifies the chromosome name of either region such that if either region satisfies the requirement the whole contact will be in the outFile.")
	var onlyInterChrom *bool = flag.Bool("onlyInterChrom", false, "When true the output file will only contains elements that describe contacts between chromosomes. When this and onlyIntraChrom options are both fase, chrom determines which contacts are kept.")
	var onlyIntraChrom *bool = flag.Bool("onlyIntraChrom", false, "When true the output file will only contains elements that describe contacts within the same chromosome. When this and onlyIntraChrom options are both fase, chrom determines which contacts are kept.")

	//var subSet *float64 = flag.Float64("subSet", 1.0, "Proportion of entries to retain in output, range from 0 to 1.")

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
		InFile:         infile,
		OutFile:        outfile,
		MinScore:       *minScore,
		MaxScore:       *maxScore,
		MinDistance:    *minDistance,
		MaxDistance:    *maxDistance,
		MinStart:       *minStart,
		MaxStart:       *maxStart,
		MinEnd:         *minEnd,
		MaxEnd:         *maxEnd,
		Chrom:          *chrom,
		OnlyInterChrom: *onlyInterChrom,
		OnlyIntraChrom: *onlyIntraChrom,
		//SubSet:                *subSet,
	}
	bedpeFilter(s)
}

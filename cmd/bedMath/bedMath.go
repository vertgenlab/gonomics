package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"strings"
)

type Operation byte

const (
	Add Operation = 0
	Subtract Operation = 1
	Multiply Operation = 2
	Divide Operation = 3
)

func bedMath(s Settings) {
	var err error
	var done = false
	var currA, currB bed.Bed
	var Op = parseOp(s.Op)
	aReader := fileio.EasyOpen(s.Afile)
	bReader := fileio.EasyOpen(s.Bfile)
	out := fileio.EasyCreate(s.OutFile)

	currA, done = bed.NextBed(aReader)
	if done {
		log.Fatalf("First bed file has no bed entries.")
	}
	currB, done = bed.NextBed(bReader)
	if done {
		log.Fatalf("Second bed file has no bed entries.")
	}

	for !done {
		if sameCoords(currA, currB) {
			bed.WriteBed(out, doBedMath(currA, currB, Op))
			currA, done = bed.NextBed(aReader)
		} else if bed.Compare(currA, currB) < 0 {//if a is less than b, we need another a
			currA, done = bed.NextBed(aReader)
		} else if bed.Compare(currA, currB) > 0 {//than a is greater than b, so we need another b
			currB, done = bed.NextBed(bReader)
		}
	}
	err = out.Close()
	exception.PanicOnErr(err)
}

func doBedMath(a bed.Bed, b bed.Bed, Op Operation) bed.Bed {
	var aFloat = common.StringToFloat64(a.Name)
	var bFloat = common.StringToFloat64(b.Name)
	if Op == Add {
		a.Name = fmt.Sprintf("%f", aFloat + bFloat)
	} else if Op == Subtract {
		a.Name = fmt.Sprintf("%f", aFloat - bFloat)
	} else if Op == Multiply {
		a.Name = fmt.Sprintf("%f", aFloat * bFloat)
	} else if Op == Divide {
		a.Name = fmt.Sprintf("%f", aFloat / bFloat)
	}
	return a
}

func sameCoords(a bed.Bed, b bed.Bed) bool {
	return a.Chrom == b.Chrom && a.ChromStart == b.ChromStart && a.ChromEnd == b.ChromEnd
}

func parseOp(op string) Operation {
	op = strings.ToLower(op)
	switch op {
	case "add":
		return Add
	case "plus":
		return Add
	case "subtract":
		return Subtract
	case "minus":
		return Subtract
	case "times":
		return Multiply
	case "multiply":
		return Multiply
	case "divide":
		return Divide
	case "divideby":
		return Divide
	default:
		log.Fatalf("Unrecognized operation: %v. Accepted operations are add, subtract, times, or divideBy.", op)
		return 0
	}
}

type Settings struct {
	Afile string
	Bfile string
	Op string
	OutFile string
}

func usage() {
	fmt.Print(
		"bedMath - Performs comparative arithmetic operations on float values in bed files.\n" +
			"Usage:\n" +
			"bedMath a.bed operation b.bed out.bed\n" +
			"operation may be one of the following options:" +
			"plus\tminus\ttimes\tdivideBy\n" +
			"Input bed files must be pre-sorted by coordinate.\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs = 4

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	aFile := flag.Arg(0)
	op := flag.Arg(1)
	bFile := flag.Arg(2)
	outFile := flag.Arg(4)

	s := Settings{
		Afile: aFile,
		Bfile: bFile,
		Op: op,
		OutFile: outFile,
	}

	bedMath(s)
}

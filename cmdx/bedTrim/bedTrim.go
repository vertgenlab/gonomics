package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/interval"
	"log"
	"math"
)

func bedTrim(trimPercent int, inBed string, outBed string) {
	var size, toRemoveInt int
	var toRemove float64

	bedChan := bed.GoReadToChan(inBed)
	out := fileio.EasyCreate(outBed)

	for b := range bedChan {
		size = interval.IntervalSize(b)
		toRemove = float64(size) * (float64(trimPercent) / 100)
		toRemoveInt = int(math.Round(toRemove))
		switch toRemoveInt % 2 {
		case 0:
			b.ChromStart += toRemoveInt / 2
			b.ChromEnd -= toRemoveInt / 2
		default:
			b.ChromStart += toRemoveInt/2 + 1
			b.ChromEnd -= toRemoveInt / 2
		}
		if interval.IntervalSize(b) > 0 {
			bed.WriteToFileHandle(out, b)
		}
	}
	exception.PanicOnErr(out.Close())
}

func usage() {
	fmt.Print("bedTrim -- Trim bed records in a file by a total of N percent. An equal amount of bases " +
		"will be removed from each side.\n" +
		"For example, a -trimPercent of 50 will remove 25% from each side of the bed region" +
		"Usage:\n" +
		"bedTrim -trimPercent [int] in.bed out.bed\n ")
	flag.PrintDefaults()
}

func main() {
	var trimPercent *int = flag.Int("trimPercent", 0, "Input a percentage for each bed record to be shortened by. "+
		"The bed region will be shorted by the same amount on either side ")

	var expectedNumArgs int = 2
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Expected %d arguments, got %d", expectedNumArgs, len(flag.Args()))
	}

	if *trimPercent < 0 || *trimPercent > 100 {
		flag.Usage()
		log.Fatalf("Error: trimPercent must be an integer between 0 and 100")
	}

	inBed := flag.Arg(0)
	outBed := flag.Arg(1)

	bedTrim(*trimPercent, inBed, outBed)

}

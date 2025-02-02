package main

import (
	"flag"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/interval"
	"github.com/vertgenlab/gonomics/numbers/parse"
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

func main() {
	var trimPercent *int = flag.Int("trimPercent", 0, "Input a percentage for each bed record to be shortened by. "+
		"The bed region will be shorted by the same amount on either side ")

	flag.Parse()
	if parse.StringToInt(flag.Arg(0)) < 0 || parse.StringToInt(flag.Arg(0)) > 100 {
		log.Fatalf("Error: trimPercent must be an integer between 0 and 100")
	}

	inBed := flag.Arg(1)
	outBed := flag.Arg(2)

	bedTrim(*trimPercent, inBed, outBed)

}

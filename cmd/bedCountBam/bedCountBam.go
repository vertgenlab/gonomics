package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/sam"
	"log"
)

type settings struct {
	inSam   string
	inBed   string
	outFile string
	norm    bool
}

func size(b bed.Bed) int {
	return numbers.AbsInt(b.ChromStart - b.ChromEnd)
}

func bedCountBam(s settings) {
	var ans []sam.Sam
	o := fileio.EasyCreate(s.outFile)
	if s.norm {
		fileio.WriteToFileHandle(o, "bedRegion\tcountsPerBP")
	} else {
		fileio.WriteToFileHandle(o, "bedRegion\tcounts")
	}
	bd := bed.Read(s.inBed)
	bm, _ := sam.OpenBam(s.inSam)
	bai := sam.ReadBai(s.inSam + ".bai")

	for i := range bd {
		ans = sam.SeekBamRegion(bm, bai, bd[i].Chrom, uint32(bd[i].ChromStart), uint32(bd[i].ChromEnd))
		if s.norm {
			fileio.WriteToFileHandle(o, fmt.Sprintf("%s\t%f", bd[i].Name, float64(len(ans))/float64(size(bd[i]))))
		} else {
			fileio.WriteToFileHandle(o, fmt.Sprintf("%s\t%d", bd[i].Name, len(ans)))
		}
	}
	err := o.Close()
	exception.PanicOnErr(err)
}

func usage() {

}

func main() {
	var norm *bool = flag.Bool("norm", false, "Report counts per base pairs instead of raw counts")

	var expectedNumArgs int = 3
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()
	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Expected %d arguments, got %d", expectedNumArgs, len(flag.Args()))
	}
	var s settings = settings{
		inSam:   flag.Arg(0),
		inBed:   flag.Arg(1),
		outFile: flag.Arg(2),
		norm:    *norm,
	}
	bedCountBam(s)
}

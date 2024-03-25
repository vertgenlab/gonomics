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
	inBam   string
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
	bm, _ := sam.OpenBam(s.inBam)
	bai := sam.ReadBai(s.inBam + ".bai")

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
	fmt.Print("bedCountBam -- read counts for bed regions from a position sorted and indexed bam file.\n" +
		"Bam file must be position sorted and indexed. Regions must have a unique identifier in the fourth field.\n" +
		"Usage:\n" +
		"bedCountBam [options] in.bam in.bed out.txt\n")
	flag.PrintDefaults()
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
		inBam:   flag.Arg(0),
		inBed:   flag.Arg(1),
		outFile: flag.Arg(2),
		norm:    *norm,
	}
	bedCountBam(s)
}

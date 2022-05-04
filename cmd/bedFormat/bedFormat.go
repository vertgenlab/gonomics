// Command Group: "BED Tools"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/convert"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
)

func bedFormat(s Settings) {
	var err error
	ch := bed.GoReadToChan(s.InFile)
	out := fileio.EasyCreate(s.OutFile)

	if s.EnsemblToUCSC && s.UCSCToEnsembl {
		log.Fatalf("Both conversions (UCSCToEnsembl and EnsemblToUCSC) are incompatable.")
	}

	for v := range ch {
		if s.Pad > 0 {
			v.ChromStart -= s.Pad
			v.ChromEnd += s.Pad
			if v.ChromStart < 0 {
				v.ChromStart = 0
			}
		}
		if s.EnsemblToUCSC {
			v.Chrom = convert.EnsemblToUCSC(v.Chrom)
		}
		if s.UCSCToEnsembl {
			v.Chrom = convert.UCSCToEnsembl(v.Chrom)
		}
		if s.ScaleNameFloat != 1 {
			v.Name = fmt.Sprintf("%.8g", s.ScaleNameFloat*common.StringToFloat64(v.Name))
		}
		bed.WriteBed(out, v)
	}

	err = out.Close()
	exception.PanicOnErr(err)
}

func usage() {
	fmt.Print(
		"bedFormat - Options alter bed formatting.\n" +
			"Usage:\n" +
			"bedFormat input.bed output.bed\n" +
			"options:\n")
	flag.PrintDefaults()
}

type Settings struct {
	InFile         string
	OutFile        string
	UCSCToEnsembl  bool
	EnsemblToUCSC  bool
	ScaleNameFloat float64
	Pad            int
}

func main() {
	var expectedNumArgs int = 2
	var ensemblToUCSC *bool = flag.Bool("ensemblToUCSC", false, "Changes chromosome format type.")
	var UCSCToEnsembl *bool = flag.Bool("UCSCToEnsembl", false, "Changes chromosome format type.")
	var scaleNameFloat *float64 = flag.Float64("scaleNameFloat", 1, "If float values are held in the name field, scale those values by this constant multiplier.")
	var pad *int = flag.Int("pad", 0, "Add # of bases to both ends of each bed record")

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
		UCSCToEnsembl:  *UCSCToEnsembl,
		EnsemblToUCSC:  *ensemblToUCSC,
		ScaleNameFloat: *scaleNameFloat,
		Pad:            *pad,
	}

	bedFormat(s)
}

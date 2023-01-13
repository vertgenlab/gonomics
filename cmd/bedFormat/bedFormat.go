// Command Group: "BED Tools"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/convert"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"log"
)

type Settings struct {
	InFile         string
	OutFile        string
	UCSCToEnsembl  bool
	EnsemblToUCSC  bool
	ScaleNameFloat float64
	EvenPadLength  int
	UpstreamPadLength int
	DownstreamPadLength int
	ChromSizeFile  string
	ToMidpoint     bool
	ToTss bool
}

func bedFormat(s Settings) {
	var err error
	var inMap bool
	var sizes map[string]chromInfo.ChromInfo
	ch := bed.GoReadToChan(s.InFile)
	out := fileio.EasyCreate(s.OutFile)

	if s.EnsemblToUCSC && s.UCSCToEnsembl {
		log.Fatalf("Both conversions (UCSCToEnsembl and EnsemblToUCSC) are incompatable.")
	}
	if s.ChromSizeFile == "" && (s.EvenPadLength > 0 || s.UpstreamPadLength > 0 || s.DownstreamPadLength > 0) {
		log.Fatalf("Must specify a chromFile to use a padLength option.")
	}
	if s.ChromSizeFile != "" && (s.EvenPadLength > 0 || s.UpstreamPadLength > 0 || s.DownstreamPadLength > 0) {
		sizes = chromInfo.ReadToMap(s.ChromSizeFile)
	}
	if s.ToTss && s.ToMidpoint {
		log.Fatalf("Cannot trim bed elements to midpoint AND to Tss.")
	}

	for v := range ch {
		if s.ToMidpoint {
			v = bed.ToMidpoint(v)
		}
		if s.ToTss {
			v = bed.ToTss(v)
		}
		if s.EvenPadLength > 0 {
			if _, inMap = sizes[v.Chrom]; !inMap {
				log.Fatalf("Chrom for current bed entry not found in chromSizes file. BedChrom: %s.", v.Chrom)
			}
			v.ChromStart = numbers.Max(v.ChromStart-s.EvenPadLength, 0)
			v.ChromEnd = numbers.Min(v.ChromEnd+s.EvenPadLength, sizes[v.Chrom].Size)
		}
		if s.UpstreamPadLength > 0 {
			if _, inMap = sizes[v.Chrom]; !inMap {
				log.Fatalf("Chrom for current bed entry not found in chromSizes file. BedChrom: %s.", v.Chrom)
			}
			switch v.Strand {
			case bed.Positive:
				v.ChromStart = numbers.Max(v.ChromStart - s.UpstreamPadLength, 0)
			case bed.Negative:
				v.ChromEnd = numbers.Min(v.ChromEnd + s.UpstreamPadLength, sizes[v.Chrom].Size)
			default:
				log.Fatalf("Bed entries must have annotated strand information to perform upstream padding.")
			}
		}
		if s.DownstreamPadLength > 0 {
			if _, inMap = sizes[v.Chrom]; !inMap {
				log.Fatalf("Chrom for current bed entry not found in chromSizes file. BedChrom: %s.", v.Chrom)
			}
			switch v.Strand {
			case bed.Positive:
				v.ChromEnd = numbers.Min(v.ChromEnd + s.DownstreamPadLength, sizes[v.Chrom].Size)
			case bed.Negative:
				v.ChromStart = numbers.Max(v.ChromStart - s.DownstreamPadLength, 0)
			default:
				log.Fatalf("Bed entries must have annotated strand information to perform downstream padding.")
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

func main() {
	var expectedNumArgs int = 2
	var evenPadLength *int = flag.Int("evenPadLength", 0, "Add # of bases to both ends of each bed record. Requires chromSizeFile.")
	var upstreamPadLength *int = flag.Int("upstreamPadLength", 0, "Add # of bases upstream of a bed record, strand-sensitive. Requires chromSizeFile.")
	var downstreamPadLength *int = flag.Int("downstreamPadLength", 0, "Add # of bases downstream of a bed record, strand-sensitive. Requires chromSizeFile.")
	var ensemblToUCSC *bool = flag.Bool("ensemblToUCSC", false, "Changes chromosome format type.")
	var UCSCToEnsembl *bool = flag.Bool("UCSCToEnsembl", false, "Changes chromosome format type.")
	var scaleNameFloat *float64 = flag.Float64("scaleNameFloat", 1, "If float values are held in the name field, scale those values by this constant multiplier.")
	var chromSizeFile *string = flag.String("chromSizeFile", "", "Specify a .chrom.sizes file for use with the padLength option. Ensures padding is truncated at chromosome ends.")
	var ToMidpoint *bool = flag.Bool("ToMidpoint", false, "Trim the output bed to single-base pair ranges at the midpoint of the input bed ranges.")
	var ToTss *bool = flag.Bool("ToTss", false, "Trim the output bed to a single-base pair range at the start of the region. Strand-sensitive.")

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
		EvenPadLength:  *evenPadLength,
		UpstreamPadLength: *upstreamPadLength,
		DownstreamPadLength: *downstreamPadLength,
		ChromSizeFile:  *chromSizeFile,
		ToMidpoint:     *ToMidpoint,
		ToTss: *ToTss,
	}

	bedFormat(s)
}

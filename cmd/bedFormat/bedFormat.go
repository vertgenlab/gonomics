// Command Group: "BED Tools"

// Options to alter bed formatting
package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/convert"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/numbers/parse"
	"log"
	"math"
	"os"
)

type Settings struct {
	InFile                   string
	OutFile                  string
	UCSCToEnsembl            bool
	EnsemblToUCSC            bool
	ScaleNameFloat           float64
	EvenPadLength            int
	UpstreamPadLength        int
	DownstreamPadLength      int
	ChromSizeFile            string
	ToMidpoint               bool
	ToTss                    bool
	FdrScoreAnnotation       bool
	RawPValueAnnotationField int
	ScoreBuffer              int
}

// FdrConverter contains a RawPValue and Rank, which are used to calculate an AdjPValue
type FdrConverter struct {
	Count     int
	RawPValue float64
	Rank      int
	AdjPValue float64
}

func bedFormat(s Settings) {
	var err error
	var inMap bool
	var tmp *fileio.EasyWriter
	var sizes map[string]chromInfo.ChromInfo
	ch := bed.GoReadToChan(s.InFile)

	fdrList := make([]FdrConverter, s.ScoreBuffer)
	var totalWindows int = 0

	if s.FdrScoreAnnotation {
		tmp = fileio.EasyCreate(s.OutFile + ".tmp")
	}

	out := fileio.EasyCreate(s.OutFile)

	if s.EnsemblToUCSC && s.UCSCToEnsembl {
		log.Fatalf("Both conversions (UCSCToEnsembl and EnsemblToUCSC) are incompatible.")
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
				v.ChromStart = numbers.Max(v.ChromStart-s.UpstreamPadLength, 0)
			case bed.Negative:
				v.ChromEnd = numbers.Min(v.ChromEnd+s.UpstreamPadLength, sizes[v.Chrom].Size)
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
				v.ChromEnd = numbers.Min(v.ChromEnd+s.DownstreamPadLength, sizes[v.Chrom].Size)
			case bed.Negative:
				v.ChromStart = numbers.Max(v.ChromStart-s.DownstreamPadLength, 0)
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
			v.Name = fmt.Sprintf("%.8g", s.ScaleNameFloat*parse.StringToFloat64(v.Name))
		}

		if s.FdrScoreAnnotation {
			totalWindows++
			if v.Score > s.ScoreBuffer {
				log.Fatalf("Error: score in bed entry: %v, is greater than the scoreBuffer. Raise the scoreBuffer and rerun.", v.Score)
			}
			fdrList[v.Score].Count++
			if s.RawPValueAnnotationField >= len(v.Annotation) {
				log.Fatalf("Error: rawPValueAnnotationField, %v, exceeds the length of the annotation slice in bed entry: %v.\n", s.RawPValueAnnotationField, len(v.Annotation))
			}
			fdrList[v.Score].RawPValue = parse.StringToFloat64(v.Annotation[s.RawPValueAnnotationField])
			bed.WriteBed(tmp, v)
		} else {
			bed.WriteBed(out, v)
		}
	}

	if s.FdrScoreAnnotation {
		err = tmp.Close()
		exception.PanicOnErr(err)

		var currRank int = 0
		for i := len(fdrList) - 1; i >= 0; i-- {
			currRank += fdrList[i].Count
			fdrList[i].Rank = currRank
			fdrList[i].RawPValue = math.Pow(10, -1*fdrList[i].RawPValue) //leaving logSpace
			fdrList[i].AdjPValue = math.Max(0, -1*math.Log10(fdrList[i].RawPValue*(float64(totalWindows)/float64(fdrList[i].Rank))))
			//fdrList[i].AdjPValue = numbers.Max(logspace.Multiply(math.Log10(float64(totalWindows)/float64(fdrList[i].Rank)), fdrList[i].RawPValue), 0)
		}

		ch = bed.GoReadToChan(s.OutFile + ".tmp")
		for v := range ch {
			v.Annotation = append(v.Annotation, fmt.Sprintf("%e", fdrList[v.Score].AdjPValue))
			bed.WriteBed(out, v)
		}

		err = os.Remove(s.OutFile + ".tmp")
		exception.PanicOnErr(err)
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
	var FdrScoreAnnotation *bool = flag.Bool("fdrScoreAnnotation", false, "Used when an annotation field stores a raw P value related to a numerical value in the score column (such as in faFindFast). "+
		"Appends an FDR-adjusted P values to the first free annotation column.")
	var rawPValueAnnotationField *int = flag.Int("rawPValueAnnotationField", 0, "Specify the annotation field where raw P values are stored for fdrScoreAnnotation.")
	var scoreBuffer *int = flag.Int("scoreBuffer", 500, "Set the size of the cache for FDR calculations. Should be greater than or equal to max score in the input bed file.")

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
		InFile:                   infile,
		OutFile:                  outfile,
		UCSCToEnsembl:            *UCSCToEnsembl,
		EnsemblToUCSC:            *ensemblToUCSC,
		ScaleNameFloat:           *scaleNameFloat,
		EvenPadLength:            *evenPadLength,
		UpstreamPadLength:        *upstreamPadLength,
		DownstreamPadLength:      *downstreamPadLength,
		ChromSizeFile:            *chromSizeFile,
		ToMidpoint:               *ToMidpoint,
		ToTss:                    *ToTss,
		FdrScoreAnnotation:       *FdrScoreAnnotation,
		RawPValueAnnotationField: *rawPValueAnnotationField,
		ScoreBuffer:              *scoreBuffer,
	}

	bedFormat(s)
}

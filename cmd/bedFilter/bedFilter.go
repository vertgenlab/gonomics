// Command Group: "BED Tools"

// Output a subset of a bed file using score, name, position, and length
package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/numbers/parse"
	"log"
	"math"
	"math/rand"
)

func bedFilter(s Settings) {
	rand.Seed(s.SetSeed)
	var length int
	var pass bool = false
	var r float64
	bedChan := bed.GoReadToChan(s.InFile)
	out := fileio.EasyCreate(s.OutFile)

	for curr := range bedChan {
		pass = true
		length = curr.ChromEnd - curr.ChromStart
		if curr.FieldsInitialized > 4 {
			if curr.Score < s.MinScore {
				pass = false
			}
			if curr.Score > s.MaxScore {
				pass = false
			}
		} else if s.MinScore != -1*numbers.MaxInt || s.MaxScore != numbers.MaxInt { //if the scores are not the default options but the entry has no score, the entry will not be retained.
			pass = false
		}
		if length < s.MinLength {
			pass = false
		}
		if length > s.MaxLength {
			pass = false
		}
		if curr.ChromStart < s.MinStart {
			pass = false
		}
		if curr.ChromStart > s.MaxStart {
			pass = false
		}
		if curr.ChromEnd < s.MinEnd {
			pass = false
		}
		if curr.ChromEnd > s.MaxEnd {
			pass = false
		}
		if s.MinNameFloat != -1*math.MaxFloat64 {
			if parse.StringToFloat64(curr.Name) < s.MinNameFloat {
				pass = false
			}
		}
		if s.MaxNameFloat != math.MaxFloat64 {
			if parse.StringToFloat64(curr.Name) > s.MaxNameFloat {
				pass = false
			}
		}
		if s.MinAnnotationFloat != -1*math.MaxFloat64 {
			if s.AnnotationFilterField >= len(curr.Annotation) {
				log.Fatalf("Error: specified annotationFilterField that exceeds the number of annotation fields in the input bed. AnnotationFilterField: %v. LenCurrAnnotation: %v. CurrBed: %v.", s.AnnotationFilterField, len(curr.Annotation), bed.ToString(curr, curr.FieldsInitialized))
			}
			if parse.StringToFloat64(curr.Annotation[s.AnnotationFilterField]) < s.MinAnnotationFloat {
				pass = false
			}
		}
		if s.MaxAnnotationFloat != math.MaxFloat64 {
			if s.AnnotationFilterField >= len(curr.Annotation) {
				log.Fatalf("Error: specified annotationFilterField that exceeds the number of annotation fields in the input bed. AnnotationFilterField: %v. LenCurrAnnotation: %v. CurrBed: %v.", s.AnnotationFilterField, len(curr.Annotation), bed.ToString(curr, curr.FieldsInitialized))
			}
			if parse.StringToFloat64(curr.Annotation[s.AnnotationFilterField]) > s.MaxAnnotationFloat {
				pass = false
			}
		}
		if s.Chrom != "" {
			if curr.Chrom != s.Chrom {
				pass = false
			}
		}
		if s.NameEquals != "" {
			if curr.Name != s.NameEquals {
				pass = false
			}
		}
		if s.NameNotEquals != "" {
			if curr.Name == s.NameNotEquals {
				pass = false
			}
		}
		if pass && s.SubSet < 1.0 {
			r = rand.Float64()
			if r > s.SubSet {
				pass = false
			}
		}
		if pass {
			bed.WriteBed(out, curr)
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

type Settings struct {
	InFile                string
	OutFile               string
	MinScore              int
	MaxScore              int
	MinLength             int
	MaxLength             int
	MinStart              int
	MaxStart              int
	MinEnd                int
	MaxEnd                int
	MinNameFloat          float64
	MaxNameFloat          float64
	NameEquals            string
	NameNotEquals         string
	MinAnnotationFloat    float64
	MaxAnnotationFloat    float64
	AnnotationFilterField int
	Chrom                 string
	SubSet                float64
	SetSeed               int64
}

func main() {
	var expectedNumArgs int = 2
	var minScore *int = flag.Int("minScore", -1*numbers.MaxInt, "Specifies the minimum score in the fourth field.")
	var maxScore *int = flag.Int("maxScore", numbers.MaxInt, "Specifies the maximum score in the fourth field.")
	var minLength *int = flag.Int("minLength", 0, "Specifies the minimum length of the region.")
	var maxLength *int = flag.Int("maxLength", numbers.MaxInt, "Specifies the maximum length of the region.")
	var minStart *int = flag.Int("minStart", 0, "Specifies the minimum starting position of the region.")
	var maxStart *int = flag.Int("maxStart", numbers.MaxInt, "Specifies the maximum starting position of the region.")
	var minEnd *int = flag.Int("minEnd", 0, "Specifies the minimum ending position of the region.")
	var maxEnd *int = flag.Int("maxEnd", numbers.MaxInt, "Specifies the maximum ending position of the region.")
	var minNameFloat *float64 = flag.Float64("minNameFloat", -1*math.MaxFloat64, "Specifies the minimum floating point number value for bed entries where floating point numbers are stored in the name field.")
	var maxNameFloat *float64 = flag.Float64("maxNameFloat", math.MaxFloat64, "Specifies the maximum floating point number value for bed entries where floating point numbers are stored in the name field.")
	var nameEquals *string = flag.String("nameEquals", "", "Returns all bed entries with a name field that matches the input.")
	var nameNotEquals *string = flag.String("nameNotEquals", "", "Returns all bed entries with a name field that doesn't match the input.")
	var minAnnotationFloat *float64 = flag.Float64("minAnnotationFloat", -1*math.MaxFloat64, "Specifies the minimum floating point number value for bed entries where floating point numbers are stored in the an annotation field. Annotation field for filtering can be specified with -annotationFilterField.")
	var maxAnnotationFloat *float64 = flag.Float64("maxAnnotationFloat", math.MaxFloat64, "Specifies the maximum floating point number value for bed entries where floating point numbers are stored in the an annotation field. Annotation field for filtering can be specified with -annotationFilterField.")
	var annotationFilterField *int = flag.Int("annotationFilterField", 0, "Specify which annotation field (0-based) will be used for filtering with min/maxAnnotationFloat.")
	var chrom *string = flag.String("chrom", "", "Specifies the chromosome name.")
	var subSet *float64 = flag.Float64("subSet", 1.0, "Proportion of entries to retain in output, range from 0 to 1.")
	var setSeed *int64 = flag.Int64("setSeed", -1, "Use a specific seed for the RNG.")

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
		InFile:                infile,
		OutFile:               outfile,
		MinScore:              *minScore,
		MaxScore:              *maxScore,
		MinLength:             *minLength,
		MaxLength:             *maxLength,
		MinStart:              *minStart,
		MaxStart:              *maxStart,
		MinEnd:                *minEnd,
		MaxEnd:                *maxEnd,
		MinNameFloat:          *minNameFloat,
		MaxNameFloat:          *maxNameFloat,
		NameEquals:            *nameEquals,
		NameNotEquals:         *nameNotEquals,
		MinAnnotationFloat:    *minAnnotationFloat,
		MaxAnnotationFloat:    *maxAnnotationFloat,
		AnnotationFilterField: *annotationFilterField,
		Chrom:                 *chrom,
		SubSet:                *subSet,
		SetSeed:               *setSeed,
	}
	bedFilter(s)
}

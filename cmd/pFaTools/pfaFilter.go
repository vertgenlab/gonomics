// Returns a filtered pfasta based on option parameters
package main

import (
	"flag"
	"fmt"
	"log"
	"path"
	"strings"

	"github.com/vertgenlab/gonomics/dna/pDna"
	"github.com/vertgenlab/gonomics/pFasta"
)

type settings struct {
	InFile         string
	OutFile        string
	RefPositions   bool
	Start          int
	End            int
	Name           string
	NotName        string
	NameContains   string
	MinSize        int
	MaxGC          float64
	MinGC          float64
	FinalNBases    int
	CutFinalNBases int
	AppendBefore   string
	AppendAfter    string
}

func appendSeq(s settings, outlist []pFasta.PFasta) []pFasta.PFasta {
	var toAppendBefore, toAppendAfter []pFasta.PFasta
	if path.Ext(s.AppendBefore) == ".fa" {
		toAppendBefore = pFasta.Read(s.AppendBefore)
	} else {
		toAppendBefore = []pFasta.PFasta{{Seq: pDna.StringToBases(s.AppendBefore)}}
	}
	if path.Ext(s.AppendAfter) == ".fa" {
		toAppendAfter = pFasta.Read(s.AppendAfter)
	} else {
		toAppendAfter = []pFasta.PFasta{{Seq: pDna.StringToBases(s.AppendAfter)}}
	}

	if len(toAppendBefore) != 1 || len(toAppendAfter) != 1 {
		log.Fatalf("ERROR: Fasta file for appending must only have 1 fasta record")
	}

	for i := range outlist {
		outlist[i].Seq = append(toAppendBefore[0].Seq, outlist[i].Seq...)
		outlist[i].Seq = append(outlist[i].Seq, toAppendAfter[0].Seq...)
	}
	return outlist
}

func faFilter(s settings) {
	records := pFasta.Read(s.InFile) //read the fasta infile
	var length int
	var outlist []pFasta.PFasta //make the variable to store the fasta records that will be written out
	var pass bool = true

	if s.Start > s.End && s.End != -1 { //throws an error if the start and end positions given as options (see more info below) are illogical
		log.Fatalf("End must be larger than Start.")
	}

	if s.RefPositions { //if the user would like the output to be written in the context of the reference genome rather than counting any gaps that could exist from a multiFasta file
		s.Start = pFasta.RefPosToAlnPos(records[0], s.Start) //adjusts the start position of the fasta record according to user's input
		s.End = pFasta.RefPosToAlnPos(records[0], s.End)     //adjusts end coodinates
	}

	for i := 0; i < len(records); i++ { //for each fasta record
		pass = true
		if s.Name != "" && records[i].Name != s.Name { //filtering based on name if a string is given to the name option by the user
			pass = false
		}
		if s.NotName != "" && records[i].Name == s.NotName { //filtering out records based on name similar to the above option, name must not match option entry
			pass = false
		}
		if s.NameContains != "" && !strings.Contains(records[i].Name, s.NameContains) {
			pass = false
		}
		if len(records[i].Seq) < s.MinSize { //filtering based on the size of the fasta record, must be larger than the value given to the minSize option
			pass = false
		}
		if pDna.GCContent(records[i].Seq) > s.MaxGC { //filtering based on GC content of the fasta record; GC content must be less than or equal to the value given to the maxGC option
			pass = false
		}
		if pDna.GCContent(records[i].Seq) < s.MinGC { //filtering based on GC content of the fasta record; GC content must be greater than or equal to the value given to the minGC option
			pass = false
		}
		if pass { //if checks passed on a record
			if s.FinalNBases > 0 {
				length = len(records[i].Seq)
				if s.FinalNBases > length {
					length = s.FinalNBases
				}
				records[i].Seq = records[i].Seq[length-s.FinalNBases:]
			} else if s.CutFinalNBases > 0 {
				length = len(records[i].Seq)
				if s.CutFinalNBases >= length {
					continue
				}
				records[i].Seq = records[i].Seq[:length-s.CutFinalNBases]
			} else {
				if s.End == -1 { //if the user didn't ask the record to stop at a specific location, append the fasta until the end of the record
					records[i].Seq = records[i].Seq[s.Start:]
				} else { //if an endis specified by the user, append from the start to the finish specified
					records[i].Seq = records[i].Seq[s.Start:s.End]
				}
			}
			outlist = append(outlist, records[i]) //write any records to the outlist
		}
	}
	if s.AppendBefore != "" || s.AppendAfter != "" {
		outlist = appendSeq(s, outlist)
	}
	pFasta.Write(s.OutFile, outlist) //write the outlist to a file
}

func usage() {
	fmt.Print(
		"faFilter - Returns a filtered fasta based on option parameters.\n" +
			"Usage:\n" +
			"faFilter input.fa output.fa\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	//here are examples of using options in a command. They have default values, which will depend on their type and the desired result. Bools can be set to true or false, ints can be any number,
	//but if there's a desired "missing value" number, then that may be best as the default as shown with the "end " option here. Strings defaults are usually empty because they refer to files that aren't needed in the defaults, typically.
	var refPositions *bool = flag.Bool("refPositions", false, "Uses reference positions for range specifications instead of alignment positions.")
	var start *int = flag.Int("start", 0, "Retains the sequence after this position.")
	var end *int = flag.Int("end", -1, "Retains the sequence before this position.")
	var name *string = flag.String("name", "", "Specifies the fasta record name.")
	var notName *string = flag.String("notName", "", "Returns all fasta records except for this input.")
	var nameContains *string = flag.String("nameContains", "", "Returns all fasta records whose name contains this input. String matching is case-sensitive.")
	var minSize *int = flag.Int("minSize", 0, "Retains all fasta records with a sequence of at least that size")
	var maxGC *float64 = flag.Float64("maxGC", 100, "Retains all fasta records with GC content less than or equal to this percentage.")
	var minGC *float64 = flag.Float64("minGC", 0, "Retains all fasta records with GC content greater than or equal to this percentage")
	var finalNBases *int = flag.Int("finalNBases", -1, "Retains the final N bases in the fasta record. Not compatible with -start or -end")
	var cutFinalNBases *int = flag.Int("cutFinalNbases", -1, "cuts the final N bases from each fasta record. Not compatible with -finalNbases, -start or -end")
	var appendBefore *string = flag.String("appendBefore", "", "Provide a fasta file with 1 sequence which will be appended in front of all fasta records in the input file. The append step will happen after any filtering steps")
	var appendAfter *string = flag.String("appendAfter", "", "Provide a fasta file with 1 sequence which will be appended atfer all fasta records in the input file. The append step will happen after any filtering steps")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	if (*finalNBases > 0 || *cutFinalNBases > 0) && (*start > 0 || *end > 0) {
		flag.Usage()
		log.Fatalf("Error: -finalNbases/-cutFinalNbases and -start/-end are not compatible with each other.")
	}

	if *finalNBases > 0 && *cutFinalNBases > 0 {
		flag.Usage()
		log.Fatalf("Error: -finalNbases and -cutFinalNbases are not compatible with each other.")
	}

	//all options exist as pointers in the arguments of the function call, since they may or may not exist when the user calls the function.
	var s settings = settings{
		InFile:         flag.Arg(0),
		OutFile:        flag.Arg(1),
		Name:           *name,
		NameContains:   *nameContains,
		NotName:        *notName,
		Start:          *start,
		End:            *end,
		MinSize:        *minSize,
		MinGC:          *minGC,
		MaxGC:          *maxGC,
		FinalNBases:    *finalNBases,
		CutFinalNBases: *cutFinalNBases,
		RefPositions:   *refPositions,
		AppendAfter:    *appendAfter,
		AppendBefore:   *appendBefore,
	}

	faFilter(s)
}

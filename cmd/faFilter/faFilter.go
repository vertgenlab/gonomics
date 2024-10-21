// Command Group: "FASTA and Multi-FASTA Tools"

// Returns a filtered fasta based on option parameters
package main

import (
	"flag"
	"fmt"
	"log"
	"strings"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
)

func faFilter(infile string, outfile string, name string, notName string, nameContains string, refPositions bool, start int, end int, minSize int, maxGC float64, minGC float64, finalNBases int, cutFinalNbases int) {
	records := fasta.Read(infile) //read the fasta infile
	var length int
	var outlist []fasta.Fasta //make the variable to store the fasta records that will be written out
	var pass bool = true

	if start > end && end != -1 { //throws an error if the start and end positions given as options (see more info below) are illogical
		log.Fatalf("End must be larger than Start.")
	}

	if refPositions { //if the user would like the output to be written in the context of the reference genome rather than counting any gaps that could exist from a multiFasta file
		start = fasta.RefPosToAlnPos(records[0], start) //adjusts the start position of the fasta record according to user's input
		end = fasta.RefPosToAlnPos(records[0], end)     //adjusts end coodinates
	}

	for i := 0; i < len(records); i++ { //for each fasta record
		pass = true
		if name != "" && records[i].Name != name { //filtering based on name if a string is given to the name option by the user
			pass = false
		}
		if notName != "" && records[i].Name == notName { //filtering out records based on name similar to the above option, name must not match option entry
			pass = false
		}
		if nameContains != "" && !strings.Contains(records[i].Name, nameContains) {
			pass = false
		}
		if len(records[i].Seq) < minSize { //filtering based on the size of the fasta record, must be larger than the value given to the minSize option
			pass = false
		}
		if dna.GCContent(records[i].Seq) > maxGC { //filtering based on GC content of the fasta record; GC content must be less than or equal to the value given to the maxGC option
			pass = false
		}
		if dna.GCContent(records[i].Seq) < minGC { //filtering based on GC content of the fasta record; GC content must be greater than or equal to the value given to the minGC option
			pass = false
		}
		if pass { //if checks passed on a record
			if finalNBases > 0 {
				length = len(records[i].Seq)
				if finalNBases > length {
					length = finalNBases
				}
				records[i].Seq = records[i].Seq[length-finalNBases:]
			} else if cutFinalNbases > 0 {
				length = len(records[i].Seq)
				if cutFinalNbases >= length {
					continue
				}
				records[i].Seq = records[i].Seq[:length-cutFinalNbases]
			} else {
				if end == -1 { //if the user didn't ask the record to stop at a specific location, append the fasta until the end of the record
					records[i].Seq = records[i].Seq[start:]
				} else { //if an endis specified by the user, append from the start to the finish specified
					records[i].Seq = records[i].Seq[start:end]
				}
			}
			outlist = append(outlist, records[i]) //write any records to the outlist
		}
	}
	fasta.Write(outfile, outlist) //write the outlist to a file
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

	inFile := flag.Arg(0)
	outFile := flag.Arg(1)

	faFilter(inFile, outFile, *name, *notName, *nameContains, *refPositions, *start, *end, *minSize, *maxGC, *minGC, *finalNBases, *cutFinalNBases) //all options exist as pointers in the arguments of the function call, since they may or may not exist when the user calls the function.
}

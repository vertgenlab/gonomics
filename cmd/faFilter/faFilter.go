// Command Group: "FASTA and Multi-FASTA Tools"

package main

import (
	"flag"
	"fmt"
	"log"

	"github.com/vertgenlab/gonomics/fasta"
)

func faFilter(infile string, outfile string, name string, notName string, refPositions bool, start int, end int, minSize int) {
	records := fasta.Read(infile) //read the fasta infile
	var outlist []fasta.Fasta     //make the variable to store the fasta records that will be written out
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
		if len(records[i].Seq) < minSize { //filtering based on the size of the fasta record, must be larger than the value given to the minSize option
			pass = false
		}
		if pass { //if checks passed on a record
			if end == -1 { //if the user didn't ask the record to stop at a specific location, append the fasta until the end of the record
				records[i].Seq = records[i].Seq[start:]
			} else { //if an endis specified by the user, append from the start to the finish specified
				records[i].Seq = records[i].Seq[start:end]
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
	//here are examples of using options in a command. They have default values, which will depend on their type and the desired result. Bools can be set to tru or false, ints can be any number,
	//but if there's a desired "missing value" number, then that may be best as the default as shown with the "end " option here. Strings defaults are usually empty because they refer to files that aren't needed in the defaults, typically.
	var refPositions *bool = flag.Bool("refPositions", false, "Uses reference positions for range specifications instead of alignment positions.")
	var start *int = flag.Int("start", 0, "Retains the sequence after this position.")
	var end *int = flag.Int("end", -1, "Retains the sequence before this position.")
	var name *string = flag.String("name", "", "Specifies the fasta record name.")
	var notName *string = flag.String("notName", "", "Returns all fasta records except for this input.")
	var minSize *int = flag.Int("minSize", 0, "Retains all fasta records with a sequence of at least that size")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)
	outFile := flag.Arg(1)

	faFilter(inFile, outFile, *name, *notName, *refPositions, *start, *end, *minSize) //all options exist as pointers in the arguments of the function call, since they may or may not exist when the user calls the function.
}

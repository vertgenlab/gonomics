// Command Group: "FASTA and Multi-FASTA Tools"

// Scan multiple Fasta alignment for a user-specified pattern ('present bases (A,C,G,T, not gap- or N)' for now) and report counts
package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
)

type Settings struct {
	InFile          string
	OutFile         string
	QueryName       string
	Both            bool
	Either          bool
	SecondQueryName string
	ThirdQueryName  string
}

func multiFaCount(s Settings) {
	aln := fasta.Read(s.InFile)
	var presentBaseCount int
	var err error

	// Initialize output file
	out := fileio.EasyCreate(s.OutFile)

	if s.Both && s.Either { // illegal
		log.Fatalf("both and either modes were both activated, but they are currently mutually exclusive")
	} else if !s.Both && !s.Either { // regular mode
		_, err = fmt.Fprintf(out, "#querySequenceName\tpresentBaseCount\n")
		exception.PanicOnErr(err)

		presentBaseCount = fasta.ScanPresentBase(aln, s.QueryName)

		_, err = fmt.Fprintf(out, "%s\t%d\n", s.QueryName, presentBaseCount)
		exception.PanicOnErr(err)
	} else if s.Both { // (&& !s.Either) both mode
		_, err = fmt.Fprintf(out, "#firstQuerySequenceName\tsecondQuerySequenceName\tbothPresentBaseCount\n")
		exception.PanicOnErr(err)

		presentBaseCount = fasta.ScanPresentBaseBoth(aln, s.QueryName, s.SecondQueryName)

		_, err = fmt.Fprintf(out, "%s\t%s\t%d\n", s.QueryName, s.SecondQueryName, presentBaseCount)
		exception.PanicOnErr(err)
	} else { // (the only possibility left is !s.Both && s.Either) either mode
		_, err = fmt.Fprintf(out, "#firstQuerySequenceName\tsecondQuerySequenceName\tthirdQuerySequenceName\teitherPresentBaseCount\n")
		exception.PanicOnErr(err)

		presentBaseCount = fasta.ScanPresentBaseEither(aln, s.QueryName, s.SecondQueryName, s.ThirdQueryName)

		_, err = fmt.Fprintf(out, "%s\t%s\t%s\t%d\n", s.QueryName, s.SecondQueryName, s.ThirdQueryName, presentBaseCount)
		exception.PanicOnErr(err)
	}

	// close output file
	err = out.Close()
	exception.PanicOnErr(err)
}

func usage() {
	fmt.Print(
		"multiFaCount - Scan multiple Fasta alignment for a user-specified pattern ('present bases (A,C,G,T, not gap- or N)' for now) and report counts.\n" +
			"queryName is set to the reference sequence by default\n" +
			"Usage:\n" +
			"multiFaCount queryName multi.fa\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 3
	var both *bool = flag.Bool("both", false, "When set to true, activates both mode. Count the number of positions where 2 sequences both match the user-specified pattern")
	var either *bool = flag.Bool("either", false, "When set to true, activates either mode. Count the number of positions where any of the specified sequences matches the user-specified pattern")
	var secondQueryName *string = flag.String("secondQueryName", "", "Specify the name of the second sequence to count for the user-specified pattern. Required in both and either modes. Set to the reference sequence by default")
	var thirdQueryName *string = flag.String("thirdQueryName", "", "Specify the name of the third sequence to count for the user-specified pattern. May be used in both or either modes. Set to the reference sequence by default")
	var s Settings

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	s = Settings{
		InFile:          flag.Arg(1),
		OutFile:         flag.Arg(2),
		QueryName:       flag.Arg(0),
		Both:            *both,
		Either:          *either,
		SecondQueryName: *secondQueryName,
		ThirdQueryName:  *thirdQueryName,
	}

	multiFaCount(s)

	// Caveat: if QueryName or SecondQueryName is misspelled, will be reference sequence by default
	//TODO: better way to do both and any modes. Maybe take queryNames as a comma-delimited slice that gets parsed. Maybe have mode option for comma-delimited list of AND,OR,etc. instead.
	// Note: Seems like "any" is a keyword so don't use that word.
	//TODO: Also, both and either modes shouldn't have to be mutually exclusive if >2 queryNames
	//TODO: exapnd either mode to not just either mode
}

// Command Group: "FASTQ TOOLS"

// Finds and replaces specified bits of sequence field

package main

import (
	"flag"
	"fmt"
	"log"
	"strings"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/fileio"
)

type FindReplaceSeq struct {
	find    []dna.Base
	replace []dna.Base
}

// readFindReplaceDNA converts a delim separated findReplace file to a []FindReplaceSeqs
func readFindReplaceDNA(filename string, delim string) []FindReplaceSeq {
	var in *fileio.EasyReader
	var findReplacePairs []FindReplaceSeq
	var currPair FindReplaceSeq
	var seqs []string
	var done bool
	var line string
	var err error

	in = fileio.EasyOpen(filename)
	for line, done = fileio.EasyNextLine(in); !done; line, done = fileio.EasyNextLine(in) {
		seqs = strings.Split(line, delim)

		if len(seqs) != 2 {
			log.Fatalf("Error: the following line:\n\"%s\"\ndoes not give two substrings when split with \"%s\"", line, delim)
		} else if len(seqs[0]) != len(seqs[1]) {
			log.Fatalf("Error: find sequence must be same length as replace sequence.\n")
		}

		currPair.find = dna.StringToBases(seqs[0])
		currPair.replace = dna.StringToBases(seqs[1])
		findReplacePairs = append(findReplacePairs, currPair)
	}
	err = in.Close()
	exception.PanicOnErr(err)

	return findReplacePairs
}

// findSeq loops through a query sequence starting at the beginning and compares it to a find sequence, base by base
// outputs true if find sequence is found in query sequence
func findSeq(seq []dna.Base, find []dna.Base, ignoreCase bool) bool {
	if len(find) > len(seq) {
		log.Fatalf("Error: Length of find sequence must be less then or equal to length of query sequence.")
	}
	for currBase := range find {
		if dna.CompareBases(seq[currBase], find[currBase], ignoreCase) != 0 {
			return false
		}
	}
	return true
}

// replacePrefix replaces bases in a query sequence with a replace sequence, starting from the beginning of the query sequence
func replacePrefix(seq []dna.Base, replace []dna.Base) []dna.Base {
	copy(seq, replace[0:])
	return seq
}

type Settings struct {
	inFile              string
	outFile             string
	findReplaceFile     string
	findReplaceDelim    string
	ignoreCase          bool
	replacedRecordsOnly bool
}

func usage() {
	fmt.Print(
		"fastqReplaceSeq - finds a sequence in a fastq file (starting from the beginning of the line) and replaces it with a new sequence of the same length. \n" +
			"Usage:\n" +
			"fastqReplaceSeq input.fastq findReplaceFile output.fastq\n" +
			"options:\n")
	flag.PrintDefaults()
}

func fastqReplaceSeq(s Settings) {
	var findReplacePairs []FindReplaceSeq
	var currPair FindReplaceSeq
	var find []dna.Base
	var replace []dna.Base
	var replaceCounter int = 0
	fq := fastq.GoReadToChan(s.inFile)
	out := fileio.EasyCreate(s.outFile)
	var err error

	findReplacePairs = readFindReplaceDNA(s.findReplaceFile, s.findReplaceDelim)

	for currRecord := range fq {
		found := false
		for _, currPair = range findReplacePairs {
			find = currPair.find
			replace = currPair.replace
			found = findSeq(currRecord.Seq, find, s.ignoreCase)
			if found {
				currRecord.Seq = replacePrefix(currRecord.Seq, replace)
				replaceCounter++
				fastq.WriteToFileHandle(out, currRecord)
				break
			}
		}
		if !found && !s.replacedRecordsOnly {
			fastq.WriteToFileHandle(out, currRecord)
		}
	}
	if replaceCounter == 0 {
		log.Fatalf("Error: No pattern(s) found in input file.\n")
	}
	err = out.Close()
	exception.PanicOnErr(err)
}

func main() {
	var expectedNumArgs int = 3
	var replaceDelim *string = flag.String("replaceDelim", "\t", "Lines in the findReplaceFile will be broken into columns based on this substring.  This should result in two substrings for each line where the first is the find and the second is the replace.")
	var ignoreCase *bool = flag.Bool("ignoreCase", true, "Ignore case when finding matches. Default is true.")
	var replacedRecordsOnly *bool = flag.Bool("replacedRecordsOnly", false, "Only output fastq records with a replacement in the sequence field. Default is false (outputs all records, regardless of whether a replacement was made).")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)
	findReplaceFile := flag.Arg(1)
	outFile := flag.Arg(2)

	s := Settings{
		inFile,
		outFile,
		findReplaceFile,
		*replaceDelim,
		*ignoreCase,
		*replacedRecordsOnly}

	fastqReplaceSeq(s)
}

package fastqReplaceSeq

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

func readFindReplaceDNA(filename string, delim string) map[dna.Base]dna.Base {
	// TODO: convert to DNA base from string
	var in *fileio.EasyReader
	var findReplaceMap = make(map[dna.Base]dna.Base)
	var words []string
	var found, done bool
	var line string
	var err error

	in = fileio.EasyOpen(filename)
	for line, done = fileio.EasyNextLine(in); !done; line, done = fileio.EasyNextLine(in) {
		words = strings.Split(line, delim)
		if len(words) != 2 {
			log.Fatalf("Error: the following line:\n\"%s\"\ndoes not give two substrings when split with \"%s\"", line, delim)
		}
		_, found = findReplaceMap[words[0]]
		if found {
			log.Fatalf("Error: this key:\"%s\" is found more than once in the findReplaceFile.\n", words[0])
		}
		findReplaceMap[words[0]] = words[1]
	}
	err = in.Close()
	exception.PanicOnErr(err)

	return findReplaceMap
}

// compareBase returns an integer related to the lexographical order of nucleotides.
// i.e. A < C < a < c < Dot < Gap.
func compareBase(alpha dna.Base, beta dna.Base, ignoreCase bool) int {
	if ignoreCase {
		alpha = dna.ToUpper(alpha)
		beta = dna.ToUpper(beta)
	}
	if alpha < beta {
		return -1
	} else if alpha > beta {
		return 1
	} else {
		return 0
	}
}

// findSeq loops through a query sequence starting at the beginning and compares it to a find sequence, base by base
// outputs true if find sequence is found in query sequence
func findSeq(seq []dna.Base, find []dna.Base, ignoreCase bool) bool {
	if len(find) > len(seq) {
		log.Fatalf("Error: Length of find sequence must be less then or equal to length of query sequence.")
	}
	found := true
	for currBase := range find {
		if found {
			if compareBase(seq[currBase], find[currBase], ignoreCase) != 0 {
				found = false
			}
		}
	}
	return found
}

// replacePrefix replaces bases in a query sequence with a replace sequence, starting from the beginning of the query sequence
func replacePrefix(seq []dna.Base, replace []dna.Base) []dna.Base {
	copy(seq, replace[0:])
	return seq
}

type Settings struct {
	inFile                   string
	outFile                  string
	findReplaceFile          string
	findReplaceDelim         string
	ignoreCase               bool
	onlyPrintReplacedRecords bool
}

func usage() {
	fmt.Print(
		"bedMath - Performs comparative arithmetic operations on float values in bed files.\n" +
			"Usage:\n" +
			"bedMath a.bed operation b.bed out.bed\n" +
			"operation may be one of the following options:" +
			"plus\tminus\ttimes\tdivideBy\n" +
			"Input bed files must be pre-sorted by coordinate.\n" +
			"options:\n")
	flag.PrintDefaults()
}

func fastqReplaceSeq(s Settings) {
	var findReplaceMap map[string]string
	fq := fastq.GoReadToChan(s.inFile)
	out := fileio.EasyCreate(s.outFile)

	findReplaceMap = readFindReplaceDNA(s.findReplaceFile, s.findReplaceDelim)

	for currRecord := range fq {
		found := false
		for find, replace = range findReplaceMap {
			found = findSeq(currRecord.Seq, find, s.ignoreCase)
			if found {
				currRecord.Seq = replacePrefix(currRecord.Seq, replace)
				fastq.WriteToFileHandle(out, currRecord)
			}
		}
		if !found && !s.onlyPrintReplacedRecords {
			fastq.WriteToFileHandle(out, currRecord)
		}
	}
}

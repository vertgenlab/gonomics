package fastqReplaceSeq

import(
	"flag"
	"fmt"
	"log"
	"strings"

	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
)

type Settings struct {
	inFile	string
	outFile	string
	findReplaceFile string
	findReplaceDelim string
}

func readFindReplacePairs(filename string, delim string) map[string]string {
	var in *fileio.EasyReader
	var findReplaceMap = make(map[string]string)
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

func replacePrefix(seq []dna.Base, find []dna.Base, replace []dna.Base, ignoreCase bool) []dna.Base {
	if ignoreCase{
		if dna.CompareSeqsIgnoreCase(seq[1:len(find)], find) == 0{
			copy(seq, replace[0:])
		}
	}else{
		if dna.CompareSeqsCaseSensitive(seq[1:len(find)], find) == 0{
			copy(seq, replace[0:])
		}
	}

	return seq
}

func fastqReplaceSeq(s Settings) {
	var findReplaceMap map[string]string
	fq := fastq.GoReadToChan(s.inFile)
	out := fileio.EasyCreate(s.outFile)

	findReplaceMap = readFindReplacePairs(s.findReplaceFile, s.findReplaceDelim)

	for i := range fq {
		if s.RetainNamesList != "" {
			if _, NameInMap = namesMap[i.Name]; !NameInMap {
				continue
			}
		}
		if s.DiscardNamesList != "" {
			if _, NameInMap = namesMap[i.Name]; NameInMap {
				continue
			}
		}
		fastq.WriteToFileHandle(out, i)
}
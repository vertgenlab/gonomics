package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"strings"
)

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

func findReplaceAnywhere(line string, findReplaceMap map[string]string) string {
	var find, replace string

	for find, replace = range findReplaceMap {
		line = strings.ReplaceAll(line, find, replace)
	}
	return line
}

func findReplaceAnyColumn(line string, delim string, findReplaceMap map[string]string) string {
	var words []string
	var i int
	var val string
	var found bool

	words = strings.Split(line, delim)
	for i = range words {
		val, found = findReplaceMap[words[i]]
		if found { // if this string is in the map
			words[i] = val //replace it's value
		}
	}
	line = strings.Join(words, "\t")
	return line
}

func findReplaceColumn(line string, delim string, columnIndex int, findReplaceMap map[string]string) string {
	var words []string
	var val string
	var found bool

	words = strings.Split(line, delim)
	val, found = findReplaceMap[words[columnIndex]]
	if found { // if this string is in the map
		words[columnIndex] = val //replace it's value
		line = strings.Join(words, "\t")
	}
	return line
}

func findAndReplace(inFile, inFileDelim, findReplaceFile, findReplaceDelim, outFile string, columnNumber int, ignoreColumns bool) {
	var findReplaceMap map[string]string
	var line string
	var in *fileio.EasyReader
	var out *fileio.EasyWriter
	var done bool
	var err error

	findReplaceMap = readFindReplacePairs(findReplaceFile, findReplaceDelim)

	out = fileio.EasyCreate(outFile)
	in = fileio.EasyOpen(inFile)
	for line, done = fileio.EasyNextLine(in); !done; line, done = fileio.EasyNextLine(in) {
		if ignoreColumns {
			line = findReplaceAnywhere(line, findReplaceMap)
		} else if columnNumber != -1 {
			line = findReplaceColumn(line, inFileDelim, columnNumber, findReplaceMap)
		} else {
			line = findReplaceAnyColumn(line, inFileDelim, findReplaceMap)
		}
		fmt.Fprintf(out, "%s\n", line)
	}
	err = in.Close()
	exception.PanicOnErr(err)
	err = out.Close()
	exception.PanicOnErr(err)
}

func usage() {
	fmt.Print(
		"findAndReplace - finds values in a file and replaces them, processes the input as a string.\n" +
			"Usage:\n" +
			"  findAndReplace inFile.txt findReplaceFile.txt outFile.txt\n" +
			"\n" +
			"findReplace file: column one is what to find, column 2 is what to replace it with\n" +
			"\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 3
	var replaceDelim *string = flag.String("replaceDelim", "\t", "Lines in the findReplaceFile will be broken into columns based on this substring.  This should result in two substrings for each line where the first is the find and the second is the replace.")
	var inDelim *string = flag.String("inDelim", "\t", "Lines in the input file should be broken into columns based on this substring and substrings in the findReplaceFile must match the entire column")
	var columnNumber *int = flag.Int("columnNumber", -1, "The index of the column to use when finding and replacing. This is zero-based (ie for column 1 use 0. Default is to search all columns.  This must be used with inDelim to know where the columns are located.")
	var ignoreColumns = flag.Bool("ignoreColumns", false, "Ignore parsing the input file into columns and replace the found substrings wherever they are in the file.")

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

	findAndReplace(inFile, *inDelim, findReplaceFile, *replaceDelim, outFile, *columnNumber, *ignoreColumns)

}

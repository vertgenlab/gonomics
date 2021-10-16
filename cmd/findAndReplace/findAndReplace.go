package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"strings"
)

//TODO
// handle both tab delimited and space delimited files, currently only tsv.
// process changing internal strings and numbers in a column. ie, change all "5"s to "9"s including when value = 5 --> 9, 65 --> 69, or 256 --> 296

func findAndReplace(inFile string, findReplaceFile string, outFile string, columnNumber int) {
	inputFile := fileio.Read(inFile) //reads input file and returns each line in file as a string
	findReplace := fileio.Read(findReplaceFile)
	//var err error
	var findReplaceMap = make(map[string]string)
	var i int
	var j int
	var words []string
	var outString string
	var outFileStrings []string
	for i = range findReplace { //making the map
		words = strings.Split(findReplace[i], "\t")
		findReplaceMap[words[0]] = words[1]
	}
	if columnNumber != -1 { //if not the default of -1, then only find and replace within the specified columnNumber.
		for i = range inputFile {
			words = strings.Split(inputFile[i], "\t")
			if val, ok := findReplaceMap[words[columnNumber]]; ok { // if this string is in the map at specified column number...
				words[columnNumber] = val
			}
			outString = strings.Join(words, "\t")
			outFileStrings = append(outFileStrings, outString)
		}
		fileio.Write(outFile, outFileStrings)

	} else { // find and replace throughout the whole file, not column specific.
		for i = range inputFile {
			words = strings.Split(inputFile[i], "\t")
			for j = range words {
				if val, ok := findReplaceMap[words[j]]; ok { // if this string is in the map
					words[j] = val //replace it's value
				}
			}
			outString = strings.Join(words, "\t")
			outFileStrings = append(outFileStrings, outString)
		}
		fileio.Write(outFile, outFileStrings)
	}
}

func usage() {
	fmt.Print(
		"findAndReplace - finds values in a file and replaces them, processes the input as a string.\n" +
			"Usage:\n" +
			"findAndReplace inFile.tsv findReplaceFile.tsv outFile.tsv\n" +
			"inFile must be a tab-separated file\n" +
			"findReplace file: column one is what to find, column 2 is what to replace it with\n" +
			"\t and must be a table separated file\n" +
			"third arg: give outfile name eg outFile.tsv\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 3
	var columnNumber *int = flag.Int("columnNumber", -1, "For findAndReplace, sets which column to search through when finding and replacing. This is zero-based (ie for column 1 use 0. Default is the whole file (columnNumber = -1).")

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

	findAndReplace(inFile, findReplaceFile, outFile, *columnNumber)

}

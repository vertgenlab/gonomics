// Command Group: "General Tools"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"strings"
)

func findReplaceOneLine(line string, find []string, replace []string) string {
	var i int

	if len(find) != len(replace) {
		log.Panic("Error: the find and replace slices were not the same length\n")
	}

	for i = range find {
		line = strings.ReplaceAll(line, find[i], replace[i])
	}
	return line
}

func findReplace(inFilename, findMe, replaceWithMe, outFilename string, findReplaceAreFiles bool) {
	var findItems, replaceItems []string
	var in *fileio.EasyReader
	var out *fileio.EasyWriter
	var line string
	var done bool
	var err error

	in = fileio.EasyOpen(inFilename)
	out = fileio.EasyCreate(outFilename)

	if findReplaceAreFiles {
		findItems = fileio.Read(findMe)
		replaceItems = fileio.Read(replaceWithMe)
	} else {
		findItems = []string{findMe}
		replaceItems = []string{replaceWithMe}
	}

	for line, done = fileio.EasyNextLine(in); !done; line, done = fileio.EasyNextLine(in) {
		line = findReplaceOneLine(line, findItems, replaceItems)
		fmt.Fprintf(out, "%s\n", line)
	}

	err = in.Close()
	exception.PanicOnErr(err)
	err = out.Close()
	exception.PanicOnErr(err)
}

func usage() {
	fmt.Print(
		"findReplace - find and replace one or more substrings in a file\n" +
			"Usage:\n" +
			" findReplace input.txt findMe replaceWithMe out.txt\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 4
	var findReplaceAreFiles *bool = flag.Bool("findReplaceAreFiles", false, "The findMe and replaceWithMe strings on the command line are actually filenames.  The files should have the same number of lines, with the string on line X of findMe, being replaced with the string on line X of replaceWithMe.")

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)
	findMe := flag.Arg(1)
	replaceWithMe := flag.Arg(2)
	outFile := flag.Arg(3)

	findReplace(inFile, findMe, replaceWithMe, outFile, *findReplaceAreFiles)
}

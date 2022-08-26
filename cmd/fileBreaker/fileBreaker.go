package fileBreaker

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
)

func FileBreaker(inFile string, chunkSize string, outFiles string) {
	var lineNum int
	var chunkNum string //zero-based number of chunks in inFile, and number of outFiles that will be made
	chunk := common.StringToInt(chunkSize)
	lines := make([]string, chunk)
	file := fileio.EasyOpen(inFile)

	for line, done := fileio.EasyNextRealLine(file); !done; line, done = fileio.EasyNextRealLine(file) {
		lineNum++
		quo := lineNum / chunk
		rem := lineNum % chunk

		if rem == 0 {
			chunkNum = fileio.IntToString(quo)
			fileio.Write(chunkNum+outFiles, lines)
			lines = make([]string, chunk)
		} else {
			lines = append(lines, line)
		}
	}
}
func usage() {
	fmt.Print(
		"fileBreaker takes a file of some length and a specified chunk size to create files with the number of " +
			"lines = chunk size. The last file will be the remaining lines. OutFiles input will be the name of all " +
			"output files with a zero-based number corresponding to the chunk in which it was found in the original file. " +
			"(Ex: 5fileName.extension if the line was found in the sixth chunk and the output file name given was 'fileName')\n" +
			"Usage:\n" +
			"fileBreaker inFile chunkSize outFile\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 3

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	infile := flag.Arg(0)
	chunkSize := flag.Arg(1)
	outFile := flag.Arg(2)

	FileBreaker(infile, chunkSize, outFile)
}

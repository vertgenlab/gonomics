package fileBreaker

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"path/filepath"
)

func FileBreaker(inFile string, chunkSize string, outFiles string) {
	var lineNum int
	var chunkNum string //zero-based number of chunks in inFile, and number of outFiles that will be made
	chunkLen := common.StringToInt(chunkSize)
	lines := make([]string, chunkLen)
	path, name := filepath.Split(outFiles)
	file := fileio.EasyOpen(inFile)

	for line, done := fileio.EasyNextRealLine(file); !done; line, done = fileio.EasyNextRealLine(file) {
		lineNum++
		quo := lineNum / chunkLen
		rem := lineNum % chunkLen

		if rem == 0 {
			lines[chunkLen-1] = line
			chunkNum = fileio.IntToString(quo)
			fileio.Write(path+chunkNum+name, lines)
			lines = make([]string, chunkLen)
		} else if rem != 0 && !done {
			lines[rem-1] = line
		} else if done { //TODO: never entering this else statement, everything else seems functional
			log.Print("in done else")
			if rem == 0 {
				lines[chunkLen-1] = line
			} else {
				lines[rem-1] = line
				log.Print("writing Oranges")
			}
			chunkNum = fileio.IntToString(quo)
			fileio.Write(path+chunkNum+name, lines)
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
	var expectedNumArgs = 3

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

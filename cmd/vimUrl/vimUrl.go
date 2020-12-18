package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
)

func usage() {
	fmt.Print(
		"vimUrl - view http url links and print data stream to stdout\n  Usage:\n  ./vimUrl link.com/file.txt\noptions:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 1
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime)
	flag.Parse()
	var toFile *string = flag.String("out", "", "provide a name to redirect data stream to a `file.txt`")
	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}
	url := flag.Arg(0)
	if *toFile != "" {
		output := fileio.EasyCreate(*toFile)
		defer output.Close()
		reader := fileio.HttpReader(url)
		for i, err := fileio.ReadLine(reader); !err; i, err = fileio.ReadLine(reader) {
			i.WriteByte('\n')
			output.Write(i.Bytes())
		}
	}

	fileio.VimUrl(url)
}

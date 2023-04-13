// Command Group: "General Tools"

package main

import (
	"flag"
	"fmt"
	"log"

	"github.com/vertgenlab/gonomics/fileio"
)

func usage() {
	fmt.Print(
		"goFetchCat - view http url links and print data stream to stdout\n  Usage:\n  ./goFetchCat link.com/file.txt\noptions:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 1
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime)
	flag.Parse()
	// TODO: add more options and functions to fetching data over http
	//var toFile *string = flag.String("out", "", "provide a name to redirect data stream to a `file.txt`")
	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}
	url := flag.Arg(0)
	str := fileio.CatUrl(url)
	fmt.Print(str)
}

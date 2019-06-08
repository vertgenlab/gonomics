package main

import(
	"flag"
	"fmt"
	"log"
	"github.com/vertgenlab/gonomics/wig"
)

func readwig(inFile string, outFile string) {
	var wigList []wig.Wig
	var err error
	wigList, err = wig.Read(inFile)
	if err != nil {
		log.Fatal(err)
	}

	wig.PrintFirst(wigList)
	wig.Write(outFile, wigList)
}

func usage() {
	fmt.Print(
		"wig_test - Tests wig package functions\n" +
			"Usage:\n" +
			" wig_test input.wig output.wig\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 2
	flag.Usage = usage
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if (len(flag.Args()) != expectedNumArgs) {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)
	outFile := flag.Arg(1)

	readwig(inFile, outFile)
}

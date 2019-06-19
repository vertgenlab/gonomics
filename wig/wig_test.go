package main

import(
	"flag"
	"fmt"
	"log"
	"github.com/vertgenlab/gonomics/wig"
)

func test_wig(inFile string, outFile string) {
	var wigList []wig.Wig
	wigList, err = wig.Read(inFile)
	if err != nil {
		return err
	}

	wig.Write(outFile, wigList)
}

func usage() {
	fmt.Print(
		"wig_test - Tests wig package functions\n" +
			"Usage:\n" +
			" wig_test input.wig output.wig\n")
	flag.PrintDefaults()
}

func main{
	var expectedNumArgs int = 2
	flag.Usage = usage
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	inFile := flag.Arg(0)
	outFile := flag.Arg(1)

	test_wig(inFile, outFile)
}

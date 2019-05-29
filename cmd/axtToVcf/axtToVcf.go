package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/axt"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
)

func usage() {
	fmt.Print(
		"AXT Alignment package - b\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 1
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}
	axtFile, _ := axt.ReadIn(flag.Arg(0))
	snps := axt.CallSnpsToVcf(axtFile)
	vcf.Sort(snps)
	vcf.PrintHeader()
	vcf.PrintVcf(snps)
}

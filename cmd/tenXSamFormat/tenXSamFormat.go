package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/sam"
	"log"
	"os"
)

func usage() {
	fmt.Print(
		"Usage:\ntenXSamFormat [options] input.sam\n")
	flag.PrintDefaults()
}

func main() {
	if len(os.Args) < 2 {
		usage()
		log.Fatalf("Error command not recognized, try something like:/ntenXSamFormat -o file $sam or\nsamtools view .bam | /ntenXSamFormat -o file /dev/stdin\n")
	}
	var threads *int = flag.Int("t", 1, "``Number of goroutines to launch")
	var output *string = flag.String("o", "/dev/stdout", "``Name for final output file")
	flag.Parse()

	//	if len(flag.Args()) != expectedNumArgs {
	//		flag.Usage()
	//		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	//	}
	sam.SwitchBxTagChannel(flag.Arg(0), *output, *threads)
}

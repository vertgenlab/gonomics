package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/sam"
	"log"
)

func usage() {
	fmt.Print(
		"Usage:\ntenXSamFormat [options] input.sam\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 1
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	var threads *int = flag.Int("threads", 4, "``Number of goroutines to launch")
	var output *string = flag.String("out", "/dev/stdout", "``Name for final output file")
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}
	sam.SwitchBxTagChannel(flag.Arg(0), *output, *threads)
}

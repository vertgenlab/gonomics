package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/giraf/binaryGiraf"
	"github.com/vertgenlab/gonomics/simpleGraph"
	"log"
	"path/filepath"
	"strings"
)

func usage() {
	fmt.Print(
		"girafCompress - GIRAF <-> GIRAF.FE conversion\n" +
			"Usage:\n" +
			" girafCompress [options] file.giraf \n" +
			"options:\n")
	flag.PrintDefaults()
}

func compress(infile string) {
	if filepath.Ext(infile) == ".giraf" {
		binaryGiraf.CompressGiraf(infile, infile+".fe")
	} else {
		log.Fatalf("ERROR: %s does not have .giraf extension", infile)
	}
}

func decompress(infile string, graph string) {
	ref := simpleGraph.Read(graph)
	if filepath.Ext(infile) == ".fe" {
		binaryGiraf.DecompressGiraf(infile, strings.TrimSuffix(infile, ".fe"), ref)
	}
}

func main() {
	var expectedNumArgs int = 1
	var decomp *bool = flag.Bool("decompress", false, "GIRAF.FE --> GIRAF, requires -ref")
	var ref *string = flag.String("ref", "", "Reference graph (.gg) used for GIRAF alignment.")
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	infile := flag.Arg(0)

	if *decomp {
		if *ref == "" {
			log.Fatalln("ERROR: -ref required for decompression")
		}
		decompress(infile, *ref)
	} else {
		if *ref != "" {
			log.Println("NOTE: -ref is not required from compression, ignoring input reference")
		}
		compress(infile)
	}
}

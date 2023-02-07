package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/bed/bedpe"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/gtf"
	"log"
)

type Settings struct {
	GtfFile string
	ContactFile string
	SizesFile string
	NearestGeneBed string
}

func great3d(s Settings) {
	genes := gtf.Read(s.GtfFile)
	contacts := bedpe.Read(s.ContactFile)
	sizes := chromInfo.ReadToMap(s.SizesFile)
	tss := gtf.GenesToTssBed(genes, sizes)
	nearestGenes := bedpe.Fill3dSpace(contacts, tss, sizes)
	if s.NearestGeneBed != "" {
		bed.Write(s.NearestGeneBed, nearestGenes)
	}
}

func usage() {
	fmt.Print(
		"great3d - Perform functional enrichment analysis by associating genomic regions\n" +
			"with their nearest gene in 3d space using bedPe contact sites.\n" +
			"Usage:\n" +
			" great3d genes.gtf contacts.bedpe genome.chrom.sizes output?\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 4
	var nearestGeneBed *string = flag.String("nearestGeneBed", "", "Write a bed representing the nearest genes in 3d space to every position in the genome.")
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}
	gtfFile := flag.Arg(0)
	contactFile := flag.Arg(1)
	sizesFile := flag.Arg(2)
	//TODO other args

	s := Settings {
		GtfFile: gtfFile,
		ContactFile: contactFile,
		SizesFile: sizesFile,
		NearestGeneBed: *nearestGeneBed,
	}

	great3d(s)
}
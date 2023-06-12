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
	GeneFile           string
	ContactFile        string
	SizesFile          string
	NearestGeneBed     string
	SinglePointGeneBed bool
}

func great3d(s Settings) {
	var tss []bed.Bed
	contacts := bedpe.Read(s.ContactFile)
	sizes := chromInfo.ReadToMap(s.SizesFile)
	if s.SinglePointGeneBed {
		tss = bed.Read(s.GeneFile)
	} else {
		genes := gtf.Read(s.GeneFile)
		tss = gtf.GenesToTssBed(genes, sizes)
	}
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
	var expectedNumArgs = 4
	var nearestGeneBed *string = flag.String("nearestGeneBed", "", "Write a bed representing the nearest genes in 3d space to every position in the genome.")
	var geneBed *bool = flag.Bool("geneBed", false, "User has provided a bed in place of a gtf for te first argument that contains information about a gene TSS already spanning a single bp distance, and therefore GTF processing can be skipped.")
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}
	geneFile := flag.Arg(0)
	contactFile := flag.Arg(1)
	sizesFile := flag.Arg(2)

	s := Settings{
		GeneFile:           geneFile,
		ContactFile:        contactFile,
		SizesFile:          sizesFile,
		NearestGeneBed:     *nearestGeneBed,
		SinglePointGeneBed: *geneBed,
	}

	great3d(s)
}

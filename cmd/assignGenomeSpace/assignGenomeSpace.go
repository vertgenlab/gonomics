// Command Group: "General Tools"

// Perform functional enrichment analysis by associating genomic regions
// with their nearest gene in 3d space using bedPe contact sites.
package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/bed/bedpe"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/gtf"
	"github.com/vertgenlab/gonomics/ontology"
	"log"
)

type Settings struct {
	GeneFile       string
	ContactFile    string
	SizesFile      string
	NearestGeneBed string
	GeneBed        bool
	Output1d       string
}

func assignGenomeSpace(s Settings) {
	var tss []bed.Bed
	sizes := chromInfo.ReadToMap(s.SizesFile)
	if s.GeneBed {
		tss = bed.Read(s.GeneFile)
	} else {
		genes := gtf.Read(s.GeneFile)
		tss = gtf.GenesToTssBed(genes, sizes, true) //always want this merged
	}
	if s.Output1d != "" {
		proximityFile := ontology.FillSpaceNoHiddenValue(tss, sizes)
		bed.Write(s.Output1d, proximityFile)
	}
	if s.ContactFile != "" {
		contacts := bedpe.Read(s.ContactFile)
		nearestGenes := ontology.Fill3dSpace(contacts, tss, sizes)
		if s.NearestGeneBed != "" {
			bed.Write(s.NearestGeneBed, nearestGenes)
		}
	}
}

func usage() {
	fmt.Print(
		"assignGenomeSpace - determine the closest gene to every base in a genome using proximity or 3D contacts\n" +
			"Usage:\n" +
			"assignGenomeSpace [options] genes.gtf genome.chrom.sizes\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs = 2
	var nearestGeneBed *string = flag.String("nearestGeneBed", "", "Write a bed representing the nearest genes in 3d space to every position in the genome.")
	var geneBed *bool = flag.Bool("geneBed", false, "If set to true user has provided a bed in place of a gtf for te first argument that contains information about a gene TSS already spanning a single bp distance, and therefore GTF processing can be skipped.")
	var proximityFile *string = flag.String("proximity", "", "If given a file name the program will output a file that contains the closest gene to a bp based only on proximity as well as any other file requested. No bedpe file necessary if this is your only desired output.")
	var contactFile *string = flag.String("contactFile", "", "If desired output is a 3D contact aware nearest gene, then provide a bedpe file with this option.")
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n", expectedNumArgs, len(flag.Args()))
	}
	geneFile := flag.Arg(0)
	sizesFile := flag.Arg(1)

	s := Settings{
		GeneFile:       geneFile,
		SizesFile:      sizesFile,
		ContactFile:    *contactFile,
		NearestGeneBed: *nearestGeneBed,
		GeneBed:        *geneBed,
		Output1d:       *proximityFile,
	}

	assignGenomeSpace(s)
}

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/bed/bedpe"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/gtf"
	"github.com/vertgenlab/gonomics/ontology"
	"github.com/vertgenlab/gonomics/ontology/gaf"
	"github.com/vertgenlab/gonomics/ontology/obo"
	"log"
)

func ontologyEnrichment(inputFile string, chromSizes string, geneFile string, annotationsFile string, oboFile string, enrichmentOut string, force bool, contactFile string) {
	var queries []bed.Bed
	var sizes map[string]chromInfo.ChromInfo
	var genes map[string]*gtf.Gene
	var contacts []bedpe.BedPe
	var annotations []gaf.Gaf
	var obos map[string]*obo.Obo

	queries = bed.Read(inputFile)
	sizes = chromInfo.ReadToMap(chromSizes)
	genes = gtf.Read(geneFile)
	if contactFile != "" {
		contacts = bedpe.Read(contactFile)
	}
	annotations, _ = gaf.Read(annotationsFile)
	obos, _ = obo.Read(oboFile, force)
	ontology.ThreeDGreat(queries, sizes, genes, contacts, annotations, obos, enrichmentOut, true, true)
}

func usage() {
	fmt.Print(
		"ontologyEnrichment will assign regions in the input data to heir closest gene provided in the gene file input either " +
			"using proximity (if no contact file is given) or by the " +
			"closest gene in 3d space if a contact file is provided. With those assigned genes, and their corresponding GO terms, " +
			"provided in the gaf and obo files, " +
			"each query region will be assigned a GO term and then enrichment for each GO term assigned will be calculated." +
			"Usage:\n" +
			"ontologyEnrichment [options] inputFile.bed chromSizes geneFile.gtf annotationsFile.gaf oboFile.obo enrichmentOut.bed\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs = 6
	var force *bool = flag.Bool("force", false, "Default is set to false, this should only be set to true if the obo should be read ignoring empty fields besides the ID.")
	var contactFile *string = flag.String("contactFile", "", "If the goal is to assign query regions to their closest 3d gene, then provide a contact file in the form of a bedpe.")
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()
	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	inputFile := flag.Arg(0)
	chromSizes := flag.Arg(1)
	geneFile := flag.Arg(2)
	annotationsFile := flag.Arg(3)
	oboFile := flag.Arg(4)
	enrichmentOut := flag.Arg(5)

	ontologyEnrichment(inputFile, chromSizes, geneFile, annotationsFile, oboFile, enrichmentOut, *force, *contactFile)
}

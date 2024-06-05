package main

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"os"
	"testing"
)

var ontologyEnrichmentTests = []struct {
	Input           string
	ChromSize       string
	GeneFile        string
	GafFile         string
	OboFile         string
	OutFile         string
	Force           bool
	ContactFile     string
	GeneEnrichments bool
	TermEnrichments bool
}{
	{
		Input:           "testdata/tiny.phastCons.bed",
		ChromSize:       "testdata/mm39.chrom.sizes", //"testdata/mm39.noGap.bed",
		GeneFile:        "testdata/tiny.tss.bed",
		GafFile:         "testdata/mgi.gaf.test.gz",
		OboFile:         "testdata/go.obo",
		OutFile:         "testdata/ontologyTest.bed",
		Force:           false,
		ContactFile:     "",
		GeneEnrichments: true,
		TermEnrichments: true,
	},
}

func TestOntologyEnrichment(t *testing.T) {
	for _, s := range ontologyEnrichmentTests {
		ontologyEnrichment(s.Input, s.ChromSize, s.GeneFile, s.GafFile, s.OboFile, s.OutFile, s.Force, s.ContactFile, s.GeneEnrichments, s.TermEnrichments)
		if bed.AllAreEqual(bed.Read(s.OutFile), bed.Read("testdata/expectedOntologies.bed")) && fileio.AreEqualIgnoreOrder("testdata/ontologyTest.geneProportions.txt", "testdata/expectedGeneProportions.txt") {
			err := os.Remove(s.OutFile)
			exception.PanicOnErr(err)
			err = os.Remove("testdata/ontologyTest.geneProportions.txt")
			exception.PanicOnErr(err)
		} else {
			log.Fatal("output for bed or gene proportions and expected didn't match")
		}
		if fileio.AreEqualIgnoreOrder("testdata/ontologyTest.termEnrichment.txt", "testdata/expectedTermEnrichment.txt") && fileio.AreEqualIgnoreOrder("testdata/ontologyTest.termProportions.txt", "testdata/expectedTermProportions.txt") {
			err := os.Remove("testdata/ontologyTest.termEnrichment.txt")
			exception.PanicOnErr(err)
			err = os.Remove("testdata/ontologyTest.termProportions.txt")
			exception.PanicOnErr(err)
			err = os.Remove("testdata/ontologyTest.inputEnrichments.txt")
		} else {
			log.Fatal("output for terms and expected didn't match")
		}
	}

}

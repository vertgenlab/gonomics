package main

import (
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
		ChromSize:       "testdata/mm39.chrom.sizes",
		GeneFile:        "mm39.all.tss.bed",
		GafFile:         "testdata/mgi.gaf.test.gz",
		OboFile:         "testdata/go.obo",
		OutFile:         "testdata/ontologyTest.bed",
		Force:           false,
		ContactFile:     "",
		GeneEnrichments: false,
		TermEnrichments: false,
	},
}

func TestOntologyEnrichment(t *testing.T) {
	for _, s := range ontologyEnrichmentTests {
		ontologyEnrichment(s.Input, s.ChromSize, s.GeneFile, s.GafFile, s.OboFile, s.OutFile, s.Force, s.ContactFile, s.GeneEnrichments, s.TermEnrichments)
		if fileio.AreEqual(s.OutFile, "testdata/expectedOntologies.bed") {
			err := os.Remove(s.OutFile)
			exception.PanicOnErr(err)
		} else {
			log.Fatal("output and expected didn't match")
		}
	}

}

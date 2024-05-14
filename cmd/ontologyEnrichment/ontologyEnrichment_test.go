package main

import (
	"testing"
)

var ontologyEnrichmentTests = []struct {
	Input       string
	ChromSize   string
	GeneFile    string
	GafFile     string
	OboFile     string
	OutFile     string
	Force       bool
	ContactFile string
}{
	{
		Input:       "mm39.phastCons.CNEEs.bed",
		ChromSize:   "testdata/mm39.chrom.sizes",
		GeneFile:    "mm39.all.tss.bed",
		GafFile:     "testdata/mgi.gaf.test.gz",
		OboFile:     "testdata/go.obo",
		OutFile:     "testdata/ontologyTest.bed",
		Force:       false,
		ContactFile: "",
	},
}

func TestOntologyEnrichment(t *testing.T) {
	for _, s := range ontologyEnrichmentTests {
		ontologyEnrichment(s.Input, s.ChromSize, s.GeneFile, s.GafFile, s.OboFile, s.OutFile, s.Force, s.ContactFile)
	}

}

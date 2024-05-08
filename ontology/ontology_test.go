package ontology

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/bed/bedpe"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/gtf"
	"github.com/vertgenlab/gonomics/ontology/gaf"
	"github.com/vertgenlab/gonomics/ontology/obo"
	"log"
	"os"
	"testing"
)

var ThreeDGreatTests = []struct {
	QueryFile       string
	ChromSizesFile  string
	GeneFile        string
	ContactsFile    string
	AnnotationsFile string
	OboFile         string
	OntOutFile      string
	//	outFile         string
	Force bool
}{
	{QueryFile: "testdata/test.bed",
		ChromSizesFile:  "testdata/hg38.chrom.sizes",
		GeneFile:        "testdata/test.gtf",
		ContactsFile:    "testdata/test.bedpe",
		AnnotationsFile: "testdata/test.gaf",
		OboFile:         "testdata/go.obo",
		OntOutFile:      "testdata/3dOntologies.bed",
		Force:           false,
	},
	{QueryFile: "testdata/test.bed",
		ChromSizesFile:  "testdata/hg38.chrom.sizes",
		GeneFile:        "testdata/test.gtf",
		ContactsFile:    "",
		AnnotationsFile: "testdata/test.gaf",
		OboFile:         "testdata/go.obo",
		OntOutFile:      "testdata/OntologiesProximity.bed",
		Force:           false,
	},
}

func TestThreeDGreat(t *testing.T) {
	var queries []bed.Bed
	var sizes map[string]chromInfo.ChromInfo
	var genes map[string]*gtf.Gene
	var contacts []bedpe.BedPe
	var annotations []gaf.Gaf
	var obos map[string]*obo.Obo
	for _, v := range ThreeDGreatTests {
		queries = bed.Read(v.QueryFile)
		sizes = chromInfo.ReadToMap(v.ChromSizesFile)
		genes = gtf.Read(v.GeneFile)
		if v.ContactsFile != "" {
			contacts = bedpe.Read(v.ContactsFile)
		}
		annotations, _ = gaf.Read(v.AnnotationsFile)
		obos, _ = obo.Read(v.OboFile, v.Force)
		ThreeDGreat(queries, sizes, genes, contacts, annotations, obos, v.OntOutFile, false, false)
	}
	if bed.AllAreEqual(bed.Read("testdata/3dOntologies.bed"), bed.Read("testdata/expected.3dOntologies.bed")) {
		err := os.Remove("testdata/3dOntologies.bed")
		exception.PanicOnErr(err)
	} else {
		log.Fatal("expected file was not a match for the output #1")
	}
	if bed.AllAreEqual(bed.Read("testdata/OntologiesProximity.bed"), bed.Read("testdata/expected.OntologiesProximity.bed")) {
		err := os.Remove("testdata/OntologiesProximity.bed")
		exception.PanicOnErr(err)
	} else {
		log.Fatal("expected file was not a match for the output #2")
	}
}

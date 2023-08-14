package ontology

import (
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/bed/bedpe"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/gtf"
	"github.com/vertgenlab/gonomics/ontology/gaf"
	"github.com/vertgenlab/gonomics/ontology/obo"
	"testing"
)

var ThreeDGreatTests = []struct {
	QueryFile       string
	ChromSizesFile  string
	GeneFile        string
	ContactsFile    string
	AnnotationsFile string
	OboFile         string
	Force           bool
}{
	{QueryFile: "testdata/haqer.bed",
		ChromSizesFile:  "testdata/hg38.chrom.sizes",
		GeneFile:        "testdata/gencode.v43.annotation.gtf",
		ContactsFile:    "testdata/GSE162819_H1_MAPS_chromatin_interactions_5kb.bedpe",
		AnnotationsFile: "testdata/goa_human.gaf",
		OboFile:         "testdata/go.obo",
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
		fmt.Printf("Reading bed.\n")
		queries = bed.Read(v.QueryFile)
		fmt.Printf("Reading chromSize.\n")
		sizes = chromInfo.ReadToMap(v.ChromSizesFile)
		fmt.Printf("Reading genes.\n")
		genes = gtf.Read(v.GeneFile)
		fmt.Printf("Reading contacts.\n")
		contacts = bedpe.Read(v.ContactsFile)
		fmt.Printf("Reading annotations.\n")
		annotations, _ = gaf.Read(v.AnnotationsFile)
		fmt.Printf("Reading obos.\n")
		obos, _ = obo.Read(v.OboFile, v.Force)
		ThreeDGreat(queries, sizes, genes, contacts, annotations, obos)
	}
}

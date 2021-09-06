package sam

import (
	"github.com/vertgenlab/gonomics/dna"
	"log"
	"strings"
)

type SingleCellAlignment struct {
	Aln Sam
	Bx  []dna.Base
	Umi []dna.Base
}

func ToSingleCellAlignment(s Sam) SingleCellAlignment {
	Bx, Umi := parseBxAndUmiFromAln(s)
	return SingleCellAlignment{
		Aln: s,
		Bx:  Bx,
		Umi: Umi,
	}
}

func parseBxAndUmiFromAln(s Sam) ([]dna.Base, []dna.Base) {
	var Bx, Umi []dna.Base
	var foundBx, foundUmi bool = false, false
	nameFields := strings.Split(s.QName, "_")
	var words []string
	for i := range nameFields {
		if strings.HasPrefix(nameFields[i], "UMI:") {
			foundUmi = true
			words = strings.Split(nameFields[i], ":")
			Umi = dna.StringToBases(words[1])
		}
		if strings.HasPrefix(nameFields[i], "BX:") {
			foundBx = true
			words = strings.Split(nameFields[i], ":")
			Bx = dna.StringToBases(words[1])
		}
	}
	if !foundBx {
		log.Fatalf("Failed to parse Barcode from single-cell read name.")
	}
	if !foundUmi {
		log.Fatalf("Failed to parse Umi from single-cell read name.")
	}
	return Bx, Umi
}

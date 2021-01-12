package giraf

import (
	"fmt"
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"testing"
)

//TODO: Will write a better high throughput test when random generator is working
func TestReadAndWrite(t *testing.T) {
	var numGiraf int = 1
	girafPath := Path{TStart: 2008, Nodes: []uint32{2015, 2017, 2018}, TEnd: 2020}
	girafNotes := []Note{Note{Tag: []byte{'B', 'Z'}, Type: 'i', Value: "5"}, Note{Tag: []byte{'G', 'P'}, Type: 'Z', Value: "30,11,35,9,23"}, Note{Tag: []byte("XO"), Type: 'i', Value: "ham5"}}
	g := &Giraf{
		QName:     "goldenState",
		QStart:    2008,
		QEnd:      2020,
		PosStrand: true,
		Path:      girafPath,
		Cigar:     []cigar.ByteCigar{{RunLen: 5, Op: 'M'}},
		AlnScore:  110335,
		MapQ:      5,
		Seq:       dna.StringToBases("ATGCG"),
		Qual:      []uint8{74, 74, 74, 74, 74},
		Notes:     girafNotes}
	fmt.Printf("%s\n", GirafToString(g))

	alpha := make([]*Giraf, 0)
	for i := 0; i < numGiraf; i++ {
		alpha = append(alpha, g)
	}
	Write("testdata/giraf.tsv", alpha)
	beta := Read("testdata/giraf.tsv")
	log.Printf("len=%d, len=%d", len(alpha), len(beta))
	if AllEqual(alpha, beta) {
		for i := 0; i < len(alpha); i++ {
			fmt.Printf("%s\n", GirafToString(alpha[i]))
		}
	} else {
		log.Fatal("Error: files are not the same...\n")
	}
	fileio.EasyRemove("testdata/giraf.tsv")
}

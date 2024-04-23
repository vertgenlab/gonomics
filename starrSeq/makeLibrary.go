package starrSeq

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"strings"
)

func MakeLibrary(upHA string, testSeqs string, downHA string, captSeq string, umiLen int, outFile string) {
	var cs, constructs, oligos []fasta.Fasta
	var construct, oligo fasta.Fasta
	var umi []dna.Base
	var tmpUMI []string
	var mid int

	testFa := fasta.Read(testSeqs)
	up := fasta.Read(upHA)
	down := fasta.Read(downHA)

	if captSeq != "" {
		cs = fasta.Read(captSeq)
	}

	if umiLen != 0 {
		for i := 0; i < umiLen; i++ {
			tmpUMI = append(tmpUMI, "N")
		}
		umi = dna.StringToBases(strings.Join(tmpUMI, ""))
	}

	for i := range testFa {
		construct = fasta.Fasta{Name: testFa[i].Name, Seq: up[0].Seq}
		construct.Seq = append(construct.Seq, testFa[i].Seq...)
		construct.Seq = append(construct.Seq, umi...)
		construct.Seq = append(construct.Seq, cs[0].Seq...)
		construct.Seq = append(construct.Seq, down[0].Seq...)
		constructs = append(constructs, construct)
	}

	for i := range constructs {
		mid = len(constructs[i].Seq) / 2
		oligo = fasta.Fasta{Name: constructs[i].Name + "_fwd"}
		oligo.Seq = constructs[i].Seq[:mid+10]
		oligos = append(oligos, oligo)
		oligo = fasta.Fasta{Name: constructs[i].Name + "_rev"}
		oligo.Seq = dna.ReverseComplementAndCopy(constructs[i].Seq[mid-10:])
		oligos = append(oligos, oligo)
	}
	fasta.Write(outFile, oligos)
}

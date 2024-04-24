package starrSeq

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"strings"
)

func checkTm(fa fasta.Fasta, mid int, minTemp, maxTemp float64) (bool, float64) {
	Tm := dna.MeltingTemp(fa.Seq[mid-10 : mid+10])
	if Tm >= minTemp && Tm <= maxTemp {
		return true, Tm
	}
	return false, Tm
}

func createOligos(seq fasta.Fasta, mid int) []fasta.Fasta {
	var oligos []fasta.Fasta
	oligo := fasta.Fasta{Name: seq.Name + "_fwd"}
	oligo.Seq = seq.Seq[:mid+10]
	oligos = append(oligos, oligo)
	oligo = fasta.Fasta{Name: seq.Name + "_rev"}
	oligo.Seq = dna.ReverseComplementAndCopy(seq.Seq[mid-10:])
	oligos = append(oligos, oligo)
	return oligos
}

func optimizeMeltingTemps(constructs []fasta.Fasta, minTemp, maxTemp float64, maxOligoSize int) {
	var mid, c int
	var pass, odd, inRange bool = false, true, false
	var oligos []fasta.Fasta
	var Tm float64

	for i := range constructs {
		mid = len(constructs[i].Seq) / 2
		inRange, _ = checkTm(constructs[i], mid, minTemp, maxTemp)
		if inRange {
			oligos = append(oligos, createOligos(constructs[i], mid)...)
			continue
		}
		for !pass {
			c++
			if odd {
				mid += c
				odd = false
			} else {
				mid -= c
				odd = true
			}
			if (!odd && mid > maxOligoSize) || (odd && mid < len(constructs[i].Seq)-maxOligoSize) {
				oligos = append(oligos, createOligos(constructs[i], mid)...)
			}
			inRange, Tm = checkTm(constructs[i], mid, minTemp, maxTemp)
			if inRange {
				oligos = append(oligos, createOligos(constructs[i], mid)...)
				pass = true
				continue
			}
		}
	}
}

func MakeLibrary(upHA string, testSeqs string, downHA string, captSeq string, umiLen int, outFile string, minTemp float64, maxTemp float64, maxOligoSize int) {
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
	if minTemp != -1 || maxTemp != -1 {
		optimizeMeltingTemps(constructs, minTemp, maxTemp, maxOligoSize)
	} else {
		for i := range constructs {
			mid = len(constructs[i].Seq) / 2
			oligo = fasta.Fasta{Name: constructs[i].Name + "_fwd"}
			oligo.Seq = constructs[i].Seq[:mid+10]
			oligos = append(oligos, oligo)
			oligo = fasta.Fasta{Name: constructs[i].Name + "_rev"}
			oligo.Seq = dna.ReverseComplementAndCopy(constructs[i].Seq[mid-10:])
			oligos = append(oligos, oligo)
		}
	}
	fasta.Write(outFile, oligos)
}

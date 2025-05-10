package starrSeq

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
)

type MakeRefSettings struct {
	Upstream      string
	Downstream    string
	Constructs    string
	OutFilePrefix string
	SepByChrom    bool
	DualBx        int
}

func Mkref(s MakeRefSettings) {
	var bd, bxbd []bed.Bed
	var ref fasta.Fasta
	var refChroms []fasta.Fasta
	var idx int

	var umi []dna.Base = dna.StringToBases("NNNNNN")
	//var cs []dna.Base = dna.StringToBases("GCTTTAAGGCCGGTCCTAGCAA")

	ref.Name = "chrSS"
	c := fasta.Read(s.Constructs)
	u := fasta.Read(s.Upstream)
	d := fasta.Read(s.Downstream)
	for i := range c {
		if s.SepByChrom {
			ref.Name = c[i].Name
			ref.Seq = []dna.Base{}
			idx = 0
		}
		ref.Seq = append(ref.Seq, u[0].Seq...)
		idx += len(u[0].Seq)
		if s.DualBx > 0 {
			bxbd = append(bxbd, bed.Bed{Chrom: ref.Name, ChromStart: idx + 20, ChromEnd: idx + 20 + s.DualBx, Name: c[i].Name + "_a", FieldsInitialized: 4})
		}
		ref.Seq = append(ref.Seq, c[i].Seq...)
		bd = append(bd, bed.Bed{Chrom: ref.Name, ChromStart: idx, ChromEnd: idx + len(c[i].Seq), Name: c[i].Name, FieldsInitialized: 4})
		idx += len(c[i].Seq)
		if s.DualBx > 0 {
			bxbd = append(bxbd, bed.Bed{Chrom: ref.Name, ChromStart: idx - 18 - s.DualBx, ChromEnd: idx - 18, Name: c[i].Name + "_b", FieldsInitialized: 4})
		}
		ref.Seq = append(ref.Seq, umi...)
		//ref.Seq = append(ref.Seq, cs...)
		ref.Seq = append(ref.Seq, d[0].Seq...)
		idx += len(umi) + len(d[0].Seq)

		if s.SepByChrom {
			refChroms = append(refChroms, ref)
		}
	}
	if !s.SepByChrom {
		refChroms = append(refChroms, ref)
	} else {
		fasta.Write(s.OutFilePrefix+".fa", refChroms)
	}
	if s.DualBx > 0 {
		bed.Write(s.OutFilePrefix+".dualBx.bed", bxbd)
	}
	bed.Write(s.OutFilePrefix+".bed", bd)
}

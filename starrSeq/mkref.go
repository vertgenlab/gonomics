package starrSeq

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
)

type MakeRefSettings struct {
	Upstream      string
	Downstream    string
	Constructs    string
	OutFilePrefix string
}

func Mkref(s MakeRefSettings) {
	var bd []bed.Bed
	var ref fasta.Fasta
	var idx int

	var umi []dna.Base = dna.StringToBases("NNNNN")
	//var cs []dna.Base = dna.StringToBases("GCTTTAAGGCCGGTCCTAGCAA")

	o := fileio.EasyCreate(s.OutFilePrefix + ".fa")
	ref.Name = "chrSS"
	c := fasta.Read(s.Constructs)
	u := fasta.Read(s.Upstream)
	d := fasta.Read(s.Downstream)
	for i := range c {
		ref.Seq = append(ref.Seq, u[0].Seq...)
		idx += len(u[0].Seq)
		ref.Seq = append(ref.Seq, c[i].Seq...)
		bd = append(bd, bed.Bed{Chrom: "chrSS", ChromStart: idx, ChromEnd: idx + len(c[i].Seq), Name: c[i].Name, FieldsInitialized: 4})
		ref.Seq = append(ref.Seq, umi...)
		//ref.Seq = append(ref.Seq, cs...)
		ref.Seq = append(ref.Seq, d[0].Seq...)
		idx += len(c[i].Seq) + len(d[0].Seq)
	}
	fasta.WriteFasta(o, ref, 50)
	err := o.Close()
	exception.PanicOnErr(err)
	bed.Write(s.OutFilePrefix+".bed", bd)
}

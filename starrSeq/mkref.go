package starrSeq

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
)

func Mkref(upstream, constructs, downstream, outputPrefix string) {
	var bd []bed.Bed
	var ref fasta.Fasta
	var idx int
	o := fileio.EasyCreate(outputPrefix + ".fa")
	ref.Name = "chrSS"
	c := fasta.Read(constructs)
	u := fasta.Read(upstream)
	d := fasta.Read(downstream)
	for i := range c {
		ref.Seq = append(ref.Seq, u[0].Seq...)
		idx += len(u[0].Seq)
		ref.Seq = append(ref.Seq, c[i].Seq...)
		bd = append(bd, bed.Bed{Chrom: "chrSS", ChromStart: idx, ChromEnd: idx + len(c[i].Seq), Name: c[i].Name, FieldsInitialized: 4})
		ref.Seq = append(ref.Seq, d[0].Seq...)
		idx += len(c[i].Seq) + len(d[0].Seq)
	}
	fasta.WriteFasta(o, ref, 50)
	err := o.Close()
	exception.PanicOnErr(err)
	bed.Write(outputPrefix+".bed", bd)
}

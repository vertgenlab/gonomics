package fastq
//TODO work in progress to deal with sequences that contain Ns
/*
import (
	"github.com/vertgenlab/gonomics/dna"
)

func TrimNs(fq *Fastq) {
	var i, j int
	//trim front
	for i = 0; i < len(fq.Seq); {
		if fq.Seq[i] == dna.N {
			i++
		} else {
			fq.Seq = fq.Seq[i:]
			fq.Qual = fq.Qual[i:]
		}
	}
	for j = len(fq.Seq) - 1; j >= 0; {
		if fq.Seq[j] == dna.N {
			j--
		} else {
			fq.Seq = fq.Seq[:j]
			fq.Qual = fq.Qual[:j]
		}
	}
}*/

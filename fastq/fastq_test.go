package fastq

import (
	"fmt"
	"github.com/vertgenlab/gonomics/fasta"
	"testing"
)

var readWriteTests = []struct {
	filename string // input
}{
	{"testdata/test.fastq"},
	//{"testdata/CL12-3_16w_test.fastq"},
}

func TestRead(t *testing.T) {
	for _, test := range readWriteTests {
		_ = Read(test.filename)
		fastq := Read(test.filename)

		//testing a few functions implemented
		//fmt.Println("Value: ", string(fastq[0].Qual[0]), "ASCII: ", fastq[0].Qual[0])
		//fmt.Println(ErrorRate(fastq[0].Qual))
		//fmt.Println(fastq[0].Seq[0], string(fastq[0].Qual[0]))
		//fmt.Println(FromBase(fastq[0].Seq, PhredToPError(fastq[0].Qual)))
		//PrintFastq(fastq)
		//rev := ReverseComplementAll(fastq)
		//PrintFastq(rev)
		//qq := FromDna(fastq[0].Seq, ErrorRate(fastq[0].Qual))
		//for i :=0; i < len(qq);i++ {
		//	fmt.Println(qq[i])
		//}
		fastaFile := fasta.Read("testdata/multiCHr_test.fa")
		samTest := GSW(fastaFile, fastq)
		for i := 0; i < len(samTest); i++ {
			fmt.Println(samTest[i])
		}
	}

}

package mHiC

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/fileio"
)

func CutSiteTrim(unMappedFastq, outFastq string, cutSite []dna.Base, minLen int) {
	var testBases []dna.Base
	var cutPos int
	var found bool = false

	o := fileio.EasyCreate(outFastq)
	inFq := fastq.GoReadToChan(unMappedFastq)

	for fq := range inFq {
		for b := 0; b <= len(fq.Seq)-len(cutSite); b++ {
			testBases = fq.Seq[b : b+len(cutSite)]
			if dna.CompareSeqsIgnoreCase(testBases, cutSite) == 0 { //if there are multiple cutsites it will update to the second one. I think this is the behavior that we want
				cutPos = b
				found = true
			}
		}
		if found && cutPos+len(cutSite)/2 >= minLen {
			fastq.WriteToFileHandle(o, fastq.Fastq{
				Name: fq.Name,
				Seq:  fq.Seq[0 : cutPos+len(cutSite)/2],
				Qual: fq.Qual[0 : cutPos+len(cutSite)/2],
			})
		}
		found = false
	}
	err := o.Close()
	exception.PanicOnErr(err)
}

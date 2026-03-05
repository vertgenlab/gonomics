package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"os"
	"testing"

	"github.com/vertgenlab/gonomics/axt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
)

var axtStatsTests = []struct {
	inAXT      string
	outStats   string
	bedregions string
	expected   string
}{
	{inAXT: "testdata/test.axt", outStats: "testdata/out.stats.txt", bedregions: "", expected: "testdata/exp.stats.txt"},
	{inAXT: "testdata/test.axt", outStats: "testdata/out.stats.bed.txt", bedregions: "testdata/bed1.bed", expected: "testdata/exp.stats.bed.txt"},
}

func TestAxtToFasta(t *testing.T) {
	log.SetFlags(0)
	var testAxt = axt.Axt{
		RName:      "chr19",
		RStart:     11,
		REnd:       40,
		QName:      "chr11",
		QStart:     2,
		QEnd:       31,
		QStrandPos: false,
		Score:      3500,
		RSeq:       dna.StringToBases("TCTGNNNNTAAACCACCNNNNATGACAAGC"),
		QSeq:       dna.StringToBases("TCAGCTCATAAATCACCTCCCATGACAAGC"),
	}
	var testFa = &fasta.Fasta{
		Name: "chr19",
		Seq:  dna.StringToBases("AAAAAAAAAATCTGNNNNTAAACCACCNNNNATGACAAGC"),
	}

	seq := axtSeq(testAxt, testFa.Seq)
	log.Printf("%s %s\n", seq.Name, dna.BasesToString(seq.Seq))

	for _, v := range axtStatsTests {
		axtStats(v.inAXT, v.outStats, &v.bedregions)
		if !fileio.AreEqual(v.outStats, v.expected) {
			t.Errorf("ERROR in axTools -stats: The output (%s) and expected (%s) files are not equal", v.outStats, v.expected)
		} else {
			exception.PanicOnErr(os.Remove(v.outStats))
		}
	}

}

package convert

import (
	"os"
	"strings"
	"testing"

	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/sam"

	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/vcf"
)

var seqA []dna.Base = dna.StringToBases("--TTTC--ATGAATAATCA")
var seqB []dna.Base = dna.StringToBases("CCATTCCAA--CAGAATNA")
var inputFa []fasta.Fasta = []fasta.Fasta{{"eggplant", seqA}, {"squash", seqB}}
var var1 vcf.Vcf = vcf.Vcf{Chr: "chr1", Pos: 1, Id: ".", Ref: "T", Alt: strings.Split("A", ","), Qual: 100.0, Filter: "PASS", Info: ".", Format: []string{"."}}
var var2 vcf.Vcf = vcf.Vcf{Chr: "chr1", Pos: 4, Id: ".", Ref: "C", Alt: strings.Split("CCA", ","), Qual: 100.0, Filter: "PASS", Info: ".", Format: []string{"."}}
var var3 vcf.Vcf = vcf.Vcf{Chr: "chr1", Pos: 5, Id: ".", Ref: "ATG", Alt: strings.Split("A", ","), Qual: 100.0, Filter: "PASS", Info: ".", Format: []string{"."}}
var var4 vcf.Vcf = vcf.Vcf{Chr: "chr1", Pos: 8, Id: ".", Ref: "A", Alt: strings.Split("C", ","), Qual: 100.0, Filter: "PASS", Info: ".", Format: []string{"."}}
var var5 vcf.Vcf = vcf.Vcf{Chr: "chr1", Pos: 10, Id: ".", Ref: "T", Alt: strings.Split("G", ","), Qual: 100.0, Filter: "PASS", Info: ".", Format: []string{"."}}
var var6 vcf.Vcf = vcf.Vcf{Chr: "chr1", Pos: 14, Id: ".", Ref: "C", Alt: strings.Split("N", ","), Qual: 100.0, Filter: "PASS", Info: ".", Format: []string{"."}}
var expected []vcf.Vcf = []vcf.Vcf{var1, var2, var3, var4, var5}
var expectedSubOnly []vcf.Vcf = []vcf.Vcf{var1, var4, var5}
var expectedRetainN []vcf.Vcf = []vcf.Vcf{var1, var2, var3, var4, var5, var6}

func TestSamToBedWithDeletions(t *testing.T) {
	var outBeds []bed.Bed
	var j bed.Bed
	var err error
	inSam, _ := sam.Read("testdata/test1.sam")
	out := fileio.EasyCreate("testdata/out.samToBedDel.bed")
	for _, i := range inSam {
		outBeds = SamToBedWithDeletions(i)
		for _, j = range outBeds {
			bed.WriteToFileHandle(out, j)
		}
	}
	err = out.Close()
	exception.PanicOnErr(err)
	if !fileio.AreEqual("testdata/out.samToBedDel.bed", "testdata/exp.samToBedDel.bed") {
		t.Errorf("ERROR in SamToBedWithDeletions: expected != actual.")
	} else {
		err = os.Remove("testdata/out.samToBedDel.bed")
		exception.PanicOnErr(err)
	}
}

func TestSamToBed(t *testing.T) {
	//TODO
}

func TestPairwiseFaToVcf(t *testing.T) { //this test is for the default settings.
	var err error
	out := fileio.EasyCreate("tmp.txt")
	PairwiseFaToVcf(inputFa, "chr1", out, false, false)
	err = out.Close()
	exception.PanicOnErr(err)
	input, _ := vcf.Read("tmp.txt")
	if !vcf.AllEqual(input, expected) {
		t.Errorf("Pairwise VCF results do not match.")
	}
	err = os.Remove("tmp.txt")
	if err != nil {
		exception.PanicOnErr(err)
	}
}

var ThreeWayFaToVcfTests = []struct {
	InFile       string
	OutFile      string
	ExpectedFile string
	Chrom        string
}{
	{InFile: "testdata/threeWayFaToVcf/input.fa", OutFile: "tmp.threeWayFaToVcf.txt", ExpectedFile: "testdata/threeWayFaToVcf/expected.vcf", Chrom: "chr1"},
}

func TestThreeWayFaToVcf(t *testing.T) {
	var err error
	var out *fileio.EasyWriter
	for _, v := range ThreeWayFaToVcfTests {
		threeWayInputFa := fasta.Read(v.InFile)
		out = fileio.EasyCreate(v.OutFile)
		ThreeWayFaToVcf(threeWayInputFa, v.Chrom, out)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Threeway VCF results do not match.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}

func TestPairwiseFaToVcfRetainN(t *testing.T) {
	var err error
	out := fileio.EasyCreate("tmpRetainN.txt")
	PairwiseFaToVcf(inputFa, "chr1", out, false, true)
	err = out.Close()
	exception.PanicOnErr(err)
	input, _ := vcf.Read("tmpRetainN.txt")
	if !vcf.AllEqual(input, expectedRetainN) {
		t.Errorf("Pairwise VCF results do not match in retainN test.")
	}
	err = os.Remove("tmpRetainN.txt")
	if err != nil {
		exception.PanicOnErr(err)
	}
}

func TestPairwiseFaToVcfSubstitutionsOnly(t *testing.T) {
	var err error
	out := fileio.EasyCreate("tmpSub.txt")
	PairwiseFaToVcf(inputFa, "chr1", out, true, false)
	err = out.Close()
	exception.PanicOnErr(err)
	input, _ := vcf.Read("tmpSub.txt")
	if !vcf.AllEqual(input, expectedSubOnly) {
		t.Errorf("Pairwise VCF results do not match in subsitutionsOnly test.")
	}
	err = os.Remove("tmpSub.txt")
	if err != nil {
		exception.PanicOnErr(err)
	}
}

var BedValuesToWigTests = []struct {
	InFile          string
	ChromSizeFile   string
	Missing         float64
	OutFile         string
	ExpectedFile    string
	Method          string
	UseRange        bool
	AnnotationField int
}{
	{InFile: "testdata/test.bed",
		ChromSizeFile:   "testdata/ref.chrom.sizes",
		Missing:         0,
		OutFile:         "testdata/name.tmp.wig",
		ExpectedFile:    "testdata/name.Expected.wig",
		Method:          "Name",
		UseRange:        false,
		AnnotationField: 0},
	{InFile: "testdata/test.bed",
		ChromSizeFile:   "testdata/ref.chrom.sizes",
		Missing:         0,
		OutFile:         "testdata/score.tmp.wig",
		ExpectedFile:    "testdata/score.Expected.wig",
		Method:          "Score",
		UseRange:        false,
		AnnotationField: 0},
	{InFile: "testdata/test.bed",
		ChromSizeFile:   "testdata/ref.chrom.sizes",
		Missing:         -1,
		OutFile:         "testdata/name.missing.tmp.wig",
		ExpectedFile:    "testdata/name.missing.Expected.wig",
		Method:          "Name",
		UseRange:        false,
		AnnotationField: 0},
	{InFile: "testdata/test.annotation.bed",
		ChromSizeFile:   "testdata/annotation.chrom.sizes",
		Missing:         -1,
		OutFile:         "testdata/tmp.annotation.wig",
		ExpectedFile:    "testdata/expected.annotation.wig",
		Method:          "Annotation",
		UseRange:        false,
		AnnotationField: 0},
}

/* TODO: Fix BedToWig
func TestBedValuesToWig(t *testing.T) {
	var err error
	var reference map[string]chromInfo.ChromInfo
	var wigs []wig.Wig
	for _, v := range BedValuesToWigTests {
		reference = chromInfo.ReadToMap(v.ChromSizeFile)
		wigs = BedValuesToWig(v.InFile, reference, v.Missing, v.Method, v.UseRange, v.AnnotationField)
		wig.Write(v.OutFile, wigs)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in BedValuesToWig. Output was not as expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}

var BedGraphToWigTests = []struct {
	InFile        string
	ChromSizeFile string
	Missing       float64
	OutFile       string
	ExpectedFile  string
}{
	{"testdata/test.bedGraph", "testdata/ref.chrom.sizes", -10, "testdata/bedGraphToWig.tmp.wig", "testdata/bedGraphToWig.expected.wig"},
}

func TestBedGraphToWig(t *testing.T) {
	var err error
	var reference map[string]chromInfo.ChromInfo
	var wigs []wig.Wig
	for _, v := range BedGraphToWigTests {
		reference = chromInfo.ReadToMap(v.ChromSizeFile)
		wigs = BedGraphToWig(v.InFile, reference, v.Missing)
		wig.Write(v.OutFile, wigs)
		if !fileio.AreEqual(v.OutFile, v.ExpectedFile) {
			t.Errorf("Error in BedGraphToWig. Output was not as expected.")
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
*/

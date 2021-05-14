package vcf

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/variant"
	"strings"
	"testing"
)

var testSub = Vcf{
	Chr: "subTest",
	Pos: 1,
	Ref: "A",
	Alt: []string{"T", "G", "C"},
}

var expSub0 = variant.Substitution{
	Chr: "subTest", Pos: 0, Ref: dna.A, Alt: dna.T,
}

var expSub1 = variant.Substitution{
	Chr: "subTest", Pos: 0, Ref: dna.A, Alt: dna.G,
}

var expSub2 = variant.Substitution{
	Chr: "subTest", Pos: 0, Ref: dna.A, Alt: dna.C,
}

var testIns = Vcf{
	Chr: "insTest",
	Pos: 1,
	Ref: "A",
	Alt: []string{"ATG"},
}

var expIns = variant.Insertion{
	Chr: "insTest", Pos: 1, Seq: []dna.Base{dna.T, dna.G},
}

var testDel = Vcf{
	Chr: "delTest",
	Pos: 1,
	Ref: "ATG",
	Alt: []string{"A"},
}

var expDel = variant.Deletion{
	Chr: "delTest", Start: 1, End: 3,
}

var testDelins = Vcf{
	Chr: "delinsTest",
	Pos: 1,
	Ref: "ATG",
	Alt: []string{"AC", "A", "ATGC"},
}

var expDelins0 = variant.Delins{
	Chr: "delinsTest", Start: 1, End: 3, InsSeq: []dna.Base{dna.C},
}

var expDelins1 = variant.Deletion{
	Chr: "delinsTest", Start: 1, End: 3,
}

var expDelins2 = variant.Insertion{
	Chr: "delinsTest", Pos: 3, Seq: []dna.Base{dna.C},
}

func sendTestVcf() <-chan Vcf {
	c := make(chan Vcf, 10)
	c <- testSub
	c <- testIns
	c <- testDel
	c <- testDelins
	close(c)
	return c
}

func TestGoChanToVariants(t *testing.T) {
	inputVcfChan := sendTestVcf()
	answer := GoChanToVariants(inputVcfChan)

	var actualSub variant.Substitution
	var actualIns variant.Insertion
	var actualDel variant.Deletion
	var actualDelins variant.Delins
	var recdVcf Vcf

	actualSub = <-answer.Substitutions
	recdVcf = <-answer.Records
	if actualSub != expSub0 || !isEqual(testSub, recdVcf) {
		t.Errorf("problem with expSub0")
	}
	actualSub = <-answer.Substitutions
	recdVcf = <-answer.Records
	if actualSub != expSub1 || !isEqual(testSub, recdVcf) {
		t.Errorf("problem with expSub1")
	}

	actualSub = <-answer.Substitutions
	recdVcf = <-answer.Records
	if actualSub != expSub2 || !isEqual(testSub, recdVcf) {
		t.Errorf("problem with expSub2")
	}

	actualIns = <-answer.Insertions
	recdVcf = <-answer.Records
	if !equalIns(actualIns, expIns) || !isEqual(testIns, recdVcf) {
		t.Errorf("problem with expIns")
	}

	actualDel = <-answer.Deletions
	recdVcf = <-answer.Records
	if actualDel != expDel || !isEqual(testDel, recdVcf) {
		t.Errorf("problem with expDel")
	}

	actualDelins = <-answer.Delins
	recdVcf = <-answer.Records
	if !equalDelins(actualDelins, expDelins0) || !isEqual(testDelins, recdVcf) {
		t.Errorf("problem with expDelins0")
	}

	actualDel = <-answer.Deletions
	recdVcf = <-answer.Records
	if actualDel != expDelins1 || !isEqual(testDelins, recdVcf) {
		t.Errorf("problem with expDelins1")
	}

	actualIns, moreToRead := <-answer.Insertions
	recdVcf = <-answer.Records
	if !equalIns(actualIns, expDelins2) || !isEqual(testDelins, recdVcf) {
		t.Errorf("problem with expDelins2")
	}

	_, moreToRead = <-answer.Substitutions
	if moreToRead {
		t.Errorf("chan did not close")
	}

	_, moreToRead = <-answer.Insertions
	if moreToRead {
		t.Errorf("chan did not close")
	}

	_, moreToRead = <-answer.Deletions
	if moreToRead {
		t.Errorf("chan did not close")
	}

	_, moreToRead = <-answer.Delins
	if moreToRead {
		t.Errorf("chan did not close")
	}

	_, moreToRead = <-answer.Records
	if moreToRead {
		t.Errorf("chan did not close")
	}
}

func equalIns(a variant.Insertion, b variant.Insertion) bool {
	return a.Chr == b.Chr && a.Pos == b.Pos && dna.CompareSeqsCaseSensitive(a.Seq, b.Seq) == 0
}

func equalDelins(a variant.Delins, b variant.Delins) bool {
	return a.Chr == b.Chr && a.Start == b.Start && a.End == b.End && dna.CompareSeqsCaseSensitive(a.InsSeq, b.InsSeq) == 0
}

var testParse1 string = "."
var testParse2 string = "ACGT"
var testParse3 string = "AC]GT["

func BenchmarkCanParseSymbol(b *testing.B) {
	for i := 0; i < b.N; i++ {
		canParseSymbol(testParse1)
		canParseSymbol(testParse2)
		canParseSymbol(testParse3)
	}
}

func BenchmarkCanParseBase(b *testing.B) {
	for i := 0; i < b.N; i++ {
		canParseACGT(testParse1)
		canParseACGT(testParse2)
		canParseACGT(testParse3)
	}
}

func canParseSymbol(s string) bool {
	return !strings.ContainsAny(s, ":>[.")
}

func canParseACGT(s string) bool {
	for i := range s {
		if !validBase(s[i]) {
			return false
		}
	}
	return true
}

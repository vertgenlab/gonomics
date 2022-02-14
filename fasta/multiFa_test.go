package fasta

import (
	"github.com/vertgenlab/gonomics/dna"
	"testing"
)

var seqThreeA = dna.StringToBases("ACGTacgTCATCATCATTACTACTAC")
var seqThreeB = dna.StringToBases("acgtACGTACGT")
var seqThreeC = dna.StringToBases("ACGTACGTACGTT")
var rcSeqThreeA = dna.StringToBases("GTAGTAGTAATGATGATGAcgtACGT")
var rcSeqThreeB = dna.StringToBases("ACGTACGTacgt")
var rcSeqThreeC = dna.StringToBases("AACGTACGTACGT")
var allRevCompTests = []struct {
	input    []Fasta
	expected []Fasta
}{
	{[]Fasta{{"apple", seqThreeA}, {"banana", seqThreeB}, {"carrot", seqThreeC}}, []Fasta{{"apple", rcSeqThreeA}, {"banana", rcSeqThreeB}, {"carrot", rcSeqThreeC}}},
}
var segSite1 = dna.StringToBases("AATCCTATTCA")
var segSite2 = dna.StringToBases("AATCCTAATCA")
var segSite3 = dna.StringToBases("AATCCTATTCG")
var expectedSegSite1 = dna.StringToBases("TA")
var expectedSegSite2 = dna.StringToBases("AA")
var expectedSegSite3 = dna.StringToBases("TG")
var segSiteInput = []Fasta{{"apple", segSite1}, {"banana", segSite2}, {"carrot", segSite3}}
var segSiteExpected = []Fasta{{"apple", expectedSegSite1}, {"banana", expectedSegSite2}, {"carrot", expectedSegSite3}}

func TestSegregatingSites(t *testing.T) {
	input := SegregatingSites(segSiteInput)
	if !AllAreEqual(input, segSiteExpected) {
		t.Errorf("Expected segregating sites to give %v columns, but got %v.", len(segSiteExpected[0].Seq), len(input[0].Seq))
	}
}

func TestNumSegregatingSites(t *testing.T) {
	currSites := NumSegregatingSites(segSiteInput)
	if currSites != 2 {
		t.Errorf("Error in NumSegregatingSites. Found: %v. Expected: %v.", currSites, 2)
	}
}

func TestPairwiseMutationDistanceReferenceWindow(t *testing.T) {
	tmpDistance, incompleteWindow, alnEnd := PairwiseMutationDistanceReferenceWindow(Fasta{"apple", segSite1}, Fasta{"banana", segSite2}, 0, 10)
	if tmpDistance != 1 {
		t.Errorf("error in PairwiseMutationDistanceReferenceWindow. Distance: %v not as expected: %v.", tmpDistance, 1)
	}
	if incompleteWindow {
		t.Errorf("Error in PairwiseMutationDistanceReferenceWindow. Window expected to be complete.")
	}
	if alnEnd != 10 {
		t.Errorf("Error in pairwiseMutationDistanceReferenceWindow. AlnEnd: %v. Expected: %v.", alnEnd, 10)
	}
}

func TestPairwiseMutationDistanceInRange(t *testing.T) {
	tmpDistance := PairwiseMutationDistanceInRange(Fasta{"apple", segSite1}, Fasta{"banana", segSite2}, 0, 10)
	if tmpDistance != 1 {
		t.Errorf("error in PairwiseMutationDistanceReferenceWindow. Distance: %v not as expected: %v.", tmpDistance, 1)
	}
}

var DistFaInput = []Fasta{
	Fasta{"apple", dna.StringToBases("AAAAAAAAAAAAAAA")},
	Fasta{"grape", dna.StringToBases("AAaaAAAAATAAAAA")},
	Fasta{"fruit", dna.StringToBases("AAAAAAAA-AAA--A")},
}

var DistFaExpected = []Fasta{
	Fasta{"apple", dna.StringToBases("AAAAAAAAAA")},
	Fasta{"grape", dna.StringToBases("AAAAAATAAA")},
	Fasta{"fruit", dna.StringToBases("AAAAAAAAAA")},
}

func TestDistColumn(t *testing.T) {
	tmp := DistColumn(DistFaInput)
	if !AllAreEqual(DistFaExpected, tmp) {
		t.Errorf("Error in DistColumn.\nExpected:\n %v\nFound:\n%v.", DistFaExpected, tmp)
	}
}

var MissingMultInput = []Fasta{
	Fasta{"apple", dna.StringToBases("AAAAAAAAAAAAAAA")},
	Fasta{"grape", dna.StringToBases("AATCAAAAATAAAAA")},
	Fasta{"fruit", dna.StringToBases("---------------")},
}

var MissingMultExpected = []Fasta{
	Fasta{"apple", dna.StringToBases("AAAAAAAAAAAAAAA")},
	Fasta{"grape", dna.StringToBases("AATCAAAAATAAAAA")},
}

func TestRemoveMissingMult(t *testing.T) {
	tmp := RemoveMissingMult(MissingMultInput)
	if !AllAreEqual(tmp, MissingMultExpected) {
		t.Errorf("Error in RemoveMissingMult.\nExpected:\n %v\nFound:\n%v.", MissingMultExpected, tmp)
	}
}

var CopySubSetInput = []Fasta{
	Fasta{"apple", dna.StringToBases("AAAAAAAAAA")},
	Fasta{"grape", dna.StringToBases("AAAAAATAAA")},
	Fasta{"fruit", dna.StringToBases("AAAAAAAAAA")},
}

var CopySubSetExpected = []Fasta{
	Fasta{"apple", dna.StringToBases("AAAAAA")},
	Fasta{"grape", dna.StringToBases("AAAATA")},
	Fasta{"fruit", dna.StringToBases("AAAAAA")},
}

func TestCopySubset(t *testing.T) {
	tmp := CopySubset(CopySubSetInput, 2, 8)
	if !AllAreEqual(tmp, CopySubSetExpected) {
		t.Errorf("Error in CopySubSet.\nExpected:\n %v\nFound:\n%v.", CopySubSetExpected, tmp)
	}
}

var AlnPosToRefPosInput Fasta = Fasta{"apple", dna.StringToBases("AA---AAAAA")}

func TestAlnPosToRefPos(t *testing.T) {
	tmp := AlnPosToRefPos(AlnPosToRefPosInput, 6)
	if tmp != 3 {
		t.Errorf("Error in AlnPosToRefPos. Expected: %v. Found: %v.", 3, tmp)
	}
}

var RefPosToAlnPosInput Fasta = Fasta{"apple", dna.StringToBases("AA---AAAAA")}

func TestRefPosToAlnPos(t *testing.T) {
	tmp := RefPosToAlnPos(RefPosToAlnPosInput, 5)
	if tmp != 8 {
		t.Errorf("Error in refPosToAlnPos. Expected: %v. Found: %v.", 8, tmp)
	}
}


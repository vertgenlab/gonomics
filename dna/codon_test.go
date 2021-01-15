package dna

import (
	"testing"
)

var expected []*Codon = []*Codon{{A, A, G}, {T, G, A}, {A, A, C}}
var testSeq []Base = StringToBases("AAGTGAAAC")

var test2DNA []Base = StringToBases("atgcagatttttgtgaaaaccctgaccggcaaaacctga")
var expectedShortProtein string = "MQIFVKTLTGKT*"
var expectedProtein string = "MetGlnIlePheValLysThrLeuThrGlyLysThrTer"

func TestBasesToCodon(t *testing.T) {
	input, err := BasesToCodons(testSeq)
	for i := range input {
		if *input[i] != *expected[i] || err != nil {
			t.Errorf("Do not match. Input: %s. Expected: %s.", *input[i], *expected[i])
		}
	}
}

func TestTranslateToString(t *testing.T) {
	input, err := TranslateToString(test2DNA)
	if input != expectedProtein || err != nil {
		t.Errorf("Do not match. Input : %s. Expected: %s.", input, expectedProtein)
	}
}

func TestTranslateToShortString(t *testing.T) {
	input, err := TranslateToShortString(test2DNA)
	if input != expectedShortProtein || err != nil {
		t.Errorf("Do not match. Input : %s. Expected: %s.", input, expectedShortProtein)
	}
}

func TestAminoAcidToShortString(t *testing.T) {
	input, err := TranslateToShortString(test2DNA)
	if input != expectedShortProtein || err != nil {
		t.Errorf("Do not match. Input : %s. Expected: %s.", input, expectedShortProtein)
	}
}

func TestCodonsToSeq(t *testing.T) {
	input := CodonsToSeq(expected)
	if CompareSeqsCaseSensitive(input, testSeq) != 0 {
		t.Errorf("Do not match. Input : %s. Expected: %s.", input, testSeq)
	}
}

var gdf2mrna string =
	"ATGTGTCCTGGGGCACTGTGGGTGGCCCTGCCCCTGCTGTCCCTGCTGGCTGGCTCCCTACAGGGGAAGC" +
	"CACTGCAGAGCTGGGGACGAGGGTCTGCTGGGGGAAACGCCCACAGCCCACTGGGGGTGCCTGGAGGTGG" +
	"GCTGCCTGAGCACACCTTCAACCTGAAGATGTTTCTGGAGAACGTGAAGGTGGATTTCCTGCGCAGCCTT" +
	"AACCTGAGTGGGGTCCCTTCGCAGGACAAAACCAGGGTGGAGCCGCCGCAGTACATGATTGACCTGTACA" +
	"ACAGGTACACGTCCGATAAGTCGACTACGCCAGCGTCCAACATTGTGCGGAGCTTCAGCATGGAAGATGC" +
	"CATCTCCATAACTGCCACAGAGGACTTCCCCTTCCAGAAGCACATCTTGCTCTTCAACATCTCCATTCCT" +
	"AGGCATGAGCAGATCACCAGAGCTGAGCTCCGACTCTATGTCTCCTGTCAAAATCACGTGGACCCCTCTC" +
	"ATGACCTGAAAGGAAGCGTGGTCATTTATGATGTTCTGGATGGAACAGATGCCTGGGATAGTGCTACAGA" +
	"GACCAAGACCTTCCTGGTGTCCCAGGACATTCAGGATGAGGGCTGGGAGACCTTGGAAGTGTCCAGCGCC" +
	"GTGAAGCGCTGGGTCCGGTCCGACTCCACCAAGAGCAAAAATAAGCTGGAAGTGACTGTGGAGAGCCACA" +
	"GGAAGGGCTGCGACACGCTGGACATCAGTGTCCCCCCAGGTTCCAGAAACCTGCCCTTCTTTGTTGTCTT" +
	"CTCCAATGACCACAGCAGTGGGACCAAGGAGACCAGGCTGGAGCTGAGGGAGATGATCAGCCATGAACAA" +
	"GAGAGCGTGCTCAAGAAGCTGTCCAAGGACGGCTCCACAGAGGCAGGTGAGAGCAGTCACGAGGAGGACA" +
	"CGGATGGCCACGTGGCTGCGGGGTCGACTTTAGCCAGGCGGAAAAGGAGCGCCGGGGCTGGCAGCCACTG" +
	"TCAAAAGACCTCCCTGCGGGTAAACTTCGAGGACATCGGCTGGGACAGCTGGATCATTGCACCCAAGGAG" +
	"TATGAAGCCTACGAGTGTAAGGGCGGCTGCTTCTTCCCCTTGGCTGACGATGTGACGCCGACGAAACACG" +
	"CTATCGTGCAGACCCTGGTGCATCTCAAGTTCCCCACAAAGGTGGGCAAGGCCTGCTGTGTGCCCACCAA" +
	"ACTGAGCCCCATCTCCGTCCTCTACAAGGATGACATGGGGGTGCCCACCCTCAAGTACCATTACGAGGGC" +
	"ATGAGCGTGGCAGAGTGTGGGTGCAGGTAG"

var gdf2prot string =
	"MCPGALWVALPLLSLLAGSLQGKPLQSWGRGSAGGNAHSPLGVPGGGLPEHTFNLKMFLENVKVDFLRSL" +
	"NLSGVPSQDKTRVEPPQYMIDLYNRYTSDKSTTPASNIVRSFSMEDAISITATEDFPFQKHILLFNISIP" +
	"RHEQITRAELRLYVSCQNHVDPSHDLKGSVVIYDVLDGTDAWDSATETKTFLVSQDIQDEGWETLEVSSA" +
	"VKRWVRSDSTKSKNKLEVTVESHRKGCDTLDISVPPGSRNLPFFVVFSNDHSSGTKETRLELREMISHEQ" +
	"ESVLKKLSKDGSTEAGESSHEEDTDGHVAAGSTLARRKRSAGAGSHCQKTSLRVNFEDIGWDSWIIAPKE" +
	"YEAYECKGGCFFPLADDVTPTKHAIVQTLVHLKFPTKVGKACCVPTKLSPISVLYKDDMGVPTLKYHYEG" +
	"MSVAECGCR*"

var gdf2LongProt string =
	"MetCysProGlyAlaLeuTrpValAlaLeuProLeuLeuSerLeuLeuAlaGlySerLeu" +
	"GlnGlyLysProLeuGlnSerTrpGlyArgGlySerAlaGlyGlyAsnAlaHisSerPro" +
	"LeuGlyValProGlyGlyGlyLeuProGluHisThrPheAsnLeuLysMetPheLeuGlu" +
	"AsnValLysValAspPheLeuArgSerLeuAsnLeuSerGlyValProSerGlnAspLys" +
	"ThrArgValGluProProGlnTyrMetIleAspLeuTyrAsnArgTyrThrSerAspLys" +
	"SerThrThrProAlaSerAsnIleValArgSerPheSerMetGluAspAlaIleSerIle" +
	"ThrAlaThrGluAspPheProPheGlnLysHisIleLeuLeuPheAsnIleSerIlePro" +
	"ArgHisGluGlnIleThrArgAlaGluLeuArgLeuTyrValSerCysGlnAsnHisVal" +
	"AspProSerHisAspLeuLysGlySerValValIleTyrAspValLeuAspGlyThrAsp" +
	"AlaTrpAspSerAlaThrGluThrLysThrPheLeuValSerGlnAspIleGlnAspGlu" +
	"GlyTrpGluThrLeuGluValSerSerAlaValLysArgTrpValArgSerAspSerThr" +
	"LysSerLysAsnLysLeuGluValThrValGluSerHisArgLysGlyCysAspThrLeu" +
	"AspIleSerValProProGlySerArgAsnLeuProPhePheValValPheSerAsnAsp" +
	"HisSerSerGlyThrLysGluThrArgLeuGluLeuArgGluMetIleSerHisGluGln" +
	"GluSerValLeuLysLysLeuSerLysAspGlySerThrGluAlaGlyGluSerSerHis" +
	"GluGluAspThrAspGlyHisValAlaAlaGlySerThrLeuAlaArgArgLysArgSer" +
	"AlaGlyAlaGlySerHisCysGlnLysThrSerLeuArgValAsnPheGluAspIleGly" +
	"TrpAspSerTrpIleIleAlaProLysGluTyrGluAlaTyrGluCysLysGlyGlyCys" +
	"PhePheProLeuAlaAspAspValThrProThrLysHisAlaIleValGlnThrLeuVal" +
	"HisLeuLysPheProThrLysValGlyLysAlaCysCysValProThrLysLeuSerPro" +
	"IleSerValLeuTyrLysAspAspMetGlyValProThrLeuLysTyrHisTyrGluGly" +
	"MetSerValAlaGluCysGlyCysArgTer"

func TestStress(t *testing.T) {
	expectedProt := ShortStringToPolypeptide(gdf2prot)
	actualProt, err := TranslateSeq(StringToBases(gdf2mrna))
	if !equal(actualProt, expectedProt) || err != nil {
		t.Errorf("stress test failed, expected\n%v, got\n%v", expectedProt, actualProt)
	}
	actualLongProt, err := TranslateToString(StringToBases(gdf2mrna))
	if actualLongProt != gdf2LongProt || err != nil {
		t.Errorf("stress test failed, expected\n%v, got\n%v", gdf2LongProt, actualLongProt)
	}
}

func equal(alpha []AminoAcid, beta []AminoAcid) bool {
	if len(alpha) != len(beta) {
		return false
	}
	for i := range alpha {
		if alpha[i] != beta[i] {
			return false
		}
	}
	return true
}

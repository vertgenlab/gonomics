package variant

import (
	"github.com/vertgenlab/gonomics/dna"
	"testing"
)

type EffectorTest struct {
	Effector
	ExpectedEffect CodingChange
}

func equalEffects(a, b CodingChange) bool {
	switch {
	case a.Type != b.Type:
		return false

	case a.CodingPos != b.CodingPos:
		return false

	case a.ProteinPos != b.ProteinPos:
		return false

	case dna.PolypeptideToString(a.AddedAa) !=
		dna.PolypeptideToString(b.AddedAa):
		return false

	case dna.PolypeptideToString(a.RemovedAa) !=
		dna.PolypeptideToString(b.RemovedAa):
		return false

	default:
		return true
	}
}

var substitutionEffectTests = []EffectorTest{
	{subTest1, subTest1ExpEff},
	{subTest2, subTest2ExpEff},
	{subTest3, subTest3ExpEff},
	{subTest4, subTest4ExpEff},
	{subTest5, subTest5ExpEff},
}

func TestSubstitutionEffect(t *testing.T) {
	for _, test := range substitutionEffectTests {
		actual, err := test.Effect(ref()[2:], -2, 0)
		if !equalEffects(actual, test.ExpectedEffect) {
			t.Errorf("problem with substitution effects.\n"+
				"Expect: %v\n"+
				"Actual: %v\n"+
				"Error: %s\n",
				test.ExpectedEffect, actual, err)
		}
	}
}

var insertionEffectTests = []EffectorTest{
	{insTest1, insTest1ExpEff},
	{insTest2, insTest2ExpEff},
	{insTest3, insTest3ExpEff},
	{insTest4, insTest4ExpEff},
	{insTest5, insTest5ExpEff},
}

func TestInsertionEffect(t *testing.T) {
	for _, test := range insertionEffectTests {
		actual, err := test.Effect(ref()[2:], -2, 0)
		if !equalEffects(actual, test.ExpectedEffect) {
			t.Errorf("problem with insertion effects.\n"+
				"Expect: %v\n"+
				"Actual: %v\n"+
				"Error: %s\n",
				test.ExpectedEffect, actual, err)
		}
	}
}

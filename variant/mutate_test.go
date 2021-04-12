package variant

import (
	"github.com/vertgenlab/gonomics/dna"
	"testing"
)

type MutatorTest struct {
	Mutator
	ExpectedSeq []dna.Base
}

var substitutionMutateTests = []MutatorTest{
	{subTest1, subTest1ExpSeq},
	{subTest2, subTest2ExpSeq},
	{subTest3, subTest3ExpSeq},
	{subTest4, subTest4ExpSeq},
	{subTest5, subTest5ExpSeq},
}

func TestSubstitutionMutate(t *testing.T) {
	for _, test := range substitutionMutateTests {
		actual, err := test.Mutate(ref(), 0, 0)
		if dna.CompareSeqsCaseSensitive(actual, test.ExpectedSeq) != 0 ||
			err != nil {
			t.Errorf("problem mutate via substitution.\n"+
				"Expect: %s\n"+
				"Actual: %s\n"+
				"Error: %s\n",
				dna.BasesToString(test.ExpectedSeq),
				dna.BasesToString(actual), err)
		}
	}
}

var insertionMutateTests = []MutatorTest{
	{insTest1, insTest1ExpSeq},
	{insTest2, insTest2ExpSeq},
	{insTest3, insTest3ExpSeq},
	{insTest4, insTest4ExpSeq},
	{insTest5, insTest5ExpSeq},
}

func TestInsertionMutate(t *testing.T) {
	for _, test := range insertionMutateTests {
		actual, err := test.Mutate(ref(), 0, 0)
		if dna.CompareSeqsCaseSensitive(actual, test.ExpectedSeq) != 0 ||
			err != nil {
			t.Errorf("problem mutate via insertion.\n"+
				"Expect: %s\n"+
				"Actual: %s\n"+
				"Error: %s\n",
				dna.BasesToString(test.ExpectedSeq),
				dna.BasesToString(actual), err)
		}
	}
}

var deletionMutateTests = []MutatorTest{
	{delTest1, delTest1ExpSeq},
	{delTest2, delTest2ExpSeq},
	{delTest3, delTest3ExpSeq},
	{delTest4, delTest4ExpSeq},
	{delTest5, delTest5ExpSeq},
}

func TestDeletionMutate(t *testing.T) {
	for _, test := range deletionMutateTests {
		actual, err := test.Mutate(ref(), 0, 0)
		if dna.CompareSeqsCaseSensitive(actual, test.ExpectedSeq) != 0 ||
			err != nil {
			t.Errorf("problem mutate via deletion.\n"+
				"Expect: %s\n"+
				"Actual: %s\n"+
				"Error: %s\n",
				dna.BasesToString(test.ExpectedSeq),
				dna.BasesToString(actual), err)
		}
	}
}

var delinsMutateTests = []MutatorTest{
	{delinsTest1, delinsTest1ExpSeq},
	{delinsTest2, delinsTest2ExpSeq},
	{delinsTest3, delinsTest3ExpSeq},
	{delinsTest4, delinsTest4ExpSeq},
	{delinsTest5, delinsTest5ExpSeq},
}

func TestDelinsMutate(t *testing.T) {
	for _, test := range delinsMutateTests {
		actual, err := test.Mutate(ref(), 0, 0)
		if dna.CompareSeqsCaseSensitive(actual, test.ExpectedSeq) != 0 ||
			err != nil {
			t.Errorf("problem mutate via delins.\n"+
				"Expect: %s\n"+
				"Actual: %s\n"+
				"Error: %s\n",
				dna.BasesToString(test.ExpectedSeq),
				dna.BasesToString(actual), err)
		}
	}
}

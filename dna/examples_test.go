package dna

import (
	"fmt"
)

func ExampleStringToBases() {
	var stringSeq string
	stringSeq = "ACGT"
	fmt.Println(stringSeq)

	var baseSeq []Base
	baseSeq = StringToBases(stringSeq)
	fmt.Println(baseSeq)
	// Output:
	// ACGT
	// [0 1 2 3]
}

func ExampleBasesToString() {
	var baseSeq []Base
	baseSeq = []Base{A, C, G, T}
	fmt.Println(baseSeq)

	var stringSeq string
	stringSeq = BasesToString(baseSeq)
	fmt.Println(stringSeq)
	// Output:
	// [0 1 2 3]
	// ACGT
}

func ExampleReverseComplement() {
	var baseSeq []Base
	baseSeq = []Base{A, T, G}
	fmt.Println(BasesToString(baseSeq))

	// Reverse complement modifies the slice in place so no return value
	ReverseComplement(baseSeq)

	fmt.Println(BasesToString(baseSeq))
	// Output:
	// ATG
	// CAT
}

func ExampleComplement() {
	var baseSeq []Base
	baseSeq = []Base{A, T, G}
	fmt.Println(BasesToString(baseSeq))

	// Complement modifies the slice in place so no return value
	Complement(baseSeq)

	fmt.Println(BasesToString(baseSeq))
	// Output:
	// ATG
	// TAC
}

func ExampleCountBase() {
	var seq []Base
	seq = []Base{A, A, C, T, T, T}

	fmt.Println(CountBase(seq, A))
	fmt.Println(CountBase(seq, C))
	fmt.Println(CountBase(seq, G))
	fmt.Println(CountBase(seq, T))
	fmt.Println(CountBase(seq, N))
	// Output:
	// 2
	// 1
	// 0
	// 3
	// 0
}

package variant

import "github.com/vertgenlab/gonomics/dna"

// Reference sequence used in the following tests
//
// Seq Idx : 0  1   2  3  4   5  6  7   8  9  10  11 12 13  14 15 16  17 18 19  20 21 22 23
// CDS Idx :        0  1  2   3  4  5   6  7  8   9  10 11  12 13 14  15 16 17
// Seq     : C  A   A  T  G   C  A  A   G  T  A   T  T  C   A  G  C   T  A  A   A  T  G  A
// Codons  :        -------   -------   -------   -------   -------   -------
// Protein :          Met       Gln       Val       Phe       Ser       Ter
// Prot Idx:           0         1         2         3         4         5

func ref() []dna.Base {
	return dna.StringToBases("CAATGCAAGTATTCAGCTAAATGA")
}

// Substitution Test 1: Internal Missense
// Seq Idx : 0  1   2  3  4   5  6  7   8  9  10  11 12 13  14 15 16  17 18 19  20 21 22 23
// CDS Idx :        0  1  2   3  4  5   6  7  8   9  10 11  12 13 14  15 16 17
// Seq     : C  A   A  T  G   C  A  A   G  T  A   T  T  C   A  G  C   T  A  A   A  T  G  A
// ChangeTo:                            T
// Codons  :        -------   -------   -------   -------   -------   -------
// Protein :          Met       Gln       Val       Phe       Ser       Ter
// Prot Idx:           0         1         2         3         4         5
var subTest1 = Substitution{
	Pos: 8,
	Ref: dna.G,
	Alt: dna.T,
}
var subTest1ExpSeq = dna.StringToBases("CAATGCAATTATTCAGCTAAATGA")
var subTest1ExpEff = CodingChange{
	CodingPos:  6,
	ProteinPos: 2,
	RemovedAa:  []dna.AminoAcid{dna.Val},
	AddedAa:    []dna.AminoAcid{dna.Leu},
	Type:       Missense,
}

// Substitution Test 2: Missense at stop
// Seq Idx : 0  1   2  3  4   5  6  7   8  9  10  11 12 13  14 15 16  17 18 19  20 21 22 23
// CDS Idx :        0  1  2   3  4  5   6  7  8   9  10 11  12 13 14  15 16 17
// Seq     : C  A   A  T  G   C  A  A   G  T  A   T  T  C   A  G  C   T  A  A   A  T  G  A
// ChangeTo:                                                                C
// Codons  :        -------   -------   -------   -------   -------   -------
// Protein :          Met       Gln       Val       Phe       Ser       Ter
// Prot Idx:           0         1         2         3         4         5
var subTest2 = Substitution{
	Pos: 19,
	Ref: dna.A,
	Alt: dna.C,
}
var subTest2ExpSeq = dna.StringToBases("CAATGCAAGTATTCAGCTACATGA")
var subTest2ExpEff = CodingChange{
	CodingPos:  17,
	ProteinPos: 5,
	RemovedAa:  []dna.AminoAcid{dna.Stop},
	AddedAa:    []dna.AminoAcid{dna.Tyr},
	Type:       Missense,
}

// Substitution Test 3: Missense at start
// Seq Idx : 0  1   2  3  4   5  6  7   8  9  10  11 12 13  14 15 16  17 18 19  20 21 22 23
// CDS Idx :        0  1  2   3  4  5   6  7  8   9  10 11  12 13 14  15 16 17
// Seq     : C  A   A  T  G   C  A  A   G  T  A   T  T  C   A  G  C   T  A  A   A  T  G  A
// ChangeTo:        T
// Codons  :        -------   -------   -------   -------   -------   -------
// Protein :          Met       Gln       Val       Phe       Ser       Ter
// Prot Idx:           0         1         2         3         4         5
var subTest3 = Substitution{
	Pos: 2,
	Ref: dna.A,
	Alt: dna.T,
}
var subTest3ExpSeq = dna.StringToBases("CATTGCAAGTATTCAGCTAAATGA")
var subTest3ExpEff = CodingChange{
	CodingPos:  0,
	ProteinPos: 0,
	RemovedAa:  []dna.AminoAcid{dna.Met},
	AddedAa:    []dna.AminoAcid{dna.Leu},
	Type:       Missense,
}

// Substitution Test 4: Nonsense
// Seq Idx : 0  1   2  3  4   5  6  7   8  9  10  11 12 13  14 15 16  17 18 19  20 21 22 23
// CDS Idx :        0  1  2   3  4  5   6  7  8   9  10 11  12 13 14  15 16 17
// Seq     : C  A   A  T  G   C  A  A   G  T  A   T  T  C   A  G  C   T  A  A   A  T  G  A
// ChangeTo:                  T
// Codons  :        -------   -------   -------   -------   -------   -------
// Protein :          Met       Gln       Val       Phe       Ser       Ter
// Prot Idx:           0         1         2         3         4         5
var subTest4 = Substitution{
	Pos: 5,
	Ref: dna.C,
	Alt: dna.T,
}
var subTest4ExpSeq = dna.StringToBases("CAATGTAAGTATTCAGCTAAATGA")
var subTest4ExpEff = CodingChange{
	CodingPos:  3,
	ProteinPos: 1,
	RemovedAa:  []dna.AminoAcid{dna.Gln},
	AddedAa:    []dna.AminoAcid{dna.Stop},
	Type:       Nonsense,
}

// Substitution Test 5: Silent
// Seq Idx : 0  1   2  3  4   5  6  7   8  9  10  11 12 13  14 15 16  17 18 19  20 21 22 23
// CDS Idx :        0  1  2   3  4  5   6  7  8   9  10 11  12 13 14  15 16 17
// Seq     : C  A   A  T  G   C  A  A   G  T  A   T  T  C   A  G  C   T  A  A   A  T  G  A
// ChangeTo:                                            T
// Codons  :        -------   -------   -------   -------   -------   -------
// Protein :          Met       Gln       Val       Phe       Ser       Ter
// Prot Idx:           0         1         2         3         4         5
var subTest5 = Substitution{
	Pos: 13,
	Ref: dna.C,
	Alt: dna.T,
}
var subTest5ExpSeq = dna.StringToBases("CAATGCAAGTATTTAGCTAAATGA")
var subTest5ExpEff = CodingChange{
	CodingPos:  11,
	ProteinPos: 3,
	RemovedAa:  []dna.AminoAcid{},
	AddedAa:    []dna.AminoAcid{},
	Type:       Silent,
}

// Insertion Test 1: Internal InFrameInsertion
// Seq Idx : 0  1   2  3  4   5  6  7   8  9  10  11 12 13  14 15 16  17 18 19  20 21 22 23
// CDS Idx :        0  1  2   3  4  5   6  7  8   9  10 11  12 13 14  15 16 17
// Seq     : C  A   A  T  G   C  A  A   G  T  A   T  T  C   A  G  C   T  A  A   A  T  G  A
// Ins Pos :                          |
// Ins Seq :                         AGG
// Codons  :        -------   -------   -------   -------   -------   -------
// Protein :          Met       Gln       Val       Phe       Ser       Ter
// Prot Idx:           0         1         2         3         4         5
var insTest1 = Insertion{
	Pos: 8,
	Seq: dna.StringToBases("AGG"),
}
var insTest1ExpSeq = dna.StringToBases("CAATGCAAAGGGTATTCAGCTAAATGA")
var insTest1ExpEff = CodingChange{
	CodingPos:  6,
	ProteinPos: 2,
	RemovedAa:  []dna.AminoAcid{},
	AddedAa:    []dna.AminoAcid{dna.Arg},
	Type:       InFrameInsertion,
}

// Insertion Test 2: Frameshift resulting in late termination
// Seq Idx : 0  1   2  3  4   5  6  7   8  9  10  11 12 13  14 15 16  17 18 19  20 21 22 23
// CDS Idx :        0  1  2   3  4  5   6  7  8   9  10 11  12 13 14  15 16 17
// Seq     : C  A   A  T  G   C  A  A   G  T  A   T  T  C   A  G  C   T  A  A   A  T  G  A
// Ins Pos :                                                  |
// Ins Seq :                                                  GT
// Codons  :        -------   -------   -------   -------   -------   -------
// Protein :          Met       Gln       Val       Phe       Ser       Ter
// Prot Idx:           0         1         2         3         4         5
var insTest2 = Insertion{
	Pos: 15,
	Seq: dna.StringToBases("GT"),
}
var insTest2ExpSeq = dna.StringToBases("CAATGCAAGTATTCAGTGCTAAATGA")
var insTest2ExpEff = CodingChange{
	CodingPos:  13,
	ProteinPos: 5, // insertion creates degenerate codon so pos should be incremented
	RemovedAa:  []dna.AminoAcid{dna.Stop},
	AddedAa:    []dna.AminoAcid{dna.Ala, dna.Lys, dna.Stop},
	Type:       Frameshift,
}

// Insertion Test 3: Frameshift with no stop found
// Seq Idx : 0  1   2  3  4   5  6  7   8  9  10  11 12 13  14 15 16  17 18 19  20 21 22 23
// CDS Idx :        0  1  2   3  4  5   6  7  8   9  10 11  12 13 14  15 16 17
// Seq     : C  A   A  T  G   C  A  A   G  T  A   T  T  C   A  G  C   T  A  A   A  T  G  A
// Ins Pos :                          |
// Ins Seq :                          A
// Codons  :        -------   -------   -------   -------   -------   -------
// Protein :          Met       Gln       Val       Phe       Ser       Ter
// Prot Idx:           0         1         2         3         4         5
var insTest3 = Insertion{
	Pos: 8,
	Seq: dna.StringToBases("A"),
}
var insTest3ExpSeq = dna.StringToBases("CAATGCAAAGTATTCAGCTAAATGA")
var insTest3ExpEff = CodingChange{
	CodingPos:  6,
	ProteinPos: 2,
	RemovedAa:  []dna.AminoAcid{dna.Val, dna.Phe, dna.Ser, dna.Stop},
	AddedAa:    []dna.AminoAcid{dna.Ser, dna.Ile, dna.Gln, dna.Leu, dna.Asn},
	Type:       Frameshift,
}

// Insertion Test 4: Insert premature stop
// Seq Idx : 0  1   2  3  4   5  6  7   8  9  10  11 12 13  14 15 16  17 18 19  20 21 22 23
// CDS Idx :        0  1  2   3  4  5   6  7  8   9  10 11  12 13 14  15 16 17
// Seq     : C  A   A  T  G   C  A  A   G  T  A   T  T  C   A  G  C   T  A  A   A  T  G  A
// Ins Pos :                          |
// Ins Seq :                         TAG
// Codons  :        -------   -------   -------   -------   -------   -------
// Protein :          Met       Gln       Val       Phe       Ser       Ter
// Prot Idx:           0         1         2         3         4         5
var insTest4 = Insertion{
	Pos: 8,
	Seq: dna.StringToBases("TAG"),
}
var insTest4ExpSeq = dna.StringToBases("CAATGCAATAGGTATTCAGCTAAATGA")
var insTest4ExpEff = CodingChange{
	CodingPos:  6,
	ProteinPos: 2,
	RemovedAa:  []dna.AminoAcid{},
	AddedAa:    []dna.AminoAcid{dna.Stop},
	Type:       InFrameInsertion,
}

// Insertion Test 5: Insertion disrupts start
// Seq Idx : 0  1   2  3  4   5  6  7   8  9  10  11 12 13  14 15 16  17 18 19  20 21 22 23
// CDS Idx :        0  1  2   3  4  5   6  7  8   9  10 11  12 13 14  15 16 17
// Seq     : C  A   A  T  G   C  A  A   G  T  A   T  T  C   A  G  C   T  A  A   A  T  G  A
// Ins Pos :             |
// Ins Seq :           CTACCC
// Codons  :        -------   -------   -------   -------   -------   -------
// Protein :          Met       Gln       Val       Phe       Ser       Ter
// Prot Idx:           0         1         2         3         4         5
var insTest5 = Insertion{
	Pos: 4,
	Seq: dna.StringToBases("CTACCC"),
}
var insTest5ExpSeq = dna.StringToBases("CAATCTACCCGCAAGTATTCAGCTAAATGA")
var insTest5ExpEff = CodingChange{
	CodingPos:  2,
	ProteinPos: 0,
	RemovedAa:  []dna.AminoAcid{dna.Met},
	AddedAa:    []dna.AminoAcid{dna.Ile, dna.Tyr, dna.Pro},
	Type:       InFrameInsertion,
}

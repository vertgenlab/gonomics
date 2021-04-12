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

// Deletion Test 1: Frameshift results in termination after original stop
// Seq Idx : 0  1   2  3  4   5  6  7   8  9  10  11 12 13  14 15 16  17 18 19  20 21 22 23
// CDS Idx :        0  1  2   3  4  5   6  7  8   9  10 11  12 13 14  15 16 17
// Seq     : C  A   A  T  G   C  A  A   G  T  A   T  T  C   A  G  C   T  A  A   A  T  G  A
// Deleted :												   *
// Codons  :        -------   -------   -------   -------   -------   -------
// Protein :          Met       Gln       Val       Phe       Ser       Ter
// Prot Idx:           0         1         2         3         4         5
var delTest1 = Deletion{
	Start: 15,
	End:   16,
}
var delTest1ExpSeq = dna.StringToBases("CAATGCAAGTATTCACTAAATGA")
var delTest1ExpEff = CodingChange{
	CodingPos:  13,
	ProteinPos: 4,
	RemovedAa:  []dna.AminoAcid{dna.Ser, dna.Stop},
	AddedAa:    []dna.AminoAcid{dna.Thr, dna.Lys, dna.Stop},
	Type:       Frameshift,
}

// Deletion Test 2: Frameshift with no new stop codon
// Seq Idx : 0  1   2  3  4   5  6  7   8  9  10  11 12 13  14 15 16  17 18 19  20 21 22 23
// CDS Idx :        0  1  2   3  4  5   6  7  8   9  10 11  12 13 14  15 16 17
// Seq     : C  A   A  T  G   C  A  A   G  T  A   T  T  C   A  G  C   T  A  A   A  T  G  A
// Deleted :                                                   ****
// Codons  :        -------   -------   -------   -------   -------   -------
// Protein :          Met       Gln       Val       Phe       Ser       Ter
// Prot Idx:           0         1         2         3         4         5
var delTest2 = Deletion{
	Start: 15,
	End:   17,
}
var delTest2ExpSeq = dna.StringToBases("CAATGCAAGTATTCATAAATGA")
var delTest2ExpEff = CodingChange{
	CodingPos:  13,
	ProteinPos: 4,
	RemovedAa:  []dna.AminoAcid{dna.Ser, dna.Stop},
	AddedAa:    []dna.AminoAcid{dna.Ile, dna.Asn},
	Type:       Frameshift,
}

// Deletion Test 3: InFrameDeletion
// Seq Idx : 0  1   2  3  4   5  6  7   8  9  10  11 12 13  14 15 16  17 18 19  20 21 22 23
// CDS Idx :        0  1  2   3  4  5   6  7  8   9  10 11  12 13 14  15 16 17
// Seq     : C  A   A  T  G   C  A  A   G  T  A   T  T  C   A  G  C   T  A  A   A  T  G  A
// Deleted :                  *******
// Codons  :        -------   -------   -------   -------   -------   -------
// Protein :          Met       Gln       Val       Phe       Ser       Ter
// Prot Idx:           0         1         2         3         4         5
var delTest3 = Deletion{
	Start: 5,
	End:   8,
}
var delTest3ExpSeq = dna.StringToBases("CAATGGTATTCAGCTAAATGA")
var delTest3ExpEff = CodingChange{
	CodingPos:  3,
	ProteinPos: 1,
	RemovedAa:  []dna.AminoAcid{dna.Gln},
	AddedAa:    []dna.AminoAcid{},
	Type:       InFrameDeletion,
}

// Deletion Test 4: InFrameDeletion disrupts a codon to make degenerate codon
// Seq Idx : 0  1   2  3  4   5  6  7   8  9  10  11 12 13  14 15 16  17 18 19  20 21 22 23
// CDS Idx :        0  1  2   3  4  5   6  7  8   9  10 11  12 13 14  15 16 17
// Seq     : C  A   A  T  G   C  A  A   G  T  A   T  T  C   A  G  C   T  A  A   A  T  G  A
// Deleted :                        ********
// Codons  :        -------   -------   -------   -------   -------   -------
// Protein :          Met       Gln       Val       Phe       Ser       Ter
// Prot Idx:           0         1         2         3         4         5
var delTest4 = Deletion{
	Start: 7,
	End:   10,
}
var delTest4ExpSeq = dna.StringToBases("CAATGCAATTCAGCTAAATGA")
var delTest4ExpEff = CodingChange{
	CodingPos:  5,
	ProteinPos: 2,
	RemovedAa:  []dna.AminoAcid{dna.Val},
	AddedAa:    []dna.AminoAcid{},
	Type:       InFrameDeletion,
}

// Deletion Test 5: InFrameDeletion disrupts codon to make nonsynonymous codon
// Seq Idx : 0  1   2  3  4   5  6  7   8  9  10  11 12 13  14 15 16  17 18 19  20 21 22 23
// CDS Idx :        0  1  2   3  4  5   6  7  8   9  10 11  12 13 14  15 16 17
// Seq     : C  A   A  T  G   C  A  A   G  T  A   T  T  C   A  G  C   T  A  A   A  T  G  A
// Deleted :                               ******************
// Codons  :        -------   -------   -------   -------   -------   -------
// Protein :          Met       Gln       Val       Phe       Ser       Ter
// Prot Idx:           0         1         2         3         4         5
var delTest5 = Deletion{
	Start: 9,
	End:   15,
}
var delTest5ExpSeq = dna.StringToBases("CAATGCAAGGCTAAATGA")
var delTest5ExpEff = CodingChange{
	CodingPos:  7,
	ProteinPos: 2,
	RemovedAa:  []dna.AminoAcid{dna.Val, dna.Phe, dna.Ser},
	AddedAa:    []dna.AminoAcid{dna.Gly},
	Type:       InFrameDeletion,
}

// Delins Test 1: Codon swap
// Seq Idx : 0  1   2  3  4   5  6  7   8  9  10  11 12 13  14 15 16  17 18 19  20 21 22 23
// CDS Idx :        0  1  2   3  4  5   6  7  8   9  10 11  12 13 14  15 16 17
// Seq     : C  A   A  T  G   C  A  A   G  T  A   T  T  C   A  G  C   T  A  A   A  T  G  A
// Deleted :                            *******
// Inserted:                            C  A  T
// Codons  :        -------   -------   -------   -------   -------   -------
// Protein :          Met       Gln       Val       Phe       Ser       Ter
// Prot Idx:           0         1         2         3         4         5
var delinsTest1 = Delins{
	Start:  8,
	End:    11,
	InsSeq: dna.StringToBases("CAT"),
}
var delinsTest1ExpSeq = dna.StringToBases("CAATGCAACATTTCAGCTAAATGA")
var delinsTest1ExpEff = CodingChange{
	CodingPos:  6,
	ProteinPos: 2,
	RemovedAa:  []dna.AminoAcid{dna.Val},
	AddedAa:    []dna.AminoAcid{dna.His},
	Type:       Missense,
}

// Delins Test 2: Deletion plus insertion of a degenerate codon
// Seq Idx : 0  1   2  3  4   5  6  7   8  9  10  11 12 13  14 15 16  17 18 19  20 21 22 23
// CDS Idx :        0  1  2   3  4  5   6  7  8   9  10 11  12 13 14  15 16 17
// Seq     : C  A   A  T  G   C  A  A   G  T  A   T  T  C   A  G  C   T  A  A   A  T  G  A
// Deleted :                        ***********
// Inserted:                        G
// Codons  :        -------   -------   -------   -------   -------   -------
// Protein :          Met       Gln       Val       Phe       Ser       Ter
// Prot Idx:           0         1         2         3         4         5
var delinsTest2 = Delins{
	Start:  7,
	End:    11,
	InsSeq: dna.StringToBases("G"),
}
var delinsTest2ExpSeq = dna.StringToBases("CAATGCAGTTCAGCTAAATGA")
var delinsTest2ExpEff = CodingChange{
	CodingPos:  5,
	ProteinPos: 2,
	RemovedAa:  []dna.AminoAcid{dna.Val},
	AddedAa:    []dna.AminoAcid{},
	Type:       InFrameDeletion,
}

// Delins Test 3: Delins results in InFrameInsertion
// Seq Idx : 0  1   2  3  4   5  6  7   8  9  10  11 12 13  14 15 16  17 18 19  20 21 22 23
// CDS Idx :        0  1  2   3  4  5   6  7  8   9  10 11  12 13 14  15 16 17
// Seq     : C  A   A  T  G   C  A  A   G  T  A   T  T  C   A  G  C   T  A  A   A  T  G  A
// Deleted :           ****
// Inserted:           ACAGT
// Codons  :        -------   -------   -------   -------   -------   -------
// Protein :          Met       Gln       Val       Phe       Ser       Ter
// Prot Idx:           0         1         2         3         4         5
var delinsTest3 = Delins{
	Start:  3,
	End:    5,
	InsSeq: dna.StringToBases("ACAGT"),
}
var delinsTest3ExpSeq = dna.StringToBases("CAAACAGTCAAGTATTCAGCTAAATGA")
var delinsTest3ExpEff = CodingChange{
	CodingPos:  1,
	ProteinPos: 0,
	RemovedAa:  []dna.AminoAcid{dna.Met},
	AddedAa:    []dna.AminoAcid{dna.Asn, dna.Ser},
	Type:       InFrameInsertion,
}

// Delins Test 4: Frameshift resulting in premature termination
// Seq Idx : 0  1   2  3  4   5  6  7   8  9  10  11 12 13  14 15 16  17 18 19  20 21 22 23
// CDS Idx :        0  1  2   3  4  5   6  7  8   9  10 11  12 13 14  15 16 17
// Seq     : C  A   A  T  G   C  A  A   G  T  A   T  T  C   A  G  C   T  A  A   A  T  G  A
// Deleted :                                  ***********
// Inserted:                                      GTG
// Codons  :        -------   -------   -------   -------   -------   -------
// Protein :          Met       Gln       Val       Phe       Ser       Ter
// Prot Idx:           0         1         2         3         4         5
var delinsTest4 = Delins{
	Start:  10,
	End:    14,
	InsSeq: dna.StringToBases("GTG"),
}
var delinsTest4ExpSeq = dna.StringToBases("CAATGCAAGTGTGAGCTAAATGA")
var delinsTest4ExpEff = CodingChange{
	CodingPos:  8,
	ProteinPos: 3,
	RemovedAa:  []dna.AminoAcid{dna.Phe, dna.Ser, dna.Stop},
	AddedAa:    []dna.AminoAcid{dna.Stop},
	Type:       Frameshift,
}

// Delins Test 5: Frameshift resulting in late termination
// Seq Idx : 0  1   2  3  4   5  6  7   8  9  10  11 12 13  14 15 16  17 18 19  20 21 22 23
// CDS Idx :        0  1  2   3  4  5   6  7  8   9  10 11  12 13 14  15 16 17
// Seq     : C  A   A  T  G   C  A  A   G  T  A   T  T  C   A  G  C   T  A  A   A  T  G  A
// Deleted :                                         ****
// Inserted:                                         ATAA
// Codons  :        -------   -------   -------   -------   -------   -------
// Protein :          Met       Gln       Val       Phe       Ser       Ter
// Prot Idx:           0         1         2         3         4         5
var delinsTest5 = Delins{
	Start:  12,
	End:    14,
	InsSeq: dna.StringToBases("ATAA"),
}
var delinsTest5ExpSeq = dna.StringToBases("CAATGCAAGTATATAAAGCTAAATGA")
var delinsTest5ExpEff = CodingChange{
	CodingPos:  10,
	ProteinPos: 3,
	RemovedAa:  []dna.AminoAcid{dna.Phe, dna.Ser, dna.Stop},
	AddedAa:    []dna.AminoAcid{dna.Tyr, dna.Lys, dna.Ala, dna.Lys, dna.Stop},
	Type:       Frameshift,
}

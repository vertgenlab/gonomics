package genomeGraph

// func TestIndexGenomeIntoMap(t *testing.T) {
// 	genome := []Node{
// 		{Seq: dna.StringToBases("ACGT")},
// 		{Seq: dna.StringToBases("TGCA")},
// 	}

// 	seedLen := 1
// 	seedStep := 1

// 	expected := map[uint64][]uint64{
// 		dnaToNumber(dna.StringToBases("AC"), 0, 2): {ChromAndPosToNumber(0, 0)},
// 		dnaToNumber(dna.StringToBases("CG"), 1, 3): {ChromAndPosToNumber(0, 1)},
// 		dnaToNumber(dna.StringToBases("GT"), 2, 4): {ChromAndPosToNumber(0, 2)},
// 		dnaToNumber(dna.StringToBases("TG"), 0, 2): {ChromAndPosToNumber(1, 0)},
// 		dnaToNumber(dna.StringToBases("GC"), 1, 3): {ChromAndPosToNumber(1, 1)},
// 		dnaToNumber(dna.StringToBases("CA"), 2, 4): {ChromAndPosToNumber(1, 2)},
// 	}

// 	result := IndexGenomeIntoMap(genome, seedLen, seedStep)

// 	if !reflect.DeepEqual(result, expected) {
// 		t.Errorf("IndexGenomeIntoMap() = %v, want %v", result, expected)
// 	}
// }

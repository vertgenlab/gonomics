package fasta

func BinFasta(fa []Fasta, binNum int) {
	var totalBases int
	bins := make(map[int][]Fasta, binNum)
	contigs := len(fa)

	if contigs < binNum {
		totalBases = 0
		for f := range fa {
			totalBases = totalBases + len(fa[f].Seq)
		}
		//TODO: chance that there's a situation where the final bin with leftOvers is empty or very small, what I have with bases per equal bin,
		//	or chance that the last bin could be almost double the size of everything else (totalBases/BinNum and remainder goes in last bin)
		//	could do a check to see if remainder is super large and then handle as necessary
		equalBinNum := binNum - 1
		basesPerEqualBin := totalBases / equalBinNum
		leftOvers := totalBases % equalBinNum
	} else if contigs == binNum {

	} else if contigs > binNum {

	}

}

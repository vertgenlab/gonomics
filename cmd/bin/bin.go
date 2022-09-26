package main

import(
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/bin"
	"github.com/vertgenlab/gonomics/fileio"
)

func binGenome(genome string, binNum int, minSize int){
	records := fasta.Read(genome)
	var bins map[int][]fasta.Fasta

	bins = bin.BinGenomeNoBreaks(records, binNum, minSize)

	for i := range bins {
		var name string
		var thisContig []fasta.Fasta
		if len(bins[i]) == 1 {
			name = bins[i][0].Name
			thisContig = bins[i]
		} else { //file name = first_second_...
			for j := range bins[i] {
				name = name + "_" + bins[i][j].Name
			}
			thisContig = bins[i]
		}
		fasta.Write(name, thisContig)
	}

}
package main

import (
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/simpleGraph"
	"github.com/vertgenlab/gonomics/vcf"
)

func vcfToSimpleGraph(vcfFile, faFile string) *simpleGraph.SimpleGraph {
	ref := fasta.GoReadToChan(faFile)

	hashByChrom := make(map[string][]*vcf.Vcf)
	file := fileio.EasyOpen(vcfFile)
	defer file.Close()
	vcf.ReadHeader(file)

	for curr, done := vcf.NextVcf(file); !done; curr, done = vcf.NextVcf(file) {
		hashByChrom[curr.Chr] = append(hashByChrom[curr.Chr], curr)
	}

	return simpleGraph.VariantGraph(ref, hashByChrom)
}

package main

import (
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/simpleGraph"
	"github.com/vertgenlab/gonomics/vcf"
)

func vcfToSimpleGraph(vcfFile, faFile string) *simpleGraph.SimpleGraph {
	file := fileio.EasyOpen(vcfFile)
	vcf.ReadHeader(file)
	vcfChannel := make(chan *vcf.Vcf)
	ref := make(chan *fasta.Fasta)
	go vcf.ReadToChan(file, vcfChannel)
	go fasta.ReadToChan(faFile, ref)

	hashByChrom := make(map[string][]*vcf.Vcf)
	for i := range vcfChannel {
		hashByChrom[i.Chr] = append(hashByChrom[i.Chr], i)
	}
	return simpleGraph.VariantGraph(ref, hashByChrom)
}

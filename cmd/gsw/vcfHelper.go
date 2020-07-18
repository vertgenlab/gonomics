package main

import (
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/simpleGraph"
	"github.com/vertgenlab/gonomics/vcf"
)

//set up vcf records split my chroms
func makeVcfChrMap(vcfChannel <-chan *vcf.Vcf) map[string][]*vcf.Vcf {
	chrVcfMap := make(map[string][]*vcf.Vcf)
	for i := range vcfChannel {
		chrVcfMap[i.Chr] = append(chrVcfMap[i.Chr], i)
	}
	return chrVcfMap
}

func vcfToSimpleGraph(vcfFile, faFile string) *simpleGraph.SimpleGraph {
	file := fileio.EasyOpen(vcfFile)
	vcf.ReadHeader(file)
	vcfChannel := make(chan *vcf.Vcf)
	ref := make(chan *fasta.Fasta)
	go vcf.ReadToChan(file, vcfChannel)
	go fasta.ReadToChan(faFile, ref)
	hashByChrom := makeVcfChrMap(vcfChannel)
	return simpleGraph.VariantGraph(ref, hashByChrom)
}

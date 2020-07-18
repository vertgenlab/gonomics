package main

import (
	"github.com/vertgenlab/gonomics/axt"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/simpleGraph"
	"github.com/vertgenlab/gonomics/vcf"
)

func convertAxt(axtFile, format, targetFa, output string) {
	switch format {
	case "vcf":
		vcfChannel := goChannelAxtVcf(axtFile)

		file := fileio.EasyCreate(output)
		vcf.NewWriteHeader(file, vcf.NewHeader(targetFa))
		for i := range vcfChannel {
			vcf.WriteVcf(file, i)
		}
	case "gg":
		simpleGraph.Write(output, axtToSimpleGraph(axtFile, targetFa))
	default:

	}
}

func goChannelAxtVcf(axtFile string) <-chan *vcf.Vcf {
	ans := make(chan *vcf.Vcf)
	axtChannel := make(chan *axt.Axt)
	go axt.ReadToChan(fileio.EasyOpen(axtFile), axtChannel)
	go workThreadAxtVcf(axtChannel, ans)
	return ans
}

func axtToSimpleGraph(axtFile, faFile string) *simpleGraph.SimpleGraph {
	vcfChannel := goChannelAxtVcf(axtFile)
	chrVcfMap := makeVcfChrMap(vcfChannel)
	ref := make(chan *fasta.Fasta)
	go fasta.ReadToChan(faFile, ref)
	return simpleGraph.VariantGraph(ref, chrVcfMap)
}

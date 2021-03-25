package main

import (
	"github.com/vertgenlab/gonomics/axt"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/genomeGraph"
	"github.com/vertgenlab/gonomics/vcf"
	"path/filepath"
	"strings"
)

func convertAxt(axtFile, format, targetFa, output string) {
	switch format {
	case "vcf":
		vcfChannel := goChannelAxtVcf(axtFile)
		file := fileio.EasyCreate(output)
		header := vcf.NewHeader(strings.TrimSuffix(targetFa, filepath.Ext(targetFa)))
		vcf.NewWriteHeader(file, header)
		var ans []*vcf.Vcf
		for i := range vcfChannel {
			ans = append(ans, i)
		}
		vcf.Sort(ans)
		vcf.WriteVcfToFileHandle(file, ans)
		file.Close()
	case "gg":
		genomeGraph.Write(output, axtToGenomeGraph(axtFile, targetFa))
	default:
		ggToolsUsage()
		errorMessage()
	}
}

func goChannelAxtVcf(axtFile string) <-chan *vcf.Vcf {
	axtChannel := make(chan axt.Axt, 2408)
	go axt.ReadToChan(axtFile, axtChannel)

	ans := make(chan *vcf.Vcf, 2408)
	go workThreadAxtVcf(axtChannel, ans)
	return ans
}

func axtToGenomeGraph(axtFile, faFile string) *genomeGraph.GenomeGraph {
	vcfChannel := goChannelAxtVcf(axtFile)
	chrVcfMap := make(map[string][]*vcf.Vcf)
	for i := range vcfChannel {
		chrVcfMap[i.Chr] = append(chrVcfMap[i.Chr], i)
	}
	ref := fasta.GoReadToChan(faFile)
	var gg *genomeGraph.GenomeGraph = &genomeGraph.GenomeGraph{}
	gg = genomeGraph.VariantGraph(ref, chrVcfMap)
	return gg
}

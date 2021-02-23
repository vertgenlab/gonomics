package main

import (
	"github.com/vertgenlab/gonomics/axt"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/simpleGraph"
	"github.com/vertgenlab/gonomics/vcf"
	"path/filepath"
	"strings"
	"sync"
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
		simpleGraph.Write(output, axtToSimpleGraph(axtFile, targetFa))
	default:
		ggToolsUsage()
		errorMessage()
	}
}

func goChannelAxtVcf(axtFile string) <-chan *vcf.Vcf {
	file := fileio.EasyOpen(axtFile)
	var wg sync.WaitGroup
	axtChannel := make(chan *axt.Axt, 2408)
	wg.Add(1)
	go axt.ReadToChan(file, axtChannel, &wg)

	go func() {
		wg.Wait()
		close(axtChannel)
	}()

	ans := make(chan *vcf.Vcf, 2408)
	go workThreadAxtVcf(axtChannel, ans)
	return ans
}

func axtToSimpleGraph(axtFile, faFile string) *simpleGraph.SimpleGraph {
	vcfChannel := goChannelAxtVcf(axtFile)
	chrVcfMap := make(map[string][]*vcf.Vcf)
	for i := range vcfChannel {
		chrVcfMap[i.Chr] = append(chrVcfMap[i.Chr], i)
	}
	ref := fasta.GoReadToChan(faFile)
	var gg *simpleGraph.SimpleGraph = &simpleGraph.SimpleGraph{}
	gg = simpleGraph.VariantGraph(ref, chrVcfMap)
	return gg
}

package main

import (
	"github.com/vertgenlab/gonomics/axt"
	"github.com/vertgenlab/gonomics/chain"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/simpleGraph"
	"github.com/vertgenlab/gonomics/vcf"
	"strings"
)

func convertChains(chainFile, targetFa, queryFa, format, output string) {
	switch format {
	case "axt":
		axtChannel := goChainToAxt(chainFile, targetFa, queryFa)
		file := fileio.EasyCreate(output)
		var idx int = 0
		for i := range axtChannel {
			axt.WriteToFileHandle(file, i, idx)
		}
	case "vcf":
		vcfChannel := goChainToVcf(chainFile, targetFa, queryFa)
		file := fileio.EasyCreate(output)
		vcf.NewWriteHeader(file, vcf.NewHeader(targetFa))
		for i := range vcfChannel {
			vcf.WriteVcf(file, i)
		}
		file.Close()
	case "gg":
		simpleGraph.Write(output, chainToSimpleGraph(chainFile, targetFa, queryFa))
	default:
		ggToolsUsage()
		errorMessage()
	}
}

func goChainToAxt(chainFile, targetFa, queryFa string) <-chan *axt.Axt {
	target, query := fasta.Read(targetFa), fasta.Read(queryFa)
	chainFa := chain.GoReadSeqChain(chainFile, target, query)
	ans := make(chan *axt.Axt, 2408)
	go workThreadChainAxt(chainFa, ans)
	return ans
}

func goChainToVcf(chainFile, targetFa, queryFa string) <-chan *vcf.Vcf {
	ans := make(chan *vcf.Vcf, 2408)
	axtChannel := goChainToAxt(chainFile, targetFa, queryFa)
	go workThreadAxtVcf(axtChannel, ans)
	return ans
}

func chainToSimpleGraph(chainFile, targetFa, queryFa string) *simpleGraph.SimpleGraph {
	target, query := fasta.Read(targetFa), fasta.Read(queryFa)
	chainFa := chain.GoReadSeqChain(chainFile, target, query)
	axtChannel := make(chan *axt.Axt, 2408)
	go workThreadChainAxt(chainFa, axtChannel)

	vcfChannel := make(chan *vcf.Vcf, 2408)
	go workThreadAxtVcf(axtChannel, vcfChannel)

	chrVcfMap := make(map[string][]*vcf.Vcf)
	for i := range vcfChannel {
		chrVcfMap[i.Chr] = append(chrVcfMap[i.Chr], i)
	}
	//set up fa channel
	ref := goFaChannel(target)
	//return the simple graph
	return simpleGraph.VariantGraph(ref, chrVcfMap)
}

func workThreadChainAxt(chFa *chain.SeqChain, ans chan<- *axt.Axt) {
	for i := range chFa.Chains {
		ans <- chain.ChainToAxt(i, chFa.TSeq[i.TName], chFa.QSeq[i.QName])
	}
	close(ans)
}

func workThreadAxtVcf(axtChannel <-chan *axt.Axt, ans chan<- *vcf.Vcf) {
	var j int = 0
	var curr []*vcf.Vcf
	for i := range axtChannel {
		//filter for uniq
		curr = vcf.FilterVcfPos(axt.AxtToVcf(i))
		for j = 0; j < len(curr); j++ {
			if !strings.Contains(curr[j].Ref, "N") && !strings.Contains(curr[j].Alt[0], "N") {
				ans <- curr[j]
			}
		}
	}
	close(ans)
}

func goFaChannel(ref []*fasta.Fasta) <-chan *fasta.Fasta {
	//set up faChan
	ans := make(chan *fasta.Fasta, 100)
	go faWorker(ref, ans)
	return ans
}

func faWorker(ref []*fasta.Fasta, faChan chan<- *fasta.Fasta) {
	for _, chr := range ref {
		faChan <- chr
	}
	close(faChan)
}

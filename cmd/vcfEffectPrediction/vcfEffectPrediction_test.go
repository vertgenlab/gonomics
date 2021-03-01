package main

import (
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/gtf"
	"github.com/vertgenlab/gonomics/vcf"
	"strings"
	"sync"
	"testing"
)

func TestVcfEffectPrediction(t *testing.T) {
	settings := &Settings{
		Vcf:            "../../gtf/testdata/test.vcf",
		Gtf:            "../../gtf/testdata/test.gtf",
		Fasta:          "",
		Threads:        1,
		AllTranscripts: false,
		OutFile:        "",
	}

	krit := fasta.Read("../../gtf/testdata/krit1.fa")
	cftr := fasta.Read("../../gtf/testdata/cftr.fa")

	// The fasta files were split up to reduce the filesize. The following 6 lines
	// assemble the fasta files so that everything is in the correct place
	f := []*fasta.Fasta{{"chr7", make([]dna.Base, 92198968)}}
	f[0].Seq = append(f[0].Seq, krit[0].Seq...)
	for i := 0; i < 117480024-92246100; i++ {
		f[0].Seq = append(f[0].Seq, dna.N)
	}
	f[0].Seq = append(f[0].Seq, cftr[0].Seq...)

	fasta.AllToUpper(f)
	fastaRecords := fasta.ToMap(f)
	gtfRecords := gtf.Read(settings.Gtf)
	tree := gtf.GenesToIntervalTree(gtfRecords)

	vcfChan, _ := vcf.GoReadToChan(settings.Vcf)
	answer := make(chan *vcf.Vcf, 1000)
	var wg sync.WaitGroup

	for i := 0; i < settings.Threads; i++ {
		wg.Add(1)
		go annotationWorker(&wg, tree, fastaRecords, vcfChan, answer, settings.AllTranscripts)
	}

	go func() {
		wg.Wait()
		close(answer)
	}()

	var words, newWords, newerWords, annotation []string
	var correctCDNA, correctProt, outputCDNA, outputProt string
	var errorCount int
	for val := range answer {
		annotation = strings.Split(val.Info, ";")
		words = strings.Split(annotation[0], "|")
		correctCDNA = words[0]
		correctProt = words[1]
		newWords = strings.Split(annotation[1], "|")
		newerWords = strings.Split(newWords[2], ":")
		outputCDNA = newerWords[1]
		outputProt = newWords[3]
		if (outputCDNA == correctCDNA && outputProt == correctProt) ||
			strings.HasPrefix(correctCDNA, "c.-") || strings.HasPrefix(correctCDNA, "c.*") {
			continue
		}

		errorCount++
		//fmt.Printf("\nWARNING: ANNOTATION MISMATCH\n")
		//fmt.Printf("EXPECTED: %s|%s\n", correctCDNA, correctProt)
		//fmt.Printf("RECEIVED: %s\n", annotation)

	}
	if errorCount > 7 { // from known issue
		t.Errorf("ERROR: %d variants were misannotated", errorCount)
	}
}

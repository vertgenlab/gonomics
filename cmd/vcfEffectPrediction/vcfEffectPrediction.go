// Command Group: "Variant Calling & Annotation"

package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/gtf"
	"github.com/vertgenlab/gonomics/interval"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
	"strings"
	"sync"
)

type Settings struct {
	Vcf            string
	Gtf            string
	Fasta          string
	Threads        int
	AllTranscripts bool
	OutFile        string
}

func usage() {
	fmt.Print(
		"vcfEffectPrediction - Annotate Vcf records with cDNA and protein effect predictions.\n" +
			"\tAnnotations are added to the INFO field of the Vcf with the following format: \"GoEP= g.XXX | Gene | TranscriptId:c.XXX | p.XXX | VariantType\"\n" +
			"\tVariantTypes include: Silent, Missense, Nonsense, Frameshift, Intergenic, Intronic, Splice (1-2 away from intron-exon boundary), FarSplice (3-10 away from intron-exon boundary)\n" +
			"Usage:\n" +
			" vcfEffectPrediction [options] -fasta ref.fa -gtf ref.gtf input.vcf output.vcf \n\n" +
			"options:\n")
	flag.PrintDefaults()
}

func AppendAnnotationHeader(header vcf.Header) vcf.Header {
	var columnIDs string
	if strings.HasPrefix(header.Text[len(header.Text)-1], "#CHROM\t") {
		columnIDs = header.Text[len(header.Text)-1]
		header.Text = header.Text[:len(header.Text)-1]
	}

	header.Text = append(header.Text, "##GoEffectPrediction Version=1.0\n")
	header.Text = append(header.Text, "##INFO=<ID=GoEP,Number=.,Type=String,Description=\"Functional annotations: HGVS.g | Gene | TranscriptId : HGVS.c | HGVS.p | VariantType\">\n")

	if columnIDs != "" {
		header.Text = append(header.Text, columnIDs)
	}

	return header
}

func vcfEffectPrediction(settings *Settings) (<-chan vcf.Vcf, vcf.Header) {
	f := fasta.Read(settings.Fasta)
	fasta.AllToUpper(f)
	fastaRecords := fasta.ToMap(f)
	gtfRecords := gtf.Read(settings.Gtf)
	tree := gtf.GenesToIntervalTree(gtfRecords)

	vcfChan, vcfHeader := vcf.GoReadToChan(settings.Vcf)
	answer := make(chan vcf.Vcf, 1000)
	var wg sync.WaitGroup

	for i := 0; i < settings.Threads; i++ {
		wg.Add(1)
		go annotationWorker(&wg, tree, fastaRecords, vcfChan, answer, settings.AllTranscripts)
	}

	go func() {
		wg.Wait()
		close(answer)
	}()

	vcfHeader = AppendAnnotationHeader(vcfHeader)
	return answer, vcfHeader
}

func annotationWorker(wg *sync.WaitGroup, tree map[string]*interval.IntervalNode, fasta map[string][]dna.Base, vcfChan <-chan vcf.Vcf, answer chan<- vcf.Vcf, allTranscripts bool) {
	var annotation string
	for vcfRecord := range vcfChan {
		variant, _ := gtf.VcfToVariant(vcfRecord, tree, fasta, allTranscripts) //TODO: make gtf.vcfEffectPrediciton a public struct and declare outside of loop
		annotation = gtf.VariantToAnnotation(variant, fasta)
		vcfRecord.Info = vcfRecord.Info + ";" + annotation
		answer <- vcfRecord
	}
	wg.Done()
}

func main() {
	var fasta *string = flag.String("fasta", "", "Fasta file used to generate Vcf file.")
	var gtf *string = flag.String("gtf", "", "Gtf/Gff file with coordinates corresponding to the input Fasta file. \nNote: Row order of Gtf file must match UCSC Genome Browser style.")
	var threads *int = flag.Int("threads", 1, "Number of threads to use. \nNote: if >1 thread is used, the order of output Vcf may not match the input Vcf.")
	var allTranscripts *bool = flag.Bool("allTranscripts", false, "Generate annotation for each transcript isoform. Default only annotates for the canonical transcript.\n"+
		"Note: Canonical transcript will always be reported first, predictions for additional transcripts are appended to the \"GoEP\" Info field with repeating \"| HGVS.c | HGVS.p | VariantType\"")
	flag.Parse()

	var expectedNumArgs int = 2
	flag.Usage = usage

	if *fasta == "" || *gtf == "" {
		usage()
		log.Fatalf("ERROR: Must input a fasta and a gtf")
	}

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n\n", expectedNumArgs, len(flag.Args()))
	} else {
		vcfFilename, outFile := flag.Arg(0), flag.Arg(1)
		options := &Settings{
			Vcf:            vcfFilename,
			Gtf:            *gtf,
			Fasta:          *fasta,
			Threads:        *threads,
			AllTranscripts: *allTranscripts,
			OutFile:        outFile,
		}
		outfile := fileio.EasyCreate(options.OutFile)
		answer, header := vcfEffectPrediction(options)
		vcf.NewWriteHeader(outfile, header)
		for val := range answer {
			vcf.WriteVcf(outfile, val)
		}
	}
}

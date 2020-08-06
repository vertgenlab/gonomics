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
		"vcfAnnotate - Annotate Vcf records with cDNA and protein effect predictions.\n" +
			"\tAnnotations are added to the INFO field of the Vcf with the following format: \"GoEP= g.XXX | VariantType | Gene | TranscriptId:c.XXX | p.XXX\"\n" +
			"\tVariantTypes include: Silent, Missense, Nonsense, Frameshift, Intergenic, Intronic, Splice (1-2 away from intron-exon boundary), FarSplice (3-10 away from intron-exon boundary)\n" +
			"Usage:\n" +
			" vcfAnnotate [options]\n" +
			"\t-fasta ref.fa \n" +
			"\t-gtf ref.gtf \n" +
			"\tinput.vcf \n" +
			"\toutput.vcf \n\n" +
			"options:\n")
	flag.PrintDefaults()
}

func vcfAnnotate(settings *Settings) (<-chan *vcf.Vcf, *vcf.VcfHeader) {
	f := fasta.Read(settings.Fasta)
	fasta.AllToUpper(f)
	fastaRecords := fasta.FastaMap(f)
	gtfRecords := gtf.Read(settings.Gtf)
	tree := gtf.GenesToIntervalTree(gtfRecords)

	vcfChan, vcfHeader := vcf.GoReadToChan(settings.Vcf)
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

	//gtf.AppendAnnotationHeader(vcfHeader) //TODO: Make this function to add header lines for addition data in INFO
	return answer, vcfHeader
}

func annotationWorker(wg *sync.WaitGroup, tree map[string]*interval.IntervalNode, fasta map[string][]dna.Base, vcfChan <-chan *vcf.Vcf, answer chan<- *vcf.Vcf, allTranscripts bool) {
	var annotation string
	for vcfRecord := range vcfChan { //TODO Add allTranscripts option to gtf.VcfToVariant and gtf.GenesToIntervalTree
		variant, _ := gtf.VcfToVariant(vcfRecord, tree, fasta) //TODO: make gtf.vcfEffectPrediciton a public struct and declare outside of loop
		annotation = gtf.VariantToAnnotation(variant, fasta)
		vcfRecord.Info = vcfRecord.Info + ";" + annotation
		answer <- vcfRecord
	}
	wg.Done()
}

func main() {
	var fasta *string = flag.String("fasta", "", "Fasta file used to generate Vcf file.")
	var gtf *string = flag.String("gtf", "", "Gtf/Gff file with coordinates corresponding to the input Fasta file. \nNote: Row order of Gtf file must match UCSC Genome Browser style.")
	var threads *int = flag.Int("threads", 1, "Number of threads to use. \nNOTE: if >1 thread is used, the order of output Vcf may not match the input Vcf.")
	var allTranscripts *bool = flag.Bool("allTranscripts", false, "Generate annotation for each transcript isoform. Default only annotates for the canonical transcript. `WORK IN PROGRESS`") //TODO
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
		answer, header := vcfAnnotate(options)
		vcf.NewWriteHeader(outfile, header)
		for val := range answer {
			vcf.WriteVcf(outfile, val)
		}
	}
}

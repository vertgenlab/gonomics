package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/interval"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
)

func haplotypeGenerator(genomeFile string, snpFile string, regionFile string, outdir string) {
	variants, header := vcf.Read(snpFile)
	variantIntervals := make([]interval.Interval, len(variants))
	for i := range variants {
		variantIntervals[i] = variants[i]
	}

	tree := interval.BuildTree(variantIntervals)

	regions := bed.Read(regionFile)
	genome := fasta.Read(genomeFile)
	genomeMap := helperFastaIndex(genome)

	samplesNames := vcf.SampleNamesInOrder(header)

	for i := range regions {
		overlappingVariants := interval.Query(tree, regions[i], "any")
		currIndex := genomeMap[regions[i].Chrom]
		refHaplotype := fasta.Extract(genome[currIndex], regions[i].ChromStart, regions[i].ChromEnd, regions[i].Chrom)

		outputFilename := fmt.Sprintf("%s/%s.%v.%v.fa", outdir, regions[i].Chrom, regions[i].ChromStart, regions[i].ChromEnd)
		out := fileio.EasyCreate(outputFilename)
		sampleHaplotypes := make([]fasta.Fasta, len(header.Samples)*2)

		for j := range sampleHaplotypes {
			sampleHaplotypes[j] = fasta.Copy(refHaplotype)
			currSampleIndex := j / 2
			if j%2 == 0 {
				sampleHaplotypes[j].Name = fmt.Sprintf("%s_A", samplesNames[currSampleIndex])
			} else {
				sampleHaplotypes[j].Name = fmt.Sprintf("%s_B", samplesNames[currSampleIndex])
			}

			for k := range overlappingVariants {
				// For each variant, update corresponding position in each sample haplotype (position = position - ChromStart - 1)
				// sampleHaplotypes[j]
				fmt.Println(k)

			}

			fasta.WriteFasta(out, sampleHaplotypes[j], 50)
		}

		err := out.Close()
		exception.PanicOnErr(err)

	}

}

// Helper function returns a map connecting chromosome names to their index in a fasta slice
func helperFastaIndex(genome []fasta.Fasta) map[string]int {
	var answer = make(map[string]int)
	for i := range genome {
		answer[genome[i].Name] = i
	}
	return answer
}

func usage() {
	fmt.Print(
		"haplotypeGenerator - Generate unique haplotypes for provided regions from genetic variation data.\n" +
			"Ignores structural variants, considers only substitutions\n" +
			"Usage:\n" +
			"genomeFile.fa snpFile.vcf regionFile.bed outDir\n" +
			"options:\n")
	flag.PrintDefaults()
}

func main() {
	var expectedNumArgs int = 4

	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if len(flag.Args()) != expectedNumArgs {
		flag.Usage()
		log.Fatalf("Error: expecting %d arguments, but got %d\n",
			expectedNumArgs, len(flag.Args()))
	}

	genomeFile := flag.Arg(0)
	snpFile := flag.Arg(1)
	regionFile := flag.Arg(2)
	outdir := flag.Arg(3)

	haplotypeGenerator(genomeFile, snpFile, regionFile, outdir)
}

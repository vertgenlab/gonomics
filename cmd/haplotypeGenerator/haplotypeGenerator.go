package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/interval"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
)

func haplotypeGenerator(genomeFile string, snpFile string, regionFile string, outdir string) {
	var currState int
	variants, header := vcf.Read(snpFile)
	variantIntervals := make([]interval.Interval, len(variants))

	// fmt.Println("started ok")

	for i := range variants {
		variantIntervals[i] = variants[i]
	}

	tree := interval.BuildTree(variantIntervals)

	regions := bed.Read(regionFile)
	genome := fasta.Read(genomeFile)
	genomeMap := helperFastaIndex(genome)

	samplesNames := vcf.SampleNamesInOrder(header)

	for i := range regions {

		// fmt.Printf("working on region index %v.\n", i)

		overlappingVariants := interval.Query(tree, regions[i], "any")
		currIndex := genomeMap[regions[i].Chrom]
		refHaplotype := fasta.Extract(genome[currIndex], regions[i].ChromStart, regions[i].ChromEnd, regions[i].Chrom)

		outputFilename := fmt.Sprintf("%s/%s.%v.%v.fa", outdir, regions[i].Chrom, regions[i].ChromStart, regions[i].ChromEnd)

		// fmt.Printf("output filename: %s\n", outputFilename)

		out := fileio.EasyCreate(outputFilename)
		sampleHaplotypes := make([]fasta.Fasta, len(header.Samples)*2)

		// fmt.Printf("last line of header: \n%s", header.Text[len(header.Text)-1])

		// fmt.Printf("length of sample haplotypes: %v.\n", len(sampleHaplotypes))
		// fmt.Printf("length of header.Samples: %v.\n", len(header.Samples))

		for j := range sampleHaplotypes {
			sampleHaplotypes[j] = fasta.Copy(refHaplotype)
			currSampleIndex := j / 2
			if j%2 == 0 {
				sampleHaplotypes[j].Name = fmt.Sprintf("%s_A", samplesNames[currSampleIndex])
			} else {
				sampleHaplotypes[j].Name = fmt.Sprintf("%s_B", samplesNames[currSampleIndex])
			}

			// fmt.Printf("working on haplotype: \n%s", sampleHaplotypes[j].Name)

			for k := range overlappingVariants {
				if j%2 == 0 {
					currState = int(overlappingVariants[k].(vcf.Vcf).Samples[currSampleIndex].Alleles[0])
					if currState > 0 {
						sampleHaplotypes[j].Seq[overlappingVariants[k].(vcf.Vcf).Pos-regions[i].ChromStart-1] = dna.StringToBase(overlappingVariants[k].(vcf.Vcf).Alt[currState-1])
					}
				} else {
					currState = int(overlappingVariants[k].(vcf.Vcf).Samples[currSampleIndex].Alleles[1])
					if currState > 0 {
						sampleHaplotypes[j].Seq[overlappingVariants[k].(vcf.Vcf).Pos-regions[i].ChromStart-1] = dna.StringToBase(overlappingVariants[k].(vcf.Vcf).Alt[currState-1])
					}
				}

			}
			// fmt.Println("writing haplotype")
			// fmt.Println(sampleHaplotypes[j])

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

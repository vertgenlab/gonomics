package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/alleles"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/simpleGraph"
	"github.com/vertgenlab/gonomics/vcf"
	"log"
	"path/filepath"
	"strings"
	"time"
)

func usage() {
	fmt.Print(
		"callVariants - Calls variants from input sam or giraf based on a linear or graph reference.\n" +
			"Usage:\n" +
			" callVariants [options] \n" +
			"\t-i experimental.sam \n" +
			"\t-n normal.sam \n" +
			"\t-lr reference.fasta \n" +
			"\t-o output.vcf \n\n" +
			"options:\n")
	flag.PrintDefaults()
}

func callVariants(linearRef string, graphRef string, expSamples string, normSamples string, outFile string, afThreshold float64, sigThreshold float64, minMapQ int64, memBufferSize int, minCov int) {
	var ref interface{}
	output := fileio.MustCreate(outFile)
	defer output.Close()

	if linearRef != "" {
		ref = fasta.Read(linearRef)
	} else if graphRef != "" {
		ref = simpleGraph.Read(graphRef)
	}
	alleleStream, normalIDs := startAlleleStreams(ref, expSamples, normSamples, minMapQ, memBufferSize)
	answer := alleles.FindNewVariation(alleleStream, normalIDs, afThreshold, sigThreshold, minCov)

	lastProgressReport := time.Now()

	//vcf.WriteHeader(output)
	for vcfRecord := range answer {
		if time.Since(lastProgressReport).Seconds() > 10 {
			lastProgressReport = time.Now()
			log.Printf("Current Position: %s\t%d", vcfRecord.Chr, vcfRecord.Pos)
		}
		vcf.WriteVcf(output, vcfRecord)
	}
}

// fills alleleChans and normalIDs with the appropriate data
// file can be .sam, .giraf, or .txt as a list of .sams or .girafs
func addChans(ref interface{}, file string, isNormal bool, alleleChans *[]<-chan *alleles.Allele, normalIDs map[string]bool, samFilesPresent *bool, girafFilesPresent *bool, minMapQ int64) {
	switch filepath.Ext(file) {
	case ".giraf":
		//TODO: Fix giraf to alleles function
		log.Fatalln("ERROR: giraf files are currently not supported")
		//*alleleChans = append(*alleleChans, alleles.(file))
		*girafFilesPresent = true
		if isNormal == true {
			normalIDs[file] = true
		} else {
			normalIDs[file] = false
		}
		log.Println("Started Allele Stream for", file)
	case ".sam":
		*alleleChans = append(*alleleChans, alleles.GoCountSamAlleles(file, ref.([]fasta.Fasta), minMapQ))
		*samFilesPresent = true
		if isNormal == true {
			normalIDs[file] = true
		} else {
			normalIDs[file] = false
		}
		log.Println("Started Allele Stream for", file)
	case ".txt":
		reader := fileio.EasyOpen(file)
		defer reader.Close()
		for line, done := fileio.EasyNextLine(reader); !done; line, done = fileio.EasyNextLine(reader) {
			if line == file {
				log.Fatalln("ERROR: Infinite recursion detected: Cannot call", line, "within", file)
			}
			if strings.HasPrefix(line, "#") {
				continue
			} else {
				addChans(ref, line, isNormal, alleleChans, normalIDs, samFilesPresent, girafFilesPresent, minMapQ)
			}
		}

	default:
		log.Println("ERROR: Did not recognize extension", filepath.Ext(file), "for input:", file)
	}
}

// Starts and syncs the allele streams for both experimental and normal files,
// returns channel to synced alleles and a map[filename]bool (true = normal)
func startAlleleStreams(ref interface{}, experimental string, normal string, minMapQ int64, memBufferSize int) (<-chan []*alleles.Allele, map[string]bool) {
	var alleleChans []<-chan *alleles.Allele
	normalIDs := make(map[string]bool) // map[filename] (true = normal)
	var samFilesPresent, girafFilesPresent bool

	addChans(ref, experimental, false, &alleleChans, normalIDs, &samFilesPresent, &girafFilesPresent, minMapQ)
	if normal != "" {
		addChans(ref, normal, true, &alleleChans, normalIDs, &samFilesPresent, &girafFilesPresent, minMapQ)
	}

	if samFilesPresent && girafFilesPresent {
		log.Fatalln("ERROR: Input directories contain both giraf and sam files")
	}
	if len(alleleChans) < 2 {
		log.Fatalln("ERROR: Must input at least two samples (between experimental and normal) to facilitate comparisons. Only", len(alleleChans), "sample was submitted")
	}

	syncedAllelesChan := alleles.SyncAlleleStreams(ref, memBufferSize, alleleChans...)

	return syncedAllelesChan, normalIDs
}

func main() {
	var outFile *string = flag.String("out", "", "Write output to a file [.vcf].")
	var sigThreshold *float64 = flag.Float64("p", 0.05, "Do not output variants with p value greater than this value.")
	var afThreshold *float64 = flag.Float64("af", 0.01, "Variants with allele frequency less than this value will be treated as p = 1.")
	var linearReference *string = flag.String("lr", "", "Linear reference used for alignment [.fasta].")
	var graphReference *string = flag.String("gr", "", "Graph reference used for alignment [.gg].")
	var experimentalSamples *string = flag.String("i", "", "Input experimental sample(s) [.sam, .giraf, .txt]. Can be a file or a txt file with a list (must have .txt extension) of sample paths.")
	var normalSamples *string = flag.String("n", "", "Input normal sample(s) [.sam, .giraf, .txt]. Can be a file or a txt file with a list (must have .txt extension) of sample paths. If no normal samples are given, each experimental sample will me measured against the other experimental samples.")
	var minMapQ *int64 = flag.Int64("minMapQ", 20, "Exclude all reads with mapping quality less than this value")
	var memBufferSize *int = flag.Int("memBuffer", 100, "Maximum number of allele records to store in memory at once")
	var minCov *int = flag.Int("minCoverage", 10, "Minimum number of covering reads to be considered a valid variant")
	flag.Usage = usage
	log.SetFlags(log.Ldate | log.Ltime | log.Lshortfile)
	flag.Parse()

	if *linearReference != "" && *graphReference != "" {
		log.Fatalln("ERROR: Cannot input both linear and graph references")
	}

	if *experimentalSamples == "" || *outFile == "" || (*linearReference == "" && *graphReference == "") {
		flag.Usage()
		log.Fatalf("ERROR: Must include parameters for -i, -out, (-lr or -gr)")
	}
	flag.Parse()

	callVariants(*linearReference, *graphReference, *experimentalSamples, *normalSamples, *outFile, *afThreshold, *sigThreshold, *minMapQ, *memBufferSize, *minCov)
}

package alleles

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/vcf"
	"io/ioutil"
	"log"
	"os"
	"path/filepath"
)

func CallVariants(ref interface{}, experimental string, normal string, afThreshold float64, sigThreshold float64, minMapQ int64, memBufferSize int) <-chan *vcf.Vcf {
	answer := make(chan *vcf.Vcf)
	alleleStream, normalIDs := startAlleleStreams(ref, experimental, normal, minMapQ, memBufferSize)
	go scoreAlleles(answer, alleleStream, normalIDs, afThreshold, sigThreshold)
	return answer
}

// Returns true if input path is a directory
func isDirectory(path string) bool {
	file, _ := os.Stat(path)

	switch fileInfo := file.Mode(); {
	case fileInfo.IsDir():
		return true
	case fileInfo.IsRegular():
		return false
	default:
		log.Fatalln("ERROR: Could not determine if input was file or dir")
	}
	return false
}

// Starts and syncs the allele streams for both experimental and normal files,
// returns channel to synced alleles and a map[filename]bool (true = normal)
func startAlleleStreams(ref interface{}, experimental string, normal string, minMapQ int64, memBufferSize int) (<-chan []*Allele, map[string]bool) {
	var alleleChans []<-chan *Allele
	normalIDs := make(map[string]bool) // map[filename] (true = normal)
	var samFilesPresent, girafFilesPresent bool

	if isDirectory(experimental) {
		files, _ := ioutil.ReadDir(experimental)
		for _, file := range files {
			path := fmt.Sprintf("%s/%s", experimental, file.Name())
			switch filepath.Ext(path) {
			case ".giraf":
				//TODO: Make giraf to alleles function
				alleleChans = append(alleleChans, GirafToAlleles(path))
				normalIDs[path] = false
				girafFilesPresent = true
				log.Println("Started Allele Stream for", path)
			case ".sam":
				alleleChans = append(alleleChans, SamToAlleles(path, ref, minMapQ))
				normalIDs[path] = false
				samFilesPresent = true
				log.Println("Started Allele Stream for", path)
			}
		}
	} else { // if input is only a single file
		switch filepath.Ext(experimental) {
		case ".giraf":
			alleleChans = append(alleleChans, GirafToAlleles(experimental))
			normalIDs[experimental] = false
			girafFilesPresent = true
			log.Println("Started Allele Stream for", experimental)
		case ".sam":
			alleleChans = append(alleleChans, SamToAlleles(experimental, ref, minMapQ))
			normalIDs[experimental] = false
			samFilesPresent = true
			log.Println("Started Allele Stream for", experimental)
		}
	}

	if normal != "" && isDirectory(normal) {
		files, _ := ioutil.ReadDir(normal)
		for _, file := range files {
			path := fmt.Sprintf("%s/%s", normal, file.Name())
			switch filepath.Ext(path) {
			case ".giraf":
				//TODO: Make giraf to alleles function
				alleleChans = append(alleleChans, GirafToAlleles(path))
				normalIDs[path] = true
				girafFilesPresent = true
				log.Println("Started Allele Stream for", path)
			case ".sam":
				alleleChans = append(alleleChans, SamToAlleles(path, ref, minMapQ))
				normalIDs[path] = true
				samFilesPresent = true
				log.Println("Started Allele Stream for", path)
			}
		}
	} else { // if input is only a single file
		switch filepath.Ext(normal) {
		case ".giraf":
			alleleChans = append(alleleChans, GirafToAlleles(normal))
			normalIDs[normal] = true
			girafFilesPresent = true
			log.Println("Started Allele Stream for", normal)
		case ".sam":
			alleleChans = append(alleleChans, SamToAlleles(normal, ref, minMapQ))
			normalIDs[normal] = true
			samFilesPresent = true
			log.Println("Started Allele Stream for", normal)
		}
	}

	if samFilesPresent && girafFilesPresent {
		log.Fatalln("ERROR: Input directories contain both giraf and sam files")
	}

	syncedAllelesChan := SyncAlleleStreams(ref, memBufferSize, alleleChans...)

	return syncedAllelesChan, normalIDs
}

// Designed to be run as a goroutine that accepts alleles from the alleleStream channel, computes the p value and makes a VCF,
// then sends the vcf record on the answer channel
func scoreAlleles(answer chan<- *vcf.Vcf, alleleStream <-chan []*Allele, normalIDs map[string]bool, afThreshold float64, sigThreshold float64) {

	for alleles := range alleleStream {
		if len(alleles) == 1 {continue} // Can only evaluate when at least two samples are present
		bkgd, normalPresent := calcBackground(alleles, normalIDs)

		// Begin gathering parameters for Fishers Exact Test done in the numbers package
		// test is for the matrix:
		// [a b]
		// [c d]
		// a = Samples Ref Allele Count
		// b = Background Ref Allele Count - Samples Ref Allele Count
		// c = Samples Alt Allele Count
		// d = Background Alt Allele Count - Samples Alt Allele Count
		fmt.Sprintln(bkgd)

		if normalPresent {
			// If the bkgd was based on normal samples then we do not need to
			// subtract out the values of the experimental sample

		} else {
			// If the bkgd was based on the experimental samples then we DO need
			// to subtract out the values of the experimental sample

		}

		answer <- &vcf.Vcf{Notes: fmt.Sprintln(alleles)}
	}

	close(answer)
}

func calcBackground(alleles []*Allele, normalIDs map[string]bool) (*Allele, bool) {
	var answer *Allele = &Allele{Count: &AlleleCount{alleles[0].Count.Ref, 0, 0, 0, 0, 0, 0, 0, 0, 0, make([]Indel, 0)}}
	var normalPresent bool

	for i := 0; i < len(alleles); i++ {
		if normalIDs[alleles[i].Sample] {
			addAlleles(answer, alleles[i])
			normalPresent = true
		}
	}

	if normalPresent {
		return answer, normalPresent
	} else {
		for k := 0; k < len(alleles); k++ {
			addAlleles(answer, alleles[k])
		}
	}
	return answer, normalPresent
}

func addAlleles(sum *Allele, toBeAdded *Allele) {
	sum.Count.BaseAF += toBeAdded.Count.BaseAF
	sum.Count.BaseAR += toBeAdded.Count.BaseAR
	sum.Count.BaseCF += toBeAdded.Count.BaseCF
	sum.Count.BaseCR += toBeAdded.Count.BaseCR
	sum.Count.BaseGF += toBeAdded.Count.BaseGF
	sum.Count.BaseGR += toBeAdded.Count.BaseGR
	sum.Count.BaseTF += toBeAdded.Count.BaseTF
	sum.Count.BaseTR += toBeAdded.Count.BaseTR

	for j := 0; j < len(toBeAdded.Count.Indel); j++ {
		indel := findMatchingIndel(&toBeAdded.Count.Indel[j], sum.Count.Indel)
		if indel != nil {
			indel.CountF += toBeAdded.Count.Indel[j].CountF
			indel.CountR += toBeAdded.Count.Indel[j].CountR
		} else {
			sum.Count.Indel = append(sum.Count.Indel, toBeAdded.Count.Indel[j])
		}
	}
}

func findMatchingIndel(queryIndel *Indel, subjectSlice []Indel) *Indel {
	var i int
	for i = 0; i < len(subjectSlice); i++ {
		if dna.CompareSeqsIgnoreCase(queryIndel.Alt, subjectSlice[i].Alt) == 0 &&
			dna.CompareSeqsIgnoreCase(queryIndel.Ref, subjectSlice[i].Ref) == 0 {
			return &subjectSlice[i]
		}
	}
	return nil
}

func alleleToVcf(allele *Allele, p float64) *vcf.Vcf {
	var answer *vcf.Vcf

	return answer
}

func alleleFishersExact(a int32, b int32, c int32, d int32, afThreshold float64) float64 {
	var p float64

	switch {
	// If alternate allele is zero then there is no variant and score is 1
	case c == 0:
		p = 1

	// If a = b and c = d then it is testing itself and should return 1
	case a == b && c == d:
		p = 1

	// If the allele frequency of d > c then p is 1
	case float64(c)/float64(c+a) < float64(d)/float64(d+b):
		p = 1

	// If the allele frequency is less than the threshold then p is noted as 1 so as to be excluded
	case float64(c)/float64(c+a) < afThreshold:
		p = 1

	// If no exclusion conditions are met, then calculate p value
	default:
		p = numbers.FisherExact(int(a), int(b), int(c), int(d), true)
	}
	return p
}
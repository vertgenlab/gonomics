package alleles

import (
	"fmt"
	"github.com/vertgenlab/gonomics/vcf"
	"io/ioutil"
	"log"
	"os"
	"path/filepath"
)

func CallVariants(ref interface{}, experimental string, normal string, afThreshold float64, sigThreshold float64, minMapQ int64, memBufferSize int) chan *vcf.Vcf {
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
				log.Println("Started Allele Stream for ", path)
			case ".sam":
				alleleChans = append(alleleChans, SamToAlleles(path, ref, minMapQ))
				normalIDs[path] = false
				samFilesPresent = true
				log.Println("Started Allele Stream for ", path)
			}
		}
	} else { // if input is only a single file
		switch filepath.Ext(experimental) {
		case ".giraf":
			alleleChans = append(alleleChans, GirafToAlleles(experimental))
			normalIDs[experimental] = false
			girafFilesPresent = true
			log.Println("Started Allele Stream for ", experimental)
		case ".sam":
			alleleChans = append(alleleChans, SamToAlleles(experimental, ref, minMapQ))
			normalIDs[experimental] = false
			samFilesPresent = true
			log.Println("Started Allele Stream for ", experimental)
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
				log.Println("Started Allele Stream for ", path)
			case ".sam":
				alleleChans = append(alleleChans, SamToAlleles(path, ref, minMapQ))
				normalIDs[path] = true
				samFilesPresent = true
				log.Println("Started Allele Stream for ", path)
			}
		}
	} else { // if input is only a single file
		switch filepath.Ext(normal) {
		case ".giraf":
			alleleChans = append(alleleChans, GirafToAlleles(normal))
			normalIDs[normal] = true
			girafFilesPresent = true
			log.Println("Started Allele Stream for ", normal)
		case ".sam":
			alleleChans = append(alleleChans, SamToAlleles(normal, ref, minMapQ))
			normalIDs[normal] = true
			samFilesPresent = true
			log.Println("Started Allele Stream for ", normal)
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
	for test := range alleleStream {
		answer <- &vcf.Vcf{Notes: fmt.Sprintln(test)}
	}
	close(answer)
}
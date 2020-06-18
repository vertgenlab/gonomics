package alleles

import (
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"path/filepath"
	"strings"
	"testing"
)

// TODO: build more robust test
// TODO: build test for graph genome
func TestFindNewVariation(t *testing.T) {
	ref := fasta.Read("testdata/human_chrM.fasta")
	experimental := "testdata/human_chrM.sam"
	normal := "testdata/chrM_head.sam"
	alleleStreams, normalIDs := startAlleleStreams(ref, experimental, normal, 20, 100)
	FindNewVariation(alleleStreams, normalIDs, 0, 1)
}

// fills alleleChans and normalIDs with the appropriate data
// file can be .sam, .giraf, or .txt as a list of .sams or .girafs
func addChans(ref interface{}, file string, isNormal bool, alleleChans *[]<-chan *Allele, normalIDs map[string]bool, samFilesPresent *bool, girafFilesPresent *bool, minMapQ int64) {
	switch filepath.Ext(file) {
	case ".giraf":
		//TODO: Make giraf to alleles function
		*alleleChans = append(*alleleChans, GirafToAlleles(file))
		*girafFilesPresent = true
		if isNormal == true {
			normalIDs[file] = true
		} else {
			normalIDs[file] = false
		}
		//log.Println("Started Allele Stream for", file)
	case ".sam":
		*alleleChans = append(*alleleChans, SamToAlleles(file, ref, minMapQ))
		*samFilesPresent = true
		if isNormal == true {
			normalIDs[file] = true
		} else {
			normalIDs[file] = false
		}
		//log.Println("Started Allele Stream for", file)
	case ".txt":
		reader := fileio.EasyOpen(file)
		defer reader.Close()
		for line, done := fileio.EasyNextLine(reader); !done; line, done = fileio.EasyNextLine(reader) {
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
func startAlleleStreams(ref interface{}, experimental string, normal string, minMapQ int64, memBufferSize int) (<-chan []*Allele, map[string]bool) {
	var alleleChans []<-chan *Allele
	normalIDs := make(map[string]bool) // map[filename] (true = normal)
	var samFilesPresent, girafFilesPresent bool

	addChans(ref, experimental, false, &alleleChans, normalIDs, &samFilesPresent, &girafFilesPresent, minMapQ)
	addChans(ref, normal, true, &alleleChans, normalIDs, &samFilesPresent, &girafFilesPresent, minMapQ)

	if samFilesPresent && girafFilesPresent {
		log.Fatalln("ERROR: Input directories contain both giraf and sam files")
	}
	if len(alleleChans) < 2 {
		log.Fatalln("ERROR: Must input at least two samples (between experimental and normal) to facilitate comparisons. Only", len(alleleChans), "sample was submitted")
	}

	syncedAllelesChan := SyncAlleleStreams(ref, memBufferSize, alleleChans...)

	return syncedAllelesChan, normalIDs
}

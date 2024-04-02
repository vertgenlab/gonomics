package sam

<<<<<<< HEAD
import "log"

type SamPE struct {
	R1 Sam
	R2 Sam
}

// IsRead1 takes a sam entry and return true is it is read 1 and false if it is read 2
func IsRead1(a Sam) bool {
	var r1 bool
	if a.Flag&0x40 == 0x40 {
		r1 = true
	} else if a.Flag&0x80 == 0x80 {
		r1 = false
	} else {
		log.Fatalln("ERROR in IsRead1: sam entry didn't have either the read 1 or 2 flag set. Check that your sam file is paired end.", a)
	}
	return r1
}

func SamToPeSamCheckFlag(a, b Sam) (SamPE, bool) {
	var readPair SamPE
	var poss bool
	if a.Flag&0x40 == 0x40 && b.Flag&0x80 == 0x80 {
		readPair = SamPE{
			R1: a,
			R2: b,
		}
		poss = true
	} else if a.Flag&0x80 == 0x80 && b.Flag&0x40 == 0x40 {
		readPair = SamPE{
			R1: b,
			R2: a,
		}
		poss = true
	} else {
		readPair = SamPE{}
		poss = false
		//log.Fatalf("Error in SamToPeCheckFlag: the input Sam entries did not have the proper flags indicating read pairing.\n"+
		//"R1 name: %s\t R1 flag: %d\n"+
		//"R2 name: %s\t R2 flag: %d", a.QName, a.Flag, b.QName, b.Flag)
	}
	return readPair, poss
=======
import (
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"io"
	"strings"
)

// SamPE is struct that contains 2 paired end Sam entries as well as a slice of Sam that represents all supplementary
// alignments for the read-pair
type SamPE struct { //TODO: rename from A/B to something more informative
	R1    Sam
	R2    Sam
	SupR1 []Sam
	SupR2 []Sam
}

// GoReadSamPeToChan takes a read-name sorted sam file and returns a chanel of SamPE and a Header
func GoReadSamPeToChan(filename string) (<-chan SamPE, Header) {

	//start by reading the file into a chanel of sam
	data := make(chan Sam, 1000)
	header := make(chan Header)
	peChan := make(chan SamPE, 1000)

	if strings.HasSuffix(filename, ".bam") {
		go readBamToChan(filename, data, header)
	} else {
		go readSamToChan(filename, data, header)
	}

	go combineEnds(data, peChan)

	return peChan, <-header
}

// combineEnds is a helper function for GoReadSamPeToChan that takes a channel of type Sam and parses it into a SamPE channel
func combineEnds(data <-chan Sam, peChan chan<- SamPE) {
	var full bool = true
	var tmpSlice []Sam
	var tmp Sam

	for full {
		if len(tmpSlice) == 0 {
			tmp, full = <-data
			tmpSlice = append(tmpSlice, tmp)
			continue
		}
		tmp, full = <-data
		if tmpSlice[0].QName == tmp.QName && full {
			tmpSlice = append(tmpSlice, tmp)
			continue
		} else {
			peChan <- parseReads(tmpSlice)
			tmpSlice = tmpSlice[0:0]
			tmpSlice = append(tmpSlice, tmp)
		}
	}
	close(peChan)
}

// parseReads is a helper function for GoReadSamPeToChan that parses a slice of Sam that all have the same read name into a SamPE struct
func parseReads(tmpSlice []Sam) SamPE {
	var pe SamPE
	for i := range tmpSlice {
		if tmpSlice[i].Flag&2048 == 2048 { // Flag 2048 == supplementary reads
			if tmpSlice[i].Flag&64 == 64 { // Flag 64 == Read 1
				pe.SupR1 = append(pe.SupR1, tmpSlice[i])
			} else if tmpSlice[i].Flag&128 == 128 { // Flag 128 == Read 2
				pe.SupR2 = append(pe.SupR2, tmpSlice[i])
			}
		} else if tmpSlice[i].Flag&64 == 64 { // is Read 1
			pe.R1 = tmpSlice[i]
		} else if tmpSlice[i].Flag&128 == 128 { // is Read 2
			pe.R2 = tmpSlice[i]
		}
	}
	return pe
}

// WriteSamPeToFileHandle writes a both reads and any supplementary alignments from a SamPE struct to a fileio EasyWriter
// TODO: Add bam writing functionality
func WriteSamPeToFileHandle(file io.Writer, pe SamPE) {
	_, err := fmt.Fprintln(file, ToString(pe.R1))
	exception.PanicOnErr(err)
	_, err = fmt.Fprintln(file, ToString(pe.R2))
	exception.PanicOnErr(err)
	if len(pe.SupR1) != 0 {
		for i := range pe.SupR1 {
			_, err = fmt.Fprintln(file, ToString(pe.SupR1[i]))
			exception.PanicOnErr(err)
		}
	}
	if len(pe.SupR2) != 0 {
		for i := range pe.SupR2 {
			_, err = fmt.Fprintln(file, ToString(pe.SupR2[i]))
			exception.PanicOnErr(err)
		}
	}
>>>>>>> main
}

package sam

import (
	"fmt"
	"io"
	"strings"

	"github.com/vertgenlab/gonomics/exception"
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
	_, err := fmt.Fprintln(file, pe.R1.String())
	exception.PanicOnErr(err)
	_, err = fmt.Fprintln(file, pe.R2.String())
	exception.PanicOnErr(err)
	if len(pe.SupR1) != 0 {
		for i := range pe.SupR1 {
			_, err = fmt.Fprintln(file, pe.SupR1[i].String())
			exception.PanicOnErr(err)
		}
	}
	if len(pe.SupR2) != 0 {
		for i := range pe.SupR2 {
			_, err = fmt.Fprintln(file, pe.SupR2[i].String())
			exception.PanicOnErr(err)
		}
	}
}

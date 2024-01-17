package sam

import (
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"io"
	"strings"
)

// SamPE is struct that contains 2 paired end Sam entries as well as a slice of Sam that represents all supplementary
// alignments for the read-pair
type SamPE struct {
	A    Sam
	B    Sam
	SupA []Sam
	SupB []Sam
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
			tmpSlice = []Sam{tmp}
		}
	}
	close(peChan)
}

// parseReads is a helper function for GoReadSamPeToChan that parses a slice of Sam that all have the same read name into a SamPE struct
func parseReads(tmpSlice []Sam) SamPE {
	var pe SamPE
	if len(tmpSlice) == 2 {
		pe = SamPE{A: tmpSlice[0], B: tmpSlice[1]}
		pe, _ = OrderSamPair(pe)
		return pe
	}
	for i := range tmpSlice {
		if tmpSlice[i].Flag&2048 == 2048 {
			if tmpSlice[i].Flag&64 == 64 {
				pe.SupA = append(pe.SupA, tmpSlice[i])
			} else if tmpSlice[i].Flag&128 == 128 {
				pe.SupB = append(pe.SupB, tmpSlice[i])
			}
		} else if pe.A.QName == "" {
			pe.A = tmpSlice[i]
		} else if pe.B.QName == "" {
			pe.B = tmpSlice[i]
		} else {
			pe.SupA = append(pe.SupA, tmpSlice[i])
		}
	}
	pe, _ = OrderSamPair(pe)
	return pe
}

// OrderSamPair returns R1, R2 in order based on Sam Flags and a bool if the pair had the proper flags. If the pair doesn't have the proper flags,
// the pair will be returned in the same order with false. Sup Slices will not be reordered.
func OrderSamPair(p SamPE) (SamPE, bool) {
	if p.A.Flag&64 == 64 && p.B.Flag&128 == 128 {
		return p, true
	} else if p.A.Flag&128 == 128 && p.B.Flag&64 == 64 {
		return SamPE{A: p.B, B: p.A, SupA: p.SupA, SupB: p.SupB}, true
	}
	return p, false
}

// WriteSamPeToFileHandle writes a both reads and any supplementary alignments from a SamPE struct to a fileio EasyWriter
// TODO: Add bam writing functionality
func WriteSamPeToFileHandle(file io.Writer, pe SamPE) {
	_, err := fmt.Fprintln(file, ToString(pe.A))
	exception.PanicOnErr(err)
	_, err = fmt.Fprintln(file, ToString(pe.B))
	exception.PanicOnErr(err)
	if len(pe.SupA) != 0 {
		for i := range pe.SupA {
			_, err = fmt.Fprintln(file, ToString(pe.SupA[i]))
			exception.PanicOnErr(err)
		}
	}
	if len(pe.SupB) != 0 {
		for i := range pe.SupB {
			_, err = fmt.Fprintln(file, ToString(pe.SupB[i]))
			exception.PanicOnErr(err)
		}
	}
}

/*		done = false
		a, full = <-data
		b, full = <-data
		for !done {
			switch {
			case full, a.QName == b.QName:
				peChan <- SamPE{a, b}
				done = true
			case a.QName != b.QName:
				a = b
				b, full = <-data
			}
		}*/

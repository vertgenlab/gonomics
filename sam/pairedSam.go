package sam

import "strings"

type SamPE struct {
	A Sam
	B Sam
}

// GoReadSamPeToChan takes a read-name sorted sam file and returns a chanel of SamPE and a Header
func GoReadSamPeToChan(filename string) (<-chan SamPE, Header) {

	//start by reading the file into a chanel of
	data := make(chan Sam, 1000)
	header := make(chan Header)
	peChan := make(chan SamPE, 1000)

	if strings.HasSuffix(filename, ".bam") {
		go readBamToChan(filename, data, header)
	} else {
		go readSamToChan(filename, data, header)
	}

	go combineEnds(peChan, data)

	return peChan, <-header
}

func combineEnds(peChan chan<- SamPE, data <-chan Sam) {
	var full bool = true
	var done bool
	var a, b Sam

	for full {
		done = false
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
		}
	}
	close(peChan)
}

package sam

import "github.com/vertgenlab/gonomics/convert"

func RemovePcrDups(samChan <-chan Sam, pe bool) {
	if pe {
		removePcrDupsPE(samChan)
	}
}

func removePcrDupsPE(samChan <-chan Sam) {
	var full bool = true
	var r1, r2 Sam

	for full == true {
		r1, full = <-samChan
		if r1.Flag&64 == 64 { //first in pair
			r2, full = <-samChan
			if r2.Flag&128 == 128 && r1.QName == r2.QName {
				convert.SamPairedEndToBed(r1, r2)
			}
		}
	}
}

package sam

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
		log.Fatalln("ERROR in IsRead1: sam entry didn't have either the read 1 or 2 flag set. Check that your sam file is paired end.")
	}
	return r1
}

func SamToPeSamCheckFlag(a, b Sam) SamPE {
	var readPair SamPE
	if a.Flag&0x40 == 0x40 && b.Flag&0x80 == 0x80 {
		readPair = SamPE{
			R1: a,
			R2: b,
		}
	} else if a.Flag&0x80 == 0x80 && b.Flag&0x40 == 0x40 {
		readPair = SamPE{
			R1: b,
			R2: a,
		}
	}
	log.Fatalln("Error in SamToPeCheckFlag: the input Sam entries did not have the proper flags indicating read pairing.")
	return readPair
}

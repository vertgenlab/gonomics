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
}

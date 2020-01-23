package sam

func flagTestBit(num int64, bit int64) bool {
	var dummy int64 = bit & num //bitwise AND
	return dummy == 0
}

func IsPaired(sam *SamAln) bool {
	return flagTestBit(sam.Flag, 1)
}

func ProperlyAligned(sam *SamAln) bool {
	return flagTestBit(sam.Flag, 2)
}

func IsUnmapped(sam *SamAln) bool {
	return flagTestBit(sam.Flag, 4)
}

func MateIsUnmapped(sam *SamAln) bool {
	return flagTestBit(sam.Flag, 8)
}

func IsPosStrand(sam *SamAln) bool {
	return flagTestBit(sam.Flag, 16)
}

func MateIsPosStrand(sam *SamAln) bool {
	return flagTestBit(sam.Flag, 32)
}

func IsForwardRead(sam *SamAln) bool {
	return flagTestBit(sam.Flag, 64)
}

func IsReverseRead(sam *SamAln) bool {
	return flagTestBit(sam.Flag, 128)
}

func IsNotPrimaryAlign(sam *SamAln) bool {
	return flagTestBit(sam.Flag, 265)
}

func ReadFailsQC(sam *SamAln) bool {
	return flagTestBit(sam.Flag, 512)
}

func IsDuplicate(sam *SamAln) bool {
	return flagTestBit(sam.Flag, 1024)
}

func IsSupplementaryAlignment(sam *SamAln) bool {
	return flagTestBit(sam.Flag, 2048)
}

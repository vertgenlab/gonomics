package sam

func IsPosStrand(sam *SamAln) bool {
	return flagIsPosStrand(sam.Flag)
}

func flagIsPosStrand(num int64) bool {
	var bit5 int64 = 16
	var dummy int64 = bit5 & num //bitwise AND
	return dummy == 0
}

func ProperlyAligned(sam *SamAln) bool {
	return flagIsProperlyAligned(sam.Flag)
}

func flagIsProperlyAligned(num int64) bool {
	var bit2 int64 = 2
	var dummy int64 = bit2 & num //bitwise AND
	return dummy == 0
}

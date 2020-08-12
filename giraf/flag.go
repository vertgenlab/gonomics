package giraf

func flagTestBit(num uint16, bit uint16) bool {
	var dummy uint16 = bit & num //bitwise AND
	return dummy == 0
}

func ProperlyAligned(giraf *Giraf) bool {
	return flagTestBit(giraf.Flag, 2)
}

func IsUnmapped(giraf *Giraf) bool {
	return flagTestBit(giraf.Flag, 4)
}

func IsPosStrand(giraf *Giraf) bool {
	return flagTestBit(giraf.Flag, 16)
}

func IsForwardRead(giraf *Giraf) bool {
	return flagTestBit(giraf.Flag, 64)
}

func IsReverseRead(giraf *Giraf) bool {
	return flagTestBit(giraf.Flag, 128)
}

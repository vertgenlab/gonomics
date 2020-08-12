package giraf

// flagTestBit test a single bit in a 16 bit flag
func flagTestBit(num uint16, bit uint16) bool {
	return bit & num == 0
}

// ProperlyAligned returns true if the aligner determined this read was aligned correctly
func ProperlyAligned(giraf *Giraf) bool {
	return flagTestBit(giraf.Flag, 2)
}

// IsUnmapped returns true if the read is unmapped
func IsUnmapped(giraf *Giraf) bool {
	return flagTestBit(giraf.Flag, 4)
}

// IsPosStrand returns true if the read aligns to the positive strand
func IsPosStrand(giraf *Giraf) bool {
	return flagTestBit(giraf.Flag, 16)
}

// IsForwardRead returns true if the read is the forward read in a pair of reads
func IsForwardRead(giraf *Giraf) bool {
	return flagTestBit(giraf.Flag, 64)
}

// IsReverseRead returns true if the read is the reverse read in a pair of reads
func IsReverseRead(giraf *Giraf) bool {
	return flagTestBit(giraf.Flag, 128)
}

package giraf

// flagTestBit test a single bit in a 16 bit flag.
func flagTestBit(num uint8, bit uint8) bool {
	return bit&num == 0
}

// ProperlyAligned returns true if the aligner determined this read was aligned correctly.
func ProperlyAligned(giraf *Giraf) bool {
	return flagTestBit(giraf.Flag, 1) // 2 in sam
}

// IsUnmapped returns true if the read is unmapped.
func IsUnmapped(giraf *Giraf) bool {
	return flagTestBit(giraf.Flag, 2) // 4 in sam
}

// IsPosStrand returns true if the read aligns to the positive strand.
func IsPosStrand(giraf *Giraf) bool {
	return flagTestBit(giraf.Flag, 4) // 16 in sam
}

// IsForwardRead returns true if the read is the forward read in a pair of reads.
func IsForwardRead(giraf *Giraf) bool {
	return flagTestBit(giraf.Flag, 8) // 64 in sam
}

// IsReverseRead returns true if the read is the reverse read in a pair of reads.
func IsReverseRead(giraf *Giraf) bool {
	return !flagTestBit(giraf.Flag, 8) // 128 in sam
}

// IsPaired returns true if the read has a mate.
func IsPaired(giraf *Giraf) bool {
	return flagTestBit(giraf.Flag, 16) // 1 in sam
}

// HasNamePrefix is reserved for binary giraf processing and returns true if read follows a defined prefix pattern.
func HasNamePrefix(giraf *Giraf) bool {
	return flagTestBit(giraf.Flag, 32) // for binary giraf
}

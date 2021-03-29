package sam

// flagTestBit tests a single bit in the sam Flag.
func flagTestBit(num uint16, bit uint16) bool {
	var dummy uint16 = bit & num //bitwise AND
	return dummy == bit
}

// IsPaired returns true if the input record has a mate pair.
func IsPaired(sam Sam) bool {
	return flagTestBit(sam.Flag, 1)
}

// ProperlyAligned returns true if the input record is properly
// aligned according to the aligner.
func ProperlyAligned(sam Sam) bool {
	return flagTestBit(sam.Flag, 2)
}

// IsUnmapped returns true if the input record is unmapped.
func IsUnmapped(sam Sam) bool {
	return flagTestBit(sam.Flag, 4)
}

// MateIsUnmapped returns true if the input records mate pair
// is unmapped.
func MateIsUnmapped(sam Sam) bool {
	return flagTestBit(sam.Flag, 8)
}

// IsPosStrand returns true if the input record aligns to the
// positive strand.
func IsPosStrand(sam Sam) bool {
	return flagTestBit(sam.Flag, 16)
}

// MateIsPosStrand returns true if the input records mate pair
// aligns to the positive strand.
func MateIsPosStrand(sam Sam) bool {
	return flagTestBit(sam.Flag, 32)
}

// IsForwardRead returns true if the input record is the forward
// read in a mate pair.
func IsForwardRead(sam Sam) bool {
	return flagTestBit(sam.Flag, 64)
}

// IsReverseRead returns true if the input record is the reverse
// read in a mate pair.
func IsReverseRead(sam Sam) bool {
	return flagTestBit(sam.Flag, 128)
}

// IsNotPrimaryAlign returns true if the input record is a secondary
// alignment according to the aligner.
func IsNotPrimaryAlign(sam Sam) bool {
	return flagTestBit(sam.Flag, 256)
}

// ReadFailsQc returns true if the input record failed quality control.
func ReadFailsQc(sam Sam) bool {
	return flagTestBit(sam.Flag, 512)
}

// IsDuplicate returns true if the input record was identified as a
// duplicate read by the aligner.
func IsDuplicate(sam Sam) bool {
	return flagTestBit(sam.Flag, 1024)
}

// IsSupplementaryAlign returns true if the input record was
// identified as a supplementary alignment by the aligner.
func IsSupplementaryAlign(sam Sam) bool {
	return flagTestBit(sam.Flag, 2048)
}

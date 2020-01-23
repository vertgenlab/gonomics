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

func getBitFlag(flag int64) []string {
	bitString := strconv.FormatInt(flag, 2)
	bits := strings.Split(bitString, "")
	return bits
}

// Note: input pos is in 1 base
func testFlag(flag int64, pos int) bool {
	bits := getBitFlag(flag)
	if len(bits) < pos {
		return false
	} else if bits[pos-1] == "1" {
		return true
	} else {
		return false
	}
}

// Parses flag to determine if read is paired
func (read *SamAln) IsPaired() bool {
	return testFlag(read.Flag, 12)
}

// Parses flag to determine if both the read and its mate are mapped
func (read *SamAln) PairIsMapped() bool {
	return testFlag(read.Flag, 11)
}

// Parses flag to determine if read is mapped
func (read *SamAln) IsMapped() bool {
	return !testFlag(read.Flag, 10)
}

// Parses flag to determine if mate is mapped
func (read *SamAln) MateIsMapped() bool {
	return !testFlag(read.Flag, 9)
}

// Parses flag to determine if read is mapped to the forward (true) or reverse (false) strand
func (read *SamAln) IsPosStrand() bool {
	return !testFlag(read.Flag, 8)
}

// Parses flag to determine if mate is mapped to the forward (true) or reverse (false) strand
func (read *SamAln) MateIsPosStrand() bool {
	return !testFlag(read.Flag, 7)
}

// Parses flag to determine if read is the forward read
func (read *SamAln) IsForwardRead() bool {
	return testFlag(read.Flag, 6)
}

// Parses flag to determine if read is the reverse read
func (read *SamAln) IsReverseRead() bool {
	return testFlag(read.Flag, 5)
}

// Parses flag to determine if read is the primary alignment
func (read *SamAln) IsPrimaryAlign() bool {
	return !testFlag(read.Flag, 4)
}

// Parses flag to determine if read passes the quality checks
func (read *SamAln) PassesQC() bool {
	return !testFlag(read.Flag, 3)
}

// Parses flag to determine if read is a duplicate
func (read *SamAln) IsDuplicate() bool {
	return testFlag(read.Flag, 2)
}

// Parses flag to determine if read is a supplementary alignment
func (read *SamAln) IsSupplementaryAlignment() bool {
	return testFlag(read.Flag, 1)
}
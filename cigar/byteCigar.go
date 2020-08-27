package cigar

import (
	"github.com/vertgenlab/gonomics/common"
	"log"
	"strconv"
	"strings"
)

// ByteCigar struct encodes sequence comparison operations and includes run length info.
type ByteCigar struct {
	RunLen uint16
	Op     byte
}

// Defined const for byte cigar
const (
	Match     byte = 'M'
	Insertion byte = 'I'
	Deletion  byte = 'D'
	N         byte = 'N'
	SoftClip  byte = 'S'
	HardClip  byte = 'H'
	Padded    byte = 'P'
	Equal     byte = '='
	Mismatch  byte = 'X'
	Unknown   byte = '*'
)

// LookUpCigByte is a helper function decode uint32 into a byte cigar struct.
// In the sam/bam specs: CIGAR: op len<<4|op. Hash is as follows: ‘MIDNSHP=X’→‘012345678’
func lookUpCigByte(op uint32) byte {
	switch op {
	case 0:
		return Match
	case 1:
		return Insertion
	case 2:
		return Deletion
	case 3:
		return N
	case 4:
		return SoftClip
	case 5:
		return HardClip
	case 6:
		return Padded
	case 7:
		return Equal
	case 8:
		return Mismatch
	default:
		log.Fatalf("Error: cound not identify input byte")
		return 0
	}
}

// ReadToBytesCigar will process a byte slice and define a small
func ReadToBytesCigar(cigar []byte) []ByteCigar {
	if cigar[0] == '*' {
		return nil
	}
	var ans []ByteCigar = make([]ByteCigar, 0, 1)
	var lastNum int = 0
	for i := 0; i < len(cigar); i++ {
		if IsValidCigar(cigar[i]) {
			ans = append(ans, ByteCigar{RunLen: common.StringToUint16(string(cigar[lastNum:i])), Op: cigar[i]})
			lastNum = i + 1
		}
	}
	return ans
}

// IsValidCigar will perform a check to make sure op is a valid byte.
func IsValidCigar(op byte) bool {
	switch op {
	case 'M':
		return true
	case 'I':
		return true
	case 'D':
		return true
	case 'N':
		return true
	case 'S':
		return true
	case 'H':
		return true
	case 'P':
		return true
	case '=':
		return true
	case 'X':
		return true
	case '*':
		return true
	default:
		return false
	}
	return false
}

// ByteCigarToString will process the cigar byte struct and parse and/or convert the data into a string
func ByteCigarToString(cigar []ByteCigar) string {
	if len(cigar) == 0 || cigar == nil {
		return "*"
	}
	var str strings.Builder
	var err error
	for _, c := range cigar {
		_, err = str.WriteString(strconv.Itoa(int(c.RunLen)))
		common.ExitIfError(err)
		err = str.WriteByte(c.Op)
		common.ExitIfError(err)
	}
	return str.String()
}

// ByteMatrixTrace will trace smith-waterman matrix alignment and return one of 3 cigar Op's.
// M: matches or mismatches, I: insertions, D: for deletions.
func ByteMatrixTrace(a int64, b int64, c int64) (int64, byte) {
	if a >= b && a >= c {
		return a, 'M'
	} else if b >= c {
		return b, 'I'
	} else {
		return c, 'D'
	}
}

// ReverseBytesCigar cigar will reverse the order of a cigar slice. Typically performed after matrix traceback
// from a local alignment.
func ReverseBytesCigar(alpha []ByteCigar) {
	var i, off int
	for i, off = len(alpha)/2-1, len(alpha)-1-i; i >= 0; i-- {
		alpha[i], alpha[off] = alpha[off], alpha[i]
	}
}

//QueryLength calculates the length of the query read from a slice of Cigar structs.
func QueryRunLen(c []ByteCigar) int {
	if c == nil {
		return 0
	}
	var ans uint16
	for _, v := range c {
		if v.Op == Match || v.Op == Insertion || v.Op == SoftClip || v.Op == Equal || v.Op == Mismatch {
			ans = ans + v.RunLen
		}
	}
	return int(ans)
}

// CatByteCigar will concatenate two cigar slices into one merged
func CatByteCigar(cigs []ByteCigar, newCigs []ByteCigar) []ByteCigar {
	if len(newCigs) == 0 || newCigs == nil {
		return cigs
	} else if len(cigs) == 0 {
		return newCigs
	} else {
		cigs = AddCigarByte(cigs, newCigs[0])
		cigs = append(cigs, newCigs[1:]...)
		return cigs
	}
}

// AddCigarByte will add append a cigar byte to an existing slice. The function
// will perform a check on the tail of the slice and incurment the run length if the
// cigar Op values are the same.
func AddCigarByte(cigs []ByteCigar, newCig ByteCigar) []ByteCigar {
	if len(cigs) == 0 {
		cigs = append(cigs, newCig)
	} else if cigs[len(cigs)-1].Op == newCig.Op {
		cigs[len(cigs)-1].RunLen += newCig.RunLen
	} else {
		cigs = append(cigs, newCig)
	}
	return cigs
}

// MatrixSetup will allocate memory for smith-waterman matrix to be used with byte cigar opertations and trace back.
func MatrixSetup(size int) ([][]int64, [][]byte) {
	m := make([][]int64, size)
	trace := make([][]byte, size)
	for idx := range m {
		m[idx] = make([]int64, size)
		trace[idx] = make([]byte, size)
	}
	return m, trace
}

// Uint32ToByteCigar will process a uint32 slice and decode each number into a byte cigar struct.
// Due to a limitation of the BAM format, CIGAR operation lengths are limited to 2^28-1
func Uint32ToByteCigar(cigar []uint32) []ByteCigar {
	var answer []ByteCigar = make([]ByteCigar, len(cigar))
	for i := 0; i < len(cigar); i++ {
		answer[i] = ByteCigar{RunLen: uint16(cigar[i] >> 4), Op: lookUpCigByte(cigar[i] & 0xf)}
	}
	return answer
}

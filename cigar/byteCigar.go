package cigar

import (
	"log"
	"strconv"
	"strings"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/numbers/parse"
)

// ByteCigar struct encodes sequence comparison operations and includes run length info.
type ByteCigar struct {
	RunLen uint16
	Op     byte
}

// Defined const for byte cigar.
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
// In the sam/bam specs: CIGAR: op len<<4|op. Hash is as follows: ‘MIDNSHP=X’→‘012345678’.
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

// lookUpUint32 will return a uint32 representation of a cigar Op.
func lookUpUint32(op byte) uint32 {
	switch op {
	case Match:
		return 0
	case Insertion:
		return 1
	case Deletion:
		return 2
	case N:
		return 3
	case SoftClip:
		return 4
	case HardClip:
		return 5
	case Padded:
		return 6
	case Equal:
		return 7
	case Mismatch:
		return 8
	default:
		log.Fatalf("Error: cound not identify input byte")
		return 0
	}
}

// ReadToBytesCigar will process a byte slice and define a small.
func ReadToBytesCigar(cigar []byte) []ByteCigar {
	if cigar[0] == '*' {
		return nil
	}
	var ans []ByteCigar = make([]ByteCigar, 0, 1)
	var lastNum int = 0
	for i := 0; i < len(cigar); i++ {
		if IsValidCigar(cigar[i]) {
			ans = append(ans, ByteCigar{RunLen: parse.StringToUint16(string(cigar[lastNum:i])), Op: cigar[i]})
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
}

// ByteCigarToString will process the cigar byte struct and parse and/or convert the data into a string.
func ByteCigarToString(cigar []ByteCigar) string {
	if len(cigar) == 0 || cigar == nil {
		return "*"
	}
	var str strings.Builder
	var err error
	for _, c := range cigar {
		_, err = str.WriteString(strconv.Itoa(int(c.RunLen)))
		exception.FatalOnErr(err)
		err = str.WriteByte(c.Op)
		exception.FatalOnErr(err)
	}
	return str.String()
}

// ByteMatrixTrace will trace smith-waterman matrix alignment and return one of 3 cigar Op's.
// M: matches or mismatches, I: insertions, D: for deletions.
func ByteMatrixTrace(a int64, b int64, c int64) (int64, byte) {
	if a >= b && a >= c {
		return a, Match
	} else if b >= c {
		return b, Insertion
	} else {
		return c, Deletion
	}
}

func TraceMatrixExtension(prev int64, a int64, b int64, c int64) (int64, byte) {
	if a >= b && a >= c {
		if a > prev {
			return a, Equal
		} else {
			return a, Mismatch
		}
	} else if b >= c {
		return b, Insertion
	} else {
		return c, Deletion
	}
}

// ReverseBytesCigar cigar will reverse the order of a cigar slice. Typically performed after matrix traceback
// from a local alignment.
func ReverseBytesCigar(alpha []ByteCigar) {
	var i, off int
	for i = len(alpha)/2 - 1; i >= 0; i-- {
		off = len(alpha) - 1 - i
		alpha[i], alpha[off] = alpha[off], alpha[i]
	}
}

// QueryRunLen calculates the length of the query read from a slice of ByteCigar structs.
func QueryRunLen(c []ByteCigar) int {
    if c == nil {
        return 0 // Return 0 if the slice is nil
    }

    var ans int // Use int to prevent overflow issues that might occur with uint16
    for _, v := range c {
        switch v.Op {
        case Match, Insertion, SoftClip, Equal, Mismatch:
            ans += int(v.RunLen) // Safely cast to int to avoid overflow
        }
    }
    return ans
}

// AddCigarByte appends or merges a new ByteCigar to the slice of ByteCigar.
func AddCigarByte(cigs []ByteCigar, newCig ByteCigar) []ByteCigar {
	if len(cigs) == 0 {
		return append(cigs, newCig)
	} else if cigs[len(cigs)-1].Op == newCig.Op {
		cigs[len(cigs)-1].RunLen += newCig.RunLen
		return cigs
	}
	return append(cigs, newCig)
}

// CatByteCigar will concatenate two cigar slices into one merged.
func CatByteCigar(cigs []ByteCigar, newCigs []ByteCigar) []ByteCigar {
	if len(newCigs) == 0 {
		return cigs
	}
	if len(cigs) == 0 {
		return append([]ByteCigar(nil), newCigs...) // Create a copy of newCigs to ensure immutability
	}

	// Merge the first of the newCigs with the last of the cigs if possible
	cigs = AddCigarByte(cigs, newCigs[0])
	
	// Append the rest of newCigs
	if len(newCigs) > 1 {
		cigs = append(cigs, newCigs[1:]...)
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
// CIGAR operation lengths are limited to 2^28-1 in the current sam/bam formats.
func Uint32ToByteCigar(cigar []uint32) []ByteCigar {
	var answer []ByteCigar = make([]ByteCigar, len(cigar))
	for i := 0; i < len(cigar); i++ {
		answer[i] = ByteCigar{RunLen: uint16(cigar[i] >> 4), Op: lookUpCigByte(cigar[i] & 0xf)}
	}
	return answer
}

// ByteCigarToUint32 will convert a slice of []ByteCigar to a slice of []uint32.
func ByteCigarToUint32(cigar []ByteCigar) []uint32 {
	var answer []uint32 = make([]uint32, len(cigar))
	for i := 0; i < len(cigar); i++ {
		answer[i] = lookUpUint32(cigar[i].Op) | uint32(cigar[i].RunLen)<<4
	}
	return answer
}

func SoftClipBases(front int, lengthOfRead int, cig []ByteCigar) []ByteCigar {
    runLen := QueryRunLen(cig) // Ensure this function counts only alignment-contributing operations

    // Direct return if the current run length already meets or exceeds the required read length
    if runLen >= lengthOfRead {
        return cig
    }

    // Calculate total needed soft clips
    totalNeeded := lengthOfRead - runLen
    frontClips := uint16(front)
    endClips := uint16(totalNeeded - front)

    // Check if there's any need to append new clips
    if frontClips == 0 && endClips == 0 {
        return cig
    }

    // Prepare new slice for answer
    answerCapacity := len(cig) + 2 // Potentially one soft clip at each end
    answer := make([]ByteCigar, 0, answerCapacity)

    // Append front soft clip if needed
    if frontClips > 0 {
        answer = append(answer, ByteCigar{RunLen: frontClips, Op: 'S'})
    }

    // Append existing cigar elements
    answer = append(answer, cig...)

    // Append end soft clip if needed
    if endClips > 0 {
        answer = append(answer, ByteCigar{RunLen: endClips, Op: 'S'})
    }

    return answer
}


package cigar

// import (
// 	"log"
// 	"strconv"
// 	"strings"

// 	"github.com/vertgenlab/gonomics/exception"
// 	"github.com/vertgenlab/gonomics/numbers/parse"
// )

// // ByteMatrixTrace will trace smith-waterman matrix alignment and return one of 3 cigar Op's.
// // M: matches or mismatches, I: insertions, D: for deletions.
// func ByteMatrixTrace(a int64, b int64, c int64) (int64, byte) {
// 	if a >= b && a >= c {
// 		return a, 'M'
// 	} else if b >= c {
// 		return b, 'I'
// 	} else {
// 		return c, 'D'
// 	}
// }

// // ReverseBytesCigar cigar will reverse the order of a cigar slice. Typically performed after matrix traceback
// // from a local alignment.
// func ReverseBytesCigar(alpha []ByteCigar) {
// 	var i, off int
// 	for i = len(alpha)/2 - 1; i >= 0; i-- {
// 		off = len(alpha) - 1 - i
// 		alpha[i], alpha[off] = alpha[off], alpha[i]
// 	}
// }

// // QueryLength calculates the length of the query read from a slice of Cigar structs.
// func QueryRunLen(c []ByteCigar) int {
// 	if c == nil {
// 		return 0
// 	}
// 	var ans uint16
// 	for _, v := range c {
// 		if v.Op == Match || v.Op == Insertion || v.Op == SoftClip || v.Op == Equal || v.Op == Mismatch {
// 			ans = ans + v.RunLen
// 		}
// 	}
// 	return int(ans)
// }

// // CatByteCigar will concatenate two cigar slices into one merged.
// func CatByteCigar(cigs []ByteCigar, newCigs []ByteCigar) []ByteCigar {
// 	if len(newCigs) == 0 || newCigs == nil {
// 		return cigs
// 	} else if len(cigs) == 0 {
// 		return newCigs
// 	} else {
// 		cigs = AddCigarByte(cigs, newCigs[0])
// 		cigs = append(cigs, newCigs[1:]...)
// 		return cigs
// 	}
// }

// // AddCigarByte will add append a cigar byte to an existing slice. The function
// // will perform a check on the tail of the slice and incurment the run length if the
// // cigar Op values are the same.
// func AddCigarByte(cigs []ByteCigar, newCig ByteCigar) []ByteCigar {
// 	if len(cigs) == 0 {
// 		cigs = append(cigs, newCig)
// 	} else if cigs[len(cigs)-1].Op == newCig.Op {
// 		cigs[len(cigs)-1].RunLen += newCig.RunLen
// 	} else {
// 		cigs = append(cigs, newCig)
// 	}
// 	return cigs
// }

// // MatrixSetup will allocate memory for smith-waterman matrix to be used with byte cigar opertations and trace back.
// func MatrixSetup(size int) ([][]int64, [][]byte) {
// 	m := make([][]int64, size)
// 	trace := make([][]byte, size)
// 	for idx := range m {
// 		m[idx] = make([]int64, size)
// 		trace[idx] = make([]byte, size)
// 	}
// 	return m, trace
// }

// // Uint32ToByteCigar will process a uint32 slice and decode each number into a byte cigar struct.
// // CIGAR operation lengths are limited to 2^28-1 in the current sam/bam formats.
// func Uint32ToByteCigar(cigar []uint32) []ByteCigar {
// 	var answer []ByteCigar = make([]ByteCigar, len(cigar))
// 	for i := 0; i < len(cigar); i++ {
// 		answer[i] = ByteCigar{RunLen: uint16(cigar[i] >> 4), Op: lookUpCigByte(cigar[i] & 0xf)}
// 	}
// 	return answer
// }

// // ByteCigarToUint32 will convert a slice of []ByteCigar to a slice of []uint32.
// func ByteCigarToUint32(cigar []ByteCigar) []uint32 {
// 	var answer []uint32 = make([]uint32, len(cigar))
// 	for i := 0; i < len(cigar); i++ {
// 		answer[i] = lookUpUint32(cigar[i].Op) | uint32(cigar[i].RunLen)<<4
// 	}
// 	return answer
// }
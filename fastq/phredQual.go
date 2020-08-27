package fastq

import (
	"github.com/vertgenlab/gonomics/common"
	"math"
	"strings"
)

// SAM format uses ascii offset of 33 to make everything start with individual characters
// without adding 33 you get values like spaces and newlines
const asciiOffset uint8 = 33

func ToQualUint8(qual []rune) []uint8 {
	var answer []uint8 = make([]uint8, len(qual))
	for i := 0; i < len(qual); i++ {
		answer[i] = uint8(qual[i]) - asciiOffset
	}
	return answer
}

func ReverseQualUint8Record(qualScore []uint8) {
	for i, j := 0, len(qualScore)-1; i <= j; i, j = i+1, j-1 {
		qualScore[i], qualScore[j] = qualScore[j], qualScore[i]
	}
}

func Uint8QualToString(qual []uint8) string {
	var answer []rune = make([]rune, len(qual))
	for i := 0; i < len(qual); i++ {
		answer[i] = rune(qual[i] + asciiOffset)
	}
	return string(answer)
}

func QualString(qual []uint8) string {
	var str strings.Builder
	var err error
	str.Grow(len(qual))
	for i := 0; i < len(qual); i++ {
		err = str.WriteByte(byte(qual[i] + asciiOffset))
		common.ExitIfError(err)
	}
	return str.String()
}

func PhredToPError(baseQual uint8) float32 {
	q := float64(baseQual)
	p := math.Pow(10, -q/10)
	return float32(p)
}

func ErrorRate(baseQual []uint8) []float32 {
	var answer []float32 = make([]float32, len(baseQual))
	for i := 0; i < len(baseQual); i++ {
		answer[i] = PhredToPError(baseQual[i])
	}
	return answer
}

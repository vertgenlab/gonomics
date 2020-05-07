package fastq

import(
	"math"
)
func ToQualUint8(qual []rune) []uint8 {
	var answer []uint8 = make([]uint8, len(qual))
	for i := 0; i < len(qual); i++ {
		// SAM format uses ascii offset of 33 to make everything start with individual characters
		// without adding 33 you get values like spaces and newlines
		var asciiOffset uint8 = 33
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
		// SAM format uses ascii offset of 33 to make everything start with individual characters
		// without adding 33 you get values like spaces and newlines
		var asciiOffset uint8 = 33
		answer[i] = rune(qual[i] + asciiOffset)
	}

	return string(answer)
}

func PhredToPError(ascii rune) float32 {
	q := float64(ascii) - 33
	p := math.Pow(10, -q/10)
	return float32(p)
}

func ErrorRate(ASCII []rune) []float32 {
	var answer []float32 = make([]float32, len(ASCII))
	for i := 0; i < len(ASCII); i++ {
		answer[i] = PhredToPError(ASCII[i])
	}
	return answer
}

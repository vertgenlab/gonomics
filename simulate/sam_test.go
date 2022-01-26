package simulate

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"testing"
)

func TestGenerateSam(t *testing.T) {
	fmt.Println(generateSamReadNoFlag("test", "chr1", []dna.Base{dna.A, dna.C, dna.G, dna.T}, -6, 15))
}

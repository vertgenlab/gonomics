// Command Group: "Linear Alignment Tools"

package main

import (
	"fmt"
	"github.com/vertgenlab/gonomics/align"
	"github.com/vertgenlab/gonomics/dna"
)

func main() {
	_, aln_test_lowMem := align.AffineGap_customizeCheckersize(dna.StringToBases("AA"), dna.StringToBases("GGGAATT"), align.HumanChimpTwoScoreMatrix, -600, -150, 10000, 10000)
	_, aln_test_highMem := align.AffineGap_highMem(dna.StringToBases("AA"), dna.StringToBases("GGGAATT"), align.HumanChimpTwoScoreMatrix, -600, -150)
	fmt.Printf("TEST LowMem species1 vs species2 cmd lowMem cigar route: %v\n", aln_test_lowMem)
	fmt.Printf("TEST HighMem species1 vs species2 cmd highMem cigar route: %v\n", aln_test_highMem)
}

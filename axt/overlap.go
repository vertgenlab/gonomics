package axt

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/common"
	"strings"
)

func OverlapAxtBed(alpha *Axt, beta *bed.Bed) bool {
	if (common.MaxInt64(alpha.RStart-1, beta.ChromStart) < common.MinInt64(alpha.REnd-1, beta.ChromEnd)) && strings.Compare(alpha.RName, beta.Chrom) == 0 {
		return true
	} else {
		return false
	}
}

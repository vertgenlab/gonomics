package popgen

import (
	"github.com/vertgenlab/gonomics/vcf"
	"testing"
	"fmt"
)

var VcfSampleToSegSiteTests = []struct {
	Vcf vcf.Vcf
	DivergenceAscertainment bool
	UnPolarized bool
	IncludeRef bool
	Expected SegSite
	ExpectedFlagBool bool
}{
	{vcf.Vcf{Chr: "chr1",
			Pos: 1,
			Ref: "A",
			Alt: []string{"T"},
			Info: "AA=A",
			Samples: []vcf.Sample{
				vcf.Sample{
					[]int16{1,0}, []bool{false, false}, []string{""},
				},
			},
	},
	false,
	false,
	false,
	SegSite{1, 2, 0},
	true,},
	{vcf.Vcf{Chr: "chr1",
		Pos: 1,
		Ref: "A",
		Alt: []string{"T"},
		Info: "AA=A",
		Samples: []vcf.Sample{
			vcf.Sample{
				[]int16{1,0}, []bool{false, false}, []string{""},
			},
		},
	},
		false,
		false,
		true,
		SegSite{1, 3, 0},
		true,},
	{vcf.Vcf{Chr: "chr1",
		Pos: 1,
		Ref: "A",
		Alt: []string{"T"},
		Info: "AA=T",
		Samples: []vcf.Sample{
			vcf.Sample{
				[]int16{1,0}, []bool{false, false}, []string{""},
			},
		},
	},
		false,
		false,
		true,
		SegSite{1, 3, 0},//expected frequency is still one because the segSite is inverted.
		true,},
	{vcf.Vcf{Chr: "chr1",
		Pos: 1,
		Ref: "A",
		Alt: []string{"T"},
		Info: "AA=T",
		Samples: []vcf.Sample{
			vcf.Sample{
				[]int16{1,0}, []bool{false, false}, []string{""},
			},
		},
	},
		true,
		false,
		true,
		SegSite{1, 3, 2},//expected frequency is still one because the segSite is inverted.
		true,},
	{vcf.Vcf{Chr: "chr1",
		Pos: 1,
		Ref: "A",
		Alt: []string{"T"},
		Info: "AA=A",
		Samples: []vcf.Sample{
			vcf.Sample{
				[]int16{1,0}, []bool{false, false}, []string{""},
			},
		},
	},
		true,
		false,
		true,
		SegSite{1, 3, 1},//expected frequency is still one because the segSite is inverted.
		true,},
}

func TestVcfSampleToSegSite(t *testing.T) {
	var output *SegSite
	var flagBool bool
	for _, v := range VcfSampleToSegSiteTests {
		output, flagBool = VcfSampleToSegSite(v.Vcf, v.DivergenceAscertainment, v.UnPolarized, v.IncludeRef)
		if flagBool != v.ExpectedFlagBool {
			t.Errorf("Error in VcfSampleToSegSite. Did not receive the expected error.")
		}
		if !segSitesAreEqual(v.Expected, *output) {
			fmt.Printf("%v.\n", *output)
			t.Errorf("Error in VcfSampleToSegSite. Output segSite was not as expected.")
		}
	}
}

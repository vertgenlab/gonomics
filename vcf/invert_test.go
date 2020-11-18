package vcf

import (
	"testing"
)

var FirstInputVcf Vcf = Vcf{Chr:"chr2", Pos: 4, Ref: "C", Alt: "T"}
var FirstExpectedVcf Vcf = Vcf{Chr:"chr2", Pos: 4, Ref: "T", Alt: "C"}

var InvertVcfTests = []struct {
	Input	Vcf
	Expected	Vcf
}{
	{FirstInputVcf, FirstExpectedVcf},
	{FirstExpectedVcf, FirstInputVcf},
}

func TestInvertVcf(t *testing.T) {
	for _, v := range InvertVcfTests {
		InvertVcf(v.Input)
		if !isEqual(&v.Input, &v.Expected) {
			t.Errorf("Error in InvertVCF. Input Ref: %s. Expected Ref: %s.", v.Input.Ref, v.Expected.Ref)
		}
	}
}

/*
var InvertGVcfTests = []struct {
	Input GVcf
	Expected GVcf
}{
	{},
}

func TestInvertGVcf(t *testing.T) {
	var answer GVcf
	for v := range InvertTests {
		answer = InvertGVcf(v.Input)
		if !EqualGVcf(answer, v.Expected) {
			log.Fatalf("Error in TestInvertGVcf.")
		}
	}
}
*/

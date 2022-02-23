package main

import (
	"github.com/vertgenlab/gonomics/vcf"
	"testing"
)

func TestExpressionFormat(t *testing.T) {
	data, header := vcf.Read("testdata/headerTest.vcf")
	testExp := "FormatF >= 1"

	tests := parseExpression(testExp, header, true, true)
	for _, v := range data {
		v = vcf.ParseFormat(v, header)
		if !passesTests(v, tests) {
			t.Errorf("problem with format expression")
		}
	}
}

func TestExpressionInfo(t *testing.T) {
	data, header := vcf.Read("testdata/headerTest.vcf")
	testExp := "InfoA < 20 ; InfoFlag"

	tests := parseExpression(testExp, header, false, true)
	for _, v := range data {
		v = vcf.ParseInfo(v, header)
		if !passesTests(v, tests) {
			t.Errorf("problem with info expression")
		}
	}
}

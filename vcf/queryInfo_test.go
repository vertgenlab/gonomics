package vcf

import (
	"fmt"
	"testing"
)

func TestParseFormat(t *testing.T) {
	data, header := Read("testdata/headerTest.vcf")
	v := ParseFormat(data[0], header)
	fmt.Println(v.parsedFormat)
}

func TestParseInfo(t *testing.T) {
	data, header := Read("testdata/headerTest.vcf")
	v := ParseInfo(data[0], header)
	fmt.Println(v.parsedInfo)
}
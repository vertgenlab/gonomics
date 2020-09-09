package giraf

import (
	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dnaThreeBit"
	"reflect"
	"testing"
)

var TestQual = []uint8{40, 5, 5, 5, 5, 5, 5, 5, 30, 20, 20, 20, 1}
var ExpectedQual = []cigar.ByteCigar{{RunLen: 1, Op: 40},
	{RunLen: 7, Op: 5}, {RunLen: 1, Op: 30},
	{RunLen: 3, Op: 20}, {RunLen: 1, Op: 1}}

func TestEncodeQual(t *testing.T) {
	answer := encodeQual(TestQual)
	if !reflect.DeepEqual(answer, ExpectedQual) {
		t.Errorf("Error with run-length encoding of quality scores")
	}
}

var TestSeq = dna.StringToBases("ACGTGGTCA")
var TestCigar = []cigar.ByteCigar{{RunLen: 1, Op: 'S'},
	{RunLen: 4, Op: '='}, {RunLen: 2, Op: 'I'},
	{RunLen: 1, Op: 'X'}, {RunLen: 3, Op: '='}}
var FancySeq = "AGTC"

func TestGetFancySeq(t *testing.T) {
	answer := getFancySeq(TestSeq, TestCigar)
	if dnaThreeBit.ToString(&answer) != FancySeq {
		t.Errorf("Error retrieving fancy seq from giraf")
	}
}

var TestNotes = []Note{{Tag: []byte("BC"), Type: 'I', Value: "TEST"},
	{Tag: []byte("AD"), Type: 'E', Value: "TEST2"}}
var BinNotes = []BinNote{{tag: [2]byte{'B', 'C'}, tagType: 'I', data: []byte("TEST")},
	{tag: [2]byte{'A', 'D'}, tagType: 'E', data: []byte("TEST2")}}

func TestEncodeNotes(t *testing.T) {
	answer := encodeNotes(TestNotes)
	if !reflect.DeepEqual(answer, BinNotes) {
		t.Errorf("Error encoding giraf notes")
	}
}

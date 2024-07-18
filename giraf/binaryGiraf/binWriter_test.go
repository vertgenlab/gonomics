package binaryGiraf

import (
	"io"
	"reflect"
	"testing"

	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/dna/dnaThreeBit"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/giraf"
)

var TestQual = []uint8{40, 5, 5, 5, 5, 5, 5, 5, 30, 20, 20, 20, 1}
var ExpectedQual = []cigar.Cigar{{RunLength: 1, Op: 40},
	{RunLength: 7, Op: 5}, {RunLength: 1, Op: 30},
	{RunLength: 3, Op: 20}, {RunLength: 1, Op: 1}}

func TestEncodeQual(t *testing.T) {
	answer := encodeQual(TestQual)
	if !reflect.DeepEqual(answer, ExpectedQual) {
		t.Errorf("Error with run-length encoding of quality scores")
	}
}

var TestSeq = dna.StringToBases("ACGTGGTCA")
var TestCigar = []cigar.Cigar{{RunLength: 1, Op: 'S'},
	{RunLength: 4, Op: '='}, {RunLength: 2, Op: 'I'},
	{RunLength: 1, Op: 'X'}, {RunLength: 3, Op: '='}}
var FancySeq = "AGTC"

func TestGetFancySeq(t *testing.T) {
	answer := getFancySeq(TestSeq, TestCigar)
	if dnaThreeBit.ToString(&answer) != FancySeq {
		t.Errorf("Error retrieving fancy seq from giraf")
	}
}

var TestNotes = []giraf.Note{{Tag: []byte("BC"), Type: 'Z', Value: "TEST" + "\000"},
	{Tag: []byte("AD"), Type: 'Z', Value: "TEST2" + "\000"}}
var ByteNotes = "BCZTEST" + "\000" + "ADZTEST2" + "\000"

func TestEncodeNotes(t *testing.T) {
	answer := notesToBytes(TestNotes)
	if !reflect.DeepEqual(answer, []byte(ByteNotes)) {
		t.Errorf("Error encoding giraf notes")
	}
}

func TestWrite(t *testing.T) {
	CompressGiraf("testdata/test.giraf", "testdata/test.giraf.fe")
	result := fileio.EasyOpen("testdata/test.giraf.fe")
	if _, err := result.BuffReader.ReadByte(); err == io.EOF {
		t.Errorf("Error: binary giraf was not written")
	}
}

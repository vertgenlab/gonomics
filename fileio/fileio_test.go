package fileio

import (
	"testing"
)

var testfile string = "testdata/smallTest"
var line1 string = "#shhhh this line is a secret"
var line2 string = "Hello World"
var line3 string = "I am a gopher"

func TestNextLine(t *testing.T) {
	file := EasyOpen(testfile)
	l1, done := NextLine(file.BuffReader)
	if l1 != line1 || done {
		t.Errorf("problem reading next line")
	}
	l2, done := NextLine(file.BuffReader)
	if l2 != line2 || done {
		t.Errorf("problem reading next line")
	}
	l3, done := NextLine(file.BuffReader)
	if l3 != line3 || done {
		t.Errorf("problem reading next line")
	}
	l4, done := NextLine(file.BuffReader)
	if l4 != "" || !done {
		t.Errorf("problem reading next line")
	}
}

func TestNextRealLine(t *testing.T) {
	file := EasyOpen(testfile)
	l1, done := NextRealLine(file.BuffReader)
	if l1 != line2 || done {
		t.Errorf("problem reading next line")
	}
	l2, done := NextRealLine(file.BuffReader)
	if l2 != line3 || done {
		t.Errorf("problem reading next line")
	}
	l3, done := NextRealLine(file.BuffReader)
	if l3 != "" || !done {
		t.Errorf("problem reading next line")
	}
}

func TestEqual(t *testing.T) {
	if !AreEqual(testfile, testfile) || !AreEqualIgnoreComments(testfile, testfile) {
		t.Errorf("problem with equal")
	}
}

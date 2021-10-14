package fileio

import (
	"os"
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

func TestRead(t *testing.T) {
	file := Read(testfile)
	if file[0] != line2 {
		t.Errorf("problem with read")
	}
	if file[1] != line3 {
		t.Errorf("problem with read")
	}
}

func TestWrite(t *testing.T) {
	var createdFile string
	createdFile = "test.txt"
	var fileContent []string = []string {line1, line2, line3}
	Write(createdFile, fileContent)
	if !AreEqual(testfile, createdFile) {
		t.Errorf("problem with fileio.Write()")
	}
	os.Remove("test.txt")
}
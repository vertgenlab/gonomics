package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"os"
	"testing"
)

func TestOboToDot(t *testing.T) {
	oboToDot("testdata/go.obo", "GO:0007275", "testdata/out.dot")

	if fileio.AreEqualIgnoreOrder("testdata/out.dot", "testdata/expected.dot") {
		err := os.Remove("testdata/out.dot")
		exception.PanicOnErr(err)
	} else {
		log.Fatal("Expected file and output file do not match.")
	}
}

package main

import (
	"fmt"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fasta"
	"log"
	"strings"
	"testing"
)

var seq1, _ = dna.StringToBases("ATGTGCAAA")
var seq2, _ = dna.StringToBases("ACCGTACGT")
var sliceFasta = []*fasta.Fasta{{"dna_seq", seq1}, {"anotherSequence", seq2}}
var g = NewGraph()
var gg = NewGraph()

func TestAdd(t *testing.T) {
	records, err := fasta.ReadNew("test.fa")
	if err != nil {
		log.Fatal(err)
	}
	a := fasta.DivideFastaAll(sliceFasta, 3)
	b := fasta.DivideFastaAll(records, 3)
	g = FillGraph(a)
	gg = FillGraph(b)
	g.String()
	gg.String()
	if strings.Compare(string(g.String()), string(gg.String())) != 0 {
		fmt.Println("Fasta files are not equal")
	} else {
		fmt.Println(g.String())
	}

}

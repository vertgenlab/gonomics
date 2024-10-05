package gff3

import (
	"testing"
)

func TestEqual(t *testing.T) {
	// Test case where Gff3 structs are Equal
	a := Gff3{
		Id:     "chr1",
		Source: "Ensembl",
		Type:   "gene",
		Start:  1000,
		End:    2000,
		Score:  0.5,
		Strand: '+',
		Phase:  ".",
		Attributes: []Tag{
			{Label: "ID", Value: "gene1"},
			{Label: "Name", Value: "MyGene"},
		},
	}
	b := Gff3{
		Id:     "chr1",
		Source: "Ensembl",
		Type:   "gene",
		Start:  1000,
		End:    2000,
		Score:  0.5,
		Strand: '+',
		Phase:  ".",
		Attributes: []Tag{
			{Label: "ID", Value: "gene1"},
			{Label: "Name", Value: "MyGene"},
		},
	}
	if !Equal(a, b) {
		t.Error("Equal(a, b) = false; want true")
	}

	// Test case where Gff3 structs are different
	b.Id = "chr2"
	if Equal(a, b) {
		t.Error("Equal(a, b) = true; want false")
	}

	// Test case with different number of attributes
	b.Id = "chr1"
	b.Attributes = []Tag{
		{Label: "ID", Value: "gene1"},
	}
	if Equal(a, b) {
		t.Error("Equal(a, b) = true; want false")
	}

	// Test case with different attributes
	b.Attributes = []Tag{
		{Label: "ID", Value: "gene2"},
		{Label: "Name", Value: "MyGene"},
	}
	if Equal(a, b) {
		t.Error("Equal(a, b) = true; want false")
	}
}

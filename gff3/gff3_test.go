package gff3

import (
	"os"
	"testing"
)

func TestToGff3(t *testing.T) {
	// Test case with all fields populated
	line := "chr1\tEnsembl\tgene\t1000\t2000\t0.5\t+\t.\tID=gene1;Name=MyGene"
	expected := Gff3{
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
	actual := ToGff3(line)
	if !Equal(actual, expected) {
		t.Errorf("Error: ToGff3(%s) = %s; want %s", line, Gff3ToString(actual), Gff3ToString(expected))
	}

	// Test case with missing score
	line = "chr1\tEnsembl\tgene\t1000\t2000\t.\t+\t.\tID=gene1"
	expected = Gff3{
		Id:     "chr1",
		Source: "Ensembl",
		Type:   "gene",
		Start:  1000,
		End:    2000,
		Score:  0, // Expecting 0 for missing score
		Strand: '+',
		Phase:  ".",
		Attributes: []Tag{
			{Label: "ID", Value: "gene1"},
		},
	}
	actual = ToGff3(line)
	if !Equal(actual, expected) {
		t.Errorf("Error: ToGff3(%s) = %v; want %v", line, actual, expected)
	}

	invalidLines := []string{
		"chr1\tEnsembl\tgene\t1000\t2000\t.\t+\tID=gene1",    // Missing a field
		"chr1\tEnsembl\tgene\t1000\t2000\t.\t?\t.\tID=gene1", // Invalid strand
		"chr1 Ensembl gene 1000 2000 . + . ID=gene1",         // Incorrect delimiter
	}
	for _, line := range invalidLines {
		// Create a closure to capture the current line
		func(line string) {
			defer func() {
				if r := recover(); r == nil {
					t.Errorf("Error: ToGff3 with invalid line did not panic: %s", line)
				}
			}()
			ToGff3(line) // Should panic
		}(line)
	}
}

func TestGff3ToString(t *testing.T) {
	// Test case with all fields populated
	gff3 := Gff3{
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
	expected := "chr1\tEnsembl\tgene\t1000\t2000\t0.5\t+\t.\tID=gene1;Name=MyGene"
	actual := Gff3ToString(gff3)
	if actual != expected {
		t.Errorf("Error: Gff3ToString(%v) = %s; want %s", gff3, actual, expected)
	}

	// Test case with missing score
	gff3 = Gff3{
		Id:     "chr1",
		Source: "Ensembl",
		Type:   "gene",
		Start:  1000,
		End:    2000,
		Score:  0,
		Strand: '+',
		Phase:  ".",
		Attributes: []Tag{
			{Label: "ID", Value: "gene1"},
		},
	}
	expected = "chr1\tEnsembl\tgene\t1000\t2000\t.\t+\t.\tID=gene1"
	actual = Gff3ToString(gff3)
	if actual != expected {
		t.Errorf("Error: Gff3ToString(%v) = %s; want %s", gff3, actual, expected)
	}
}
func TestGff3ReadWrite(t *testing.T) {
	// Create a temporary file
	tmpfile := "testdata/tmp.gff3"

	defer os.Remove(tmpfile)

	// Write some test Gff3 data to the file
	testGff3 := Read("testdata/stickleback_v5_ensembl_genes.gff3.gz")

	Write(tmpfile, testGff3)

	// Read the Gff3 data back from the file
	readGff3 := Read(tmpfile)

	// Compare the written and read data
	if len(readGff3) != len(testGff3) {
		t.Errorf("Error: Expected %d Gff3 records, got %d", len(testGff3), len(readGff3))
	}

	for i, expected := range testGff3 {
		got := readGff3[i]
		if Gff3ToString(got) != Gff3ToString(expected) {
			t.Errorf("Error: Record %d mismatch:\nExpected: %s\nGot: %s", i, Gff3ToString(expected), Gff3ToString(got))
		}
	}
}

func TestParseTags(t *testing.T) {
	// Test case with multiple tag-value pairs
	attrStr := "ID=gene1;Name=MyGene;Alias=G1"
	expected := []Tag{
		{Label: "ID", Value: "gene1"},
		{Label: "Name", Value: "MyGene"},
		{Label: "Alias", Value: "G1"},
	}
	actual := parseTags(attrStr)
	if tagsToString(actual) != tagsToString(expected) {
		t.Errorf("Error: parseTags(%s) = %s; want %s", attrStr, tagsToString(actual), tagsToString(expected))
	}

	// Test case with a single tag-value pair
	attrStr = "ID=gene1"
	expected = []Tag{
		{Label: "ID", Value: "gene1"},
	}
	actual = parseTags(attrStr)
	if tagsToString(actual) != tagsToString(expected) {
		t.Errorf("Error: parseTags(%s) = %s; want %s", attrStr, tagsToString(actual), tagsToString(expected))
	}

	// Test case with empty attributes string
	attrStr = ""
	expected = []Tag{}
	actual = parseTags(attrStr)
	if tagsToString(actual) != tagsToString(expected) {
		t.Errorf("Error: parseTags(%s) = %s; want %s", attrStr, tagsToString(actual), tagsToString(expected))
	}

	// Test case with an invalid attribute pair (should panic)
	defer func() {
		if r := recover(); r == nil {
			t.Errorf("Error: parseTags with invalid attribute pair did not panic")
		}
	}()
	attrStr = "ID=gene1;Name" // Missing '='
	parseTags(attrStr)
}

func TestTagsToString(t *testing.T) {
	// Test case with multiple tag-value pairs
	tags := []Tag{
		{Label: "ID", Value: "gene1"},
		{Label: "Name", Value: "MyGene"},
		{Label: "Alias", Value: "G1"},
	}
	expected := "ID=gene1;Name=MyGene;Alias=G1"
	actual := tagsToString(tags)
	if actual != expected {
		t.Errorf("Error: tagsToString(%v) = %s; want %s", tags, actual, expected)
	}

	// Test case with a single tag-value pair
	tags = []Tag{
		{Label: "ID", Value: "gene1"},
	}
	expected = "ID=gene1"
	actual = tagsToString(tags)
	if actual != expected {
		t.Errorf("Error: tagsToString(%v) = %s; want %s", tags, actual, expected)
	}

	// Test case with empty tags slice
	tags = []Tag{}
	expected = ""
	actual = tagsToString(tags)
	if actual != expected {
		t.Errorf("Error: tagsToString(%v) = %s; want %s", tags, actual, expected)
	}
}

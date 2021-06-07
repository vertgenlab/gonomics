// Package vcf contains functions for reading, writing, and manipulating VCF format files. More information on the VCF file format can be found
// in its official documentation at https://samtools.github.io/hts-specs/VCFv4.3.pdf. This file is parsed into a struct containing header information
// as well as a Vcf struct containing the information from each data line.
package vcf

const Version = "VCFv4.3"

// Vcf contains information for each line of a VCF format file, corresponding to variants at one position of a reference genome.
type Vcf struct {
	Chr     string
	Pos     int
	Id      string
	Ref     string
	Alt     []string
	Qual    float64
	Filter  string
	Info    string
	Format  []string
	Samples []GenomeSample

	// parsedInfo and parsedFormat store data in Info and Format fields keyed by ID.
	// nil until initialized with ParseInfo and ParseFormat respectively.
	// These fields should only be accessed by Query functions (e.g. QueryInt).
	parsedInfo   map[string]interface{}
	parsedFormat map[string]interface{}
}

// GenomeSample is a substruct of Vcf, and contains information about each sample represented in a VCF line.
// AlleleOne and AlleleTwo are set to -1 if no genotype data is present. AlleleTwo is set to -1 for hemizygous sites.
type GenomeSample struct {
	AlleleOne  int16    // First allele in genotype, 0 for reference, 1 for Alt[0], 2 for Alt[1], etc.
	AlleleTwo  int16    // Second allele in genotype, same Number format as above.
	Phased     bool     // True for phased genotype, false for unphased.
	FormatData []string // FormatData contains additional sample fields after the genotype, which are parsed into a slice delimited by colons. Currently contains a dummy empty string in FormatData[0] corresponding to "GT" in Format, so indices in FormatData will match the indices in Format.
}

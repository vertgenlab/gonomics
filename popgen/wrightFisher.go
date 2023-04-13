package popgen

import (
	"fmt"
	"io"
	"strconv"

	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/vcf"
)

type WrightFisherSettings struct {
	PopSize         int
	MutRate         float64
	NumGen          int
	GenomeSize      int
	RFitness        float64
	GcContent       float64
	InitFreq        string
	FitnessString   string
	SetSeed         int64
	Verbose         bool
	Fasta           bool
	Vcf             bool
	AncestralAllele string
}

type WrightFisherPopData struct {
	Fasta     []fasta.Fasta
	Vcf       vcf.Vcf
	Meta      []string
	Freq      [][][]float64
	Ancestral []string //Ancestral states of each allele on each site (aka replicates)
	Settings  WrightFisherSettings
}

/*
WriteTsv() writes a .tsv file based on WrightFisherPopData, including:
Comments starts with # that include metadata about the parameters of the simulation
Header line: Gen	Site	Freq.A	Freq.C	Freq.G	Freq.T	Ancestral
Frequencies table.
*/
func WriteTSV(outFile string, wf WrightFisherPopData) {
	file := fileio.EasyCreate(outFile)
	header := []string{
		"Gen",
		"Site",
		"Freq.A",
		"Freq.C",
		"Freq.G",
		"Freq.T",
		"Ancestral",
	}

	writeMeta(file, wf.Meta)                        // Write metadata
	writeEachLine(file, header)                     // Write the header line
	writeToFileHandle(file, floatArrayToString(wf)) // Write the main frequencies table

	exception.PanicOnErr(file.Close()) // Close the file
}

/*
floatArrayTostring() converts the all frequency array from wf to a 2D slice of string
2D array, zero-based, [generation and site][base].
*/
func floatArrayToString(wf WrightFisherPopData) [][]string {
	// nrow = product of number of generation and number of site
	nrow := (wf.Settings.NumGen + 1) * wf.Settings.GenomeSize
	ncol := 7
	answer := make([][]string, nrow)

	var t, s, r, b int // iterator for generation, site, row, and base

	for t = 0; t <= wf.Settings.NumGen; t++ {
		for s = 0; s < wf.Settings.GenomeSize; s++ {
			r = (t * wf.Settings.GenomeSize) + s // Keeps track of index based on generation and site
			answer[r] = make([]string, ncol)
			answer[r][0] = strconv.Itoa(t)
			answer[r][1] = strconv.Itoa(s)
			for b = 0; b < 4; b++ {
				answer[r][2+b] = strconv.FormatFloat(wf.Freq[t][s][b], 'f', 5, 64)
			}
			answer[r][6] = wf.Ancestral[s]
		}
	}
	return answer
}

/*
writeToFileHandle() feeds one slice of string at a time into writeEachLine().
*/
func writeToFileHandle(file io.Writer, records [][]string) {
	for _, rec := range records {
		writeEachLine(file, rec)
	}
}

/*
writeEachLine() writes a tab-separated line based on elements in []string.
*/
func writeEachLine(file io.Writer, rec []string) {
	var err error
	for i := 0; i < len(rec); i++ {
		if i == len(rec)-1 { // Ending of each line
			_, err = fmt.Fprintf(file, "%s\n", rec[i])
			exception.PanicOnErr(err)
		} else {
			_, err = fmt.Fprintf(file, "%s\t", rec[i])
			exception.PanicOnErr(err)
		}
	}
}

/*
writeMeta() writes metadata based on the []string of metadata
Separate each entry with ":" instead of a tab.
*/
func writeMeta(file io.Writer, rec []string) {
	var err error
	for i := 0; i < len(rec); i++ {
		if i == len(rec)-1 {
			_, err = fmt.Fprintf(file, "%s\n", rec[i])
			exception.PanicOnErr(err)
		} else {
			_, err = fmt.Fprintf(file, "%s:", rec[i])
			exception.PanicOnErr(err)
		}
	}
}

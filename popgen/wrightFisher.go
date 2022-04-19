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
	PopSize    int
	MutRate    float64
	NumGen     int
	GenomeSize int
	RFitness   float64
	GcContent  float64
	SetSeed    int64
	Verbose    bool
	Fasta      bool
	Vcf        bool
}

type WrightFisherPopData struct {
	Fasta     []fasta.Fasta
	Vcf       vcf.Vcf
	Meta      []string
	Freq      [][][]float64
	Ancestral []string //Ancestral states of each allele on each site (aka replicates)
	Settings  WrightFisherSettings
}

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

	writeMeta(file, wf.Meta)
	writeEachLine(file, header)

	writeToFileHandle(file, floatArrayToString(wf))

	exception.PanicOnErr(file.Close())

	//fmt.Printf("%s\n", floatArrayToString(wf.AFreq))
}

func floatArrayToString(wf WrightFisherPopData) [][]string {
	nrow := (wf.Settings.NumGen + 1) * wf.Settings.GenomeSize
	ncol := 7
	str := make([][]string, nrow)

	var t, s, r, b int // iterator for generation, site, row, and base

	for t = 0; t <= wf.Settings.NumGen; t++ {
		for s = 0; s < wf.Settings.GenomeSize; s++ {
			r = (t * wf.Settings.GenomeSize) + s
			str[r] = make([]string, ncol)
			str[r][0] = strconv.Itoa(t)
			str[r][1] = strconv.Itoa(s)
			for b = 0; b < 4; b++ {
				str[r][2+b] = strconv.FormatFloat(wf.Freq[t][s][b], 'f', 5, 64)
			}
			str[r][6] = wf.Ancestral[s]
		}
	}
	return str
}

func writeToFileHandle(file io.Writer, records [][]string) {
	for _, rec := range records {
		writeEachLine(file, rec)
	}
}

func writeEachLine(file io.Writer, rec []string) {
	var err error
	for i := 0; i < len(rec); i++ {
		if i == len(rec)-1 {
			_, err = fmt.Fprintf(file, "%s\n", rec[i])
			exception.PanicOnErr(err)
		} else {
			_, err = fmt.Fprintf(file, "%s\t", rec[i])
			exception.PanicOnErr(err)
		}

	}

}

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

// Command Group: "Variant Calling & Annotation"

package main

import (
	"bytes"
	"encoding/json"
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/vcf"
	"io"
	"net/http"
	"strings"
)

func usage() {
	fmt.Print(
		"vcfWebAnnotate - Annotate a vcf file by querying various databases via CellBase.\n\n" +
			"Usage:\n" +
			"  vcfWebAnnotate [options] in.vcf\n\n" +
			"Options:\n\n")
	flag.PrintDefaults()
}

// Example of JSON layout can be found in the link below
// http://bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/v4/hsapiens/genomic/variant/chr1%3A878884%3AC%3AT,chr1%3A878917%3AT%3AA,chr1%3A878920%3AT%3AA,chr1%3A878991%3AGTGTT%3AG,chr1%3A879229%3AA%3AT,chr1%3A879231%3AA%3AC,chr1%3A879897%3AT%3AC,chr1%3A879957%3AG%3AT/annotation?assembly=grch38

type Responses struct { // the response for each variant queried
	Responses []Response `json:"response"`
}

type Response struct { // all the results for a single variant (e.g. multiple transcripts)
	Results []Result `json:"result"`
}

type Result struct { // data on a particular variant on a particular transcript
	Chr             string          `json:"chromosome"`
	Start           int             `json:"start"`
	Ref             string          `json:"reference"`
	Alt             string          `json:"alternate"`
	SnpId           string          `json:"id"`
	ConsequenceType string          `json:"displayConsequenceType"`
	Consequences    []Consequence   `json:"consequenceTypes"`
	PopAlleleFreqs  []PopAlleleFreq `json:"populationFrequencies"`
	// Fields below exist but are not collected for this cmd
	// conservation
	// geneExpression
	// geneTraitAssociation
	// geneDrugInteraction
	// cytoband
	// repeat
}

type PopAlleleFreq struct { // pop allele frequency from gnomAD and 1kgp
	Study      string  `json:"study"`
	Population string  `json:"population"`
	RefAf      float64 `json:"refAlleleFreq"`
	AltAf      float64 `json:"altAlleleFreq"`
}

type Consequence struct { // information about variant effect
	GeneName           string            `json:"geneName"`
	GeneId             string            `json:"ensemblGeneId"`
	TranscriptId       string            `json:"ensemblTranscriptId"`
	Strand             string            `json:"strand"`
	Biotype            string            `json:"biotype"`
	ProteinAnnotations ProteinAnnotation `json:"proteinVariantAnnotation"`
}

type ProteinAnnotation struct { // variants effect on protein
	Pos                int                  `json:"position"`
	Ref                string               `json:"reference"`
	Alt                string               `json:"alternate"`
	SubstitutionScores []SubstitutionScores `json:"substitutionScores"`
}

type SubstitutionScores struct { // sift & polyphen scores
	Source      string  `json:"source"`
	Score       float64 `json:"score"`
	Description string  `json:"description"`
}

func queryWorker(filledBufChan <-chan []vcf.Vcf, emptyBufChan chan<- []vcf.Vcf, outfile io.Writer) {
	baseUrl := "http://bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/v4/hsapiens/genomic/variant/"
	queryUrl := new(strings.Builder)
	var responses Responses
	data := new(bytes.Buffer)

	for buf := range filledBufChan { // get a slice of vcfs to query
		queryUrl.Reset()
		queryUrl.WriteString(baseUrl) // start building query url
		for i := range buf {          // generates a comma separated list of variants in the url
			if i > 0 {
				queryUrl.WriteByte(',')
			}
			queryUrl.WriteString(fmt.Sprintf("%s%%3A%d%%3A%s%%3A%s", buf[i].Chr, buf[i].Pos, buf[i].Ref, buf[i].Alt[0]))
		}
		queryUrl.WriteString("/annotation?assembly=grch38")
		response, err := http.Get(queryUrl.String()) // query
		exception.PanicOnErr(err)

		data.Reset()
		_, err = data.ReadFrom(response.Body)
		exception.PanicOnErr(err)
		err = json.Unmarshal(data.Bytes(), &responses)
		exception.PanicOnErr(err)

		annotateVcfs(buf, responses)
		vcf.WriteVcfToFileHandle(outfile, buf)
		emptyBufChan <- buf // return buffer for reuse
	}
	close(emptyBufChan)
}

func annotateVcfs(vcfs []vcf.Vcf, ann Responses) {

}

func vcfWebAnnotate(data <-chan vcf.Vcf, header vcf.Header, outfile io.Writer, batchSize int, numBuffers int) {
	vcf.NewWriteHeader(outfile, header)
	filledBufChan := make(chan []vcf.Vcf, numBuffers)
	emptyBufChan := make(chan []vcf.Vcf, numBuffers)

	for i := 0; i < numBuffers-1; i++ {
		emptyBufChan <- make([]vcf.Vcf, 0, batchSize) // send all but 1 buffer to empty
	}
	buf := make([]vcf.Vcf, 0, batchSize)

	go queryWorker(filledBufChan, emptyBufChan, outfile)

	for v := range data { // read vcfs until you have a full batch then send for annotation
		if len(buf) == batchSize {
			filledBufChan <- buf
			buf = <-emptyBufChan
			buf = buf[:0]
		}
		buf = append(buf, v)
	}

	if len(buf) > 0 { // send any variants remaining in the buffer
		filledBufChan <- buf
	}
	close(filledBufChan)

	for _ = range emptyBufChan { // stall until queryWorker is finished
	}
}

func main() {
	var outfile *string = flag.String("o", "stdout", "output to vcf file")
	var batchSize *int= flag.Int("batchSize", 200, "number of variants to pool before querying web")
	var numBuffer *int= flag.Int("bufferSize", 2, "number of batchSize buffers to keep in memory")
	//TODO species
	//TODO assembly
	//TODO desired annotation fields
	flag.Parse()
	flag.Usage = usage

	var infile string = flag.Arg(0)
	if infile == "" {
		usage()
		return
	}

	vcfs, header := vcf.GoReadToChan(infile)
	out := fileio.EasyCreate(*outfile)
	vcfWebAnnotate(vcfs, header, out, *batchSize, *numBuffer)
	err := out.Close()
	exception.PanicOnErr(err)
}

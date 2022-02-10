// Command Group: "Variant Calling & Annotation"

package main

import (
	"bytes"
	"encoding/json"
	"errors"
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/vcf"
	"io"
	"log"
	"net/http"
	"os"
	"strings"
)

var errNoData error = errors.New("NA")

func usage() {
	fmt.Print(
		"vcfWebAnnotate - Annotate a vcf file by querying various databases via CellBase.\n" +
			"Note that this tool currently only works with VCFs for human hg38 (GRCh38).\n" +
			"Currently only the annotation for the first overlapping transcript is reported\n\n" +
			"Usage:\n" +
			"  vcfWebAnnotate [options] in.vcf\n\n" +
			"Options:\n\n")
	flag.PrintDefaults()
}

func vcfWebAnnotate(data <-chan vcf.Vcf, header vcf.Header, outfile io.Writer, batchSize int, numBuffers int, resume bool) {
	if !resume {
		header = addAnnotationHeader(header)
		vcf.NewWriteHeader(outfile, header)
	}
	filledBufChan := make(chan []vcf.Vcf, numBuffers)
	emptyBufChan := make(chan []vcf.Vcf, numBuffers)

	for i := 0; i < numBuffers-1; i++ {
		emptyBufChan <- make([]vcf.Vcf, 0, batchSize) // send all but 1 buffer to empty
	}
	buf := make([]vcf.Vcf, 0, batchSize)

	// queryWorker reads vcfs from filledBufChan and returns the slice to emptyBufChan
	// to be reused. queryWorker also currently handles annotation and writing.
	// TODO delegate annotation and writing to maximize query throughput
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

	for range emptyBufChan { // stall until queryWorker is finished
	}
}

func queryWorker(filledBufChan <-chan []vcf.Vcf, emptyBufChan chan<- []vcf.Vcf, outfile io.Writer) {
	baseUrl := "http://bioinfo.hpc.cam.ac.uk/cellbase/webservices/rest/v4/hsapiens/genomic/variant/annotation?assembly=grch38"
	query := new(bytes.Buffer)
	data := new(bytes.Buffer)

	for buf := range filledBufChan { // get a slice of vcfs to query
		var responses Responses
		query.Reset()
		for i := range buf { // generates a comma separated list of variants in the url
			if i > 0 {
				query.WriteByte(',')
			}
			query.WriteString(fmt.Sprintf("%s:%d:%s:%s", buf[i].Chr, buf[i].Pos, buf[i].Ref, buf[i].Alt[0]))
		}
		response, err := http.Post(baseUrl, "text/plain", query) // query
		exception.PanicOnErr(err)
		if response.StatusCode != 200 { // code 200 is a successful request
			log.Fatal(response.Status) // kill program and report failure code
		}

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

func annotateVcfs(vcfs []vcf.Vcf, res Responses) {
	var ann []string
	var consequence Consequence
	var proteinAnn ProteinAnnotation
	var maxAf float64
	var err error
	for i := range vcfs {
		ann = ann[:0] // reset

		maxAf, err = getMaxPopAf(res.Responses[i])
		if err != errNoData {
			ann = append(ann, fmt.Sprintf("MaxPopAF=%.2g", maxAf))
		}

		if len(res.Responses[i].Results[0].Consequences) == 0 {
			continue
		}
		consequence = res.Responses[i].Results[0].Consequences[0]

		if res.Responses[i].Results[0].ConsequenceType != "" {
			ann = append(ann, fmt.Sprintf("Consequence=%s", res.Responses[i].Results[0].ConsequenceType))
		}

		if consequence.GeneName != "" {
			ann = append(ann, fmt.Sprintf("Gene=%s", consequence.GeneName))
		}

		if consequence.TranscriptId != "" {
			ann = append(ann, fmt.Sprintf("Transcript=%s", consequence.TranscriptId))
		}

		if consequence.ProteinAnnotation.Ref != "" {
			proteinAnn = consequence.ProteinAnnotation
			ann = append(ann, fmt.Sprintf("ProteinEffect=%s",
				fmt.Sprintf("%s%d%s", proteinAnn.Ref, proteinAnn.Pos, proteinAnn.Alt)))
		}

		if vcfs[i].Info == "." {
			vcfs[i].Info = strings.Join(ann, ";")
		} else {
			vcfs[i].Info += ";" + strings.Join(ann, ";")
		}
	}
}

func getMaxPopAf(r Response) (float64, error) {
	var maxAf float64 = -1
	for _, p := range r.Results[0].PopAlleleFreqs {
		if p.Study == "" {
			return -1, errNoData
		}
		if p.AltAf > maxAf {
			maxAf = p.AltAf
		}
	}
	if maxAf == -1 {
		return -1, errNoData
	}
	return maxAf, nil
}

func addAnnotationHeader(header vcf.Header) vcf.Header {
	var insertLocation int
	for insertLocation = range header.Text {
		if strings.HasPrefix(header.Text[insertLocation], "##contig") {
			break
		}
	}

	savedHeader := make([]string, len(header.Text[insertLocation:]))
	copy(savedHeader, header.Text[insertLocation:])
	header.Text = header.Text[:insertLocation]

	// MaxPopAF
	header.Text = append(header.Text,
		"##INFO=<ID=MaxPopAF,Number=1,Type=Float,Description=\"Maximum allele frequency of any population in CellBase\",Source=\"bioinfo.hpc.cam.ac.uk/cellbase/webservices\",Version=\"v4\">")
	// Consequence
	header.Text = append(header.Text,
		"##INFO=<ID=Consequence,Number=1,Type=String,Description=\"Variant consequence\",Source=\"bioinfo.hpc.cam.ac.uk/cellbase/webservices\",Version=\"v4\">")
	// Gene
	header.Text = append(header.Text,
		"##INFO=<ID=Gene,Number=1,Type=String,Description=\"Nearest gene\",Source=\"bioinfo.hpc.cam.ac.uk/cellbase/webservices\",Version=\"v4\">")
	// Transcript
	header.Text = append(header.Text,
		"##INFO=<ID=Transcript,Number=1,Type=String,Description=\"Ensembl transcript id\",Source=\"bioinfo.hpc.cam.ac.uk/cellbase/webservices\",Version=\"v4\">")
	// ProteinEffect
	header.Text = append(header.Text,
		"##INFO=<ID=ProteinEffect,Number=1,Type=String,Description=\"Effect of variant on protein\",Source=\"bioinfo.hpc.cam.ac.uk/cellbase/webservices\",Version=\"v4\">")

	header.Text = append(header.Text, savedHeader...)
	return header
}

func burnRecords(input, partial <-chan vcf.Vcf) {
	var lastProcessedRecord vcf.Vcf
	for lastProcessedRecord = range partial {
	} // read through partially written file to get last record
	for v := range input {
		if v.Pos == lastProcessedRecord.Pos {
			return
		}
	}
}

func main() {
	var outfile *string = flag.String("o", "stdout", "output to vcf file")
	var batchSize *int = flag.Int("batchSize", 1000, "number of variants to pool before querying web")
	var numBuffer *int = flag.Int("bufferSize", 2, "number of batchSize buffers to keep in memory")
	var resume *bool = flag.Bool("resume", false, "resume a partially completed annotation. -o must be specified.")
	//TODO species
	//TODO assembly
	//TODO desired annotation fields
	//TODO data for all transcripts
	flag.Parse()
	flag.Usage = usage

	var infile string = flag.Arg(0)
	if infile == "" {
		usage()
		return
	}

	if *resume && *outfile == "stdout" {
		usage()
		log.Fatal("ERROR: output file (-o) must be specified to use -resume")
	}

	var out io.WriteCloser
	var err error
	vcfs, header := vcf.GoReadToChan(infile)

	if !*resume {
		out = fileio.EasyCreate(*outfile)
	} else {
		_, err = os.Stat(*outfile)
		if err != nil {
			log.Fatalf("ERROR: could not open %s. %s must exist to use -resume.", *outfile, *outfile)
		}
		partial, _ := vcf.GoReadToChan(*outfile)
		burnRecords(vcfs, partial)
		out, err = os.OpenFile(*outfile, os.O_APPEND|os.O_WRONLY, os.ModeAppend)
		exception.FatalOnErr(err)
	}

	defer out.Close() // ensure file closure, even on panic. double close is ok
	vcfWebAnnotate(vcfs, header, out, *batchSize, *numBuffer, *resume)
	err = out.Close()
	exception.PanicOnErr(err)
}

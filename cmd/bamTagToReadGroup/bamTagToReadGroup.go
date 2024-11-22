// Command Group: "SAM Tools"

// Count bases from sequencing data
package main

import (
	"flag"
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	"io"
	"log"
	"strings"
	"sync"
)

func usage() {
	fmt.Print(
		"bamTagToReadGroup - Add bam read groups based on read tags\n\n" +
			"Usage:\n" +
			"  bamTagToReadGroup [options] -tagId CB -tagValues barcodes.tsv -i input.bam > output.bam\n\n" +
			"Options:\n")
	flag.PrintDefaults()
}

func main() {
	var infile *string = flag.String("i", "", "Input bam file.")
	var outfile *string = flag.String("o", "stdout", "Output bam file.")
	var tagId *string = flag.String("tagId", "", "Bam tag ID that will be used to form read groups.")
	var tagValuesFile *string = flag.String("tagValues", "", "File of tag values that will be added to read groups, one value per line. "+
		"If bam record has a tagId value matching an entry in this file, the record will be annotated as belonging to a read group with the same name as the tag value. "+
		"If a bam record does not have a tag with tagId, or the tagId value does not match an entry in this file then the bam record will not be assigned a read group."+
		"NOTE: All existing read groups will be removed from the file.")
	flag.Parse()
	flag.Usage = usage

	if *infile == "" {
		usage()
		log.Fatalln("ERROR: must provide a bam file with -i")
	}

	if *tagId == "" || len(*tagId) != 2 {
		usage()
		log.Fatalln("Error: -tagId must be 2 characters")
	}

	if *tagValuesFile == "" {
		usage()
		log.Fatalln("ERROR: must provide a tag values file with -tagValues")
	}

	bamTagToReadGroup(*infile, *outfile, *tagId, *tagValuesFile)
}

func bamTagToReadGroup(infile, outfile, tagId, tagValuesFile string) {
	tagValues := fileio.Read(tagValuesFile)
	tagValueMap := make(map[string]bool)
	for i := range tagValues {
		tagValueMap[tagValues[i]] = true
	}

	inputRecordChan, header := sam.GoReadToChan(infile)
	addTagsToHeader(&header, tagValues)

	out := fileio.EasyCreate(outfile)
	defer safeclose(out)
	bw := sam.NewBamWriter(out, header)
	defer safeclose(bw)

	// spawn goroutine to write records after being updated
	writeChan := make(chan sam.Sam, 1000)
	wg := new(sync.WaitGroup)
	wg.Add(1)
	go func(writeChan <-chan sam.Sam, bw *sam.BamWriter, wg *sync.WaitGroup) {
		for rec := range writeChan {
			sam.WriteToBamFileHandle(bw, rec, 0)
		}
		wg.Done()
	}(writeChan, bw, wg)

	// read records and update read groups as we go
	for rec := range inputRecordChan {
		updateRecord(&rec, tagId, tagValueMap)
		writeChan <- rec
	}
	close(writeChan)

	// wait for goroutine to finish writing before return on main thread
	wg.Wait()
}

func safeclose(f io.Closer) {
	exception.PanicOnErr(f.Close())
	exception.PanicOnErr(err)
}

// addTagsToHeader adds each element in tagValues as a read group in the header
// format is `@RG	ID:$TAGVALUE	SM:$TAGVALUE	LB:$TAGVALUE`
func addTagsToHeader(header *sam.Header, tagValues []string) {
	newHeaderText := make([]string, 0, len(header.Text)+len(tagValues))
	for i := range header.Text {
		if !strings.HasPrefix(header.Text[i], "@RG") {
			newHeaderText = append(newHeaderText, header.Text[i])
		}
	}
	for i := range tagValues {
		newHeaderText = append(newHeaderText, fmt.Sprintf("@RG\tID:%s\tSM:%s\tLB:%s", tagValues[i], tagValues[i], tagValues[i]))
	}
	header.Text = newHeaderText
}

func updateRecord(rec *sam.Sam, tagId string, tagValueMap map[string]bool) {
	var err error

	// parse the extra field where tags are since we are coming from bam
	err = sam.ParseExtra(rec)
	exception.PanicOnErr(err)

	splitExtra := strings.Split(rec.Extra, "\t")
	var newExtra []string = make([]string, 0, len(splitExtra)+1)

	var tagValue string

	// first remove any exiting read group and pull the tag value at the same time
	for i := range splitExtra {
		if !strings.HasPrefix(splitExtra[i], "RG:Z:") {
			newExtra = append(newExtra, splitExtra[i])
		}

		if strings.HasPrefix(splitExtra[i], tagId+":") {
			tagValue = splitExtra[i][5:] // skips characters defining tag (e.g. CB:Z:TagValue)
		}
	}

	// add read group if tag value is present in input values file
	if tagValueMap[tagValue] {
		newExtra = append(newExtra, fmt.Sprintf("RG:Z:%s", tagValue))
	}
	rec.Extra = strings.Join(newExtra, "\t")
}

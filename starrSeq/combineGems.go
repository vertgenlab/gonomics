package starrSeq

import (
	"fmt"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/sam"
	"strings"
)

// DetectMultipleGems takes in the input filename string from cellrangerBam and returns a bool of whether multiple GEM wells were detected and the input file name, updated with a tmp file name if
// > 1 input files were detected
func DetectMultipleGems(fileName string, cellTypeAnalysis string) (string, bool) {
	var bw *sam.BamWriter
	var combineBamWriter *fileio.EasyWriter
	var err error
	var tmpFileName string
	inFiles := strings.Split(fileName, ",")
	if len(inFiles) > 1 {
		var path, tmpFileSlice []string
		if cellTypeAnalysis != "" {
			fmt.Println("*** WARNING *** You are using multiple GEM wells with -cellTypeAnalysis. Using multiple GEM wells will add an additional suffix to cell barcodes to reinforce cell barcode uniqueness. " +
				"The first bam file provided will have the '_1' suffix, the second bam file '_2' and so on. This mirrors the default behavior of both Seurat merge() and integrate() functions. " +
				"If multiple GEM wells haven't be processed in the same way in Seurat and in the scStarrSeqAnalysis programs, cell lookup for -cellTypeAnalysis will be impaired.")
		}
		for _, i := range inFiles {
			path = strings.Split(i, "/")
			tmpFileSlice = append(tmpFileSlice, path[len(path)-1])
		}
		tmpFileSlice = append(tmpFileSlice, ".bam")
		tmpFileName = strings.Join(tmpFileSlice, "_")
		combineBamWriter = fileio.EasyCreate(tmpFileName)
		p := 1
		for _, i := range inFiles {
			bw = combineBams(i, combineBamWriter, &p, bw)
		}
		err = bw.Close()
		exception.PanicOnErr(err)
		err = combineBamWriter.Close()
		exception.PanicOnErr(err)
		return tmpFileName, true
	} else {
		return fileName, false
	}
}

// DetectMultipleGfpGems takes in the input GFP file name string from cellrangerBam and returns a bool of whether multiple GEM wells were detected and the input GFP file name, updated with
// a tmp file name if > 1 input files were detected
func DetectMultipleGfpGems(filename string) (string, bool) {
	var gfpBams []string
	var bwGFP *sam.BamWriter
	var combineBamWriterGFP *fileio.EasyWriter
	var err error
	var tmpFileName string
	gfpBams = strings.Split(filename, ",")
	if len(gfpBams) > 1 {
		var path, tmpFileSlice []string
		fmt.Println("*** WARNING *** You are using multiple GEM wells with -gfpNorm. Using multiple GEM wells will add an additional suffix to cell barcodes to reinforce cell barcode uniqueness. " +
			"The first bam file provided will have the '_1' suffix, the second bam file '_2' and so on. This mirrors the default behavior of both Seurat merge() and integrate() functions. " +
			"If multiple GEM wells haven't be processed in the same way in Seurat and in the scStarrSeqAnalysis programs, cell lookup for -cellTypeAnalysis will be impaired.")
		for _, i := range gfpBams {
			path = strings.Split(i, "/")
			tmpFileSlice = append(tmpFileSlice, path[len(path)-1])
		}
		tmpFileSlice = append(tmpFileSlice, ".bam")
		tmpFileName = strings.Join(tmpFileSlice, "_")
		combineBamWriterGFP = fileio.EasyCreate(tmpFileName)
		p := 1
		for _, i := range gfpBams {
			bwGFP = combineBams(i, combineBamWriterGFP, &p, bwGFP)
		}
		err = bwGFP.Close()
		exception.PanicOnErr(err)
		err = combineBamWriterGFP.Close()
		exception.PanicOnErr(err)
		return tmpFileName, true
	} else {
		return filename, false
	}
}

// combineBams is triggered when there are more than 1 input bams in a comma delimited list. This function appends a "_int" to the back of the cell barcode correstponding to the input file order
// and then concatenates all bam files together
func combineBams(a string, out *fileio.EasyWriter, iteration *int, bw *sam.BamWriter) *sam.BamWriter {
	var err error
	var columns, columns2 []string
	var j, join string
	var n int

	inChan, head := sam.GoReadToChan(a)
	if *iteration == 1 {
		bw = sam.NewBamWriter(out, head)
	}
	for i := range inChan {
		err = sam.ParseExtra(&i)
		columns = strings.Split(i.Extra, "\t")
		for n, j = range columns {
			columns2 = strings.Split(j, ":")
			if columns2[0] == "CB" {
				columns[n] = fmt.Sprintf("%s_%d", j, *iteration)
				join = strings.Join(columns, "\t")
				i.Extra = join
			}
		}
		sam.WriteToBamFileHandle(bw, i, 0)
	}
	exception.PanicOnErr(err)
	*iteration++
	return bw
}

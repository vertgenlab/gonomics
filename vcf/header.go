package vcf

import (
	"fmt"
	"io"
	"log"
	"strconv"
	"strings"

	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fileio"
)

// Header contains all of the information present in the header section of a VCF.
// Info, Filter, Format, and Contig lines are parsed into maps keyed by ID.
type Header struct {
	FileFormat string                         // ##fileformat=VCFv4.3
	Info       map[string]InfoHeader          // key=ID
	Filter     map[string]FilterHeader        // key=ID
	Format     map[string]FormatHeader        // key=ID
	Chroms     map[string]chromInfo.ChromInfo // key=chrom name
	Samples    map[string]int                 // key=samplename val=index in Sample
	Text       []string                       // raw text
}

// InfoType stores the type of variable that a field in the Header holds.
type InfoType byte

func (t InfoType) String() string {
	switch t {
	case Integer:
		return "Integer"
	case Float:
		return "Float"
	case Flag:
		return "Flag"
	case String:
		return "String"
	case Character:
		return "Character"
	default:
		log.Panicln("unknown type")
		return ""
	}
}

const (
	Integer InfoType = iota
	Float
	Flag
	Character
	String
)

// Key is the identifying information for a given info field.
// (e.g. the genotype field in format would be {"GT", "1", Integer, true}.)
// (e.g. a read counter may be {"ReadCount", "R", Integer, true}.)
type Key struct {
	Id       string
	Number   string // numeral or 'A', 'G', 'R', '.'
	DataType InfoType
	IsFormat bool // true if this key is for a Format field, false for Info
}

// InfoHeader contains info encoded by header lines beginning in ##INFO.
type InfoHeader struct {
	Key
	Description string
	Source      string
	Version     string
}

// FilterHeader contains info encoded by header lines beginning in ##FILTER.
type FilterHeader struct {
	Id          string
	Description string
}

// FormatHeader contains info encoded by header lines beginning in ##FORMAT.
type FormatHeader struct {
	Key
	Description string
}

// ReadHeader reads and parses the header in a vcf file.
func ReadHeader(er *fileio.EasyReader) Header {
	var line string
	var err error
	var nextBytes []byte
	var headerText []string
	for nextBytes, err = er.Peek(1); err == nil && nextBytes[0] == '#'; nextBytes, err = er.Peek(1) {
		line, _ = fileio.EasyNextLine(er)
		headerText = append(headerText, line)
	}
	return parseHeader(headerText)
}

// SampleNamesInOrder takes in the header and gives back the sample names in the order in which they appear in the header
func SampleNamesInOrder(header Header) []string {
	var answer []string = make([]string, len(header.Samples))
	for sampleName, idx := range header.Samples {
		answer[idx] = sampleName
	}
	return answer
}

// newHeader allocates memory for a new vcf Header.
func newHeader() Header {
	var h Header
	h.Info = make(map[string]InfoHeader)
	h.Filter = make(map[string]FilterHeader)
	h.Format = make(map[string]FormatHeader)
	h.Chroms = make(map[string]chromInfo.ChromInfo)
	h.Samples = make(map[string]int)
	return h
}

// parseHeader parses the text in a vcf Header to a more useful format.
func parseHeader(text []string) Header {
	if len(text) == 0 { // vcf w/o header is valid
		return Header{}
	}

	header := newHeader()
	header.Text = text
	var chromIdx int
	for i := range text {
		switch strings.Split(text[i][2:], "=")[0] { // pulls line label (e.g. "INFO" from "##INFO=<...>")
		case "fileformat":
			header.FileFormat = parseFileFormatFromHeader(text[i])

		case "contig":
			parseChromsFromHeader(text[i], header.Chroms, chromIdx)
			chromIdx++

		case "INFO":
			parseInfoFromHeader(text[i], header.Info)

		case "FILTER":
			parseFilterFromHeader(text[i], header.Filter)

		case "FORMAT":
			parseFormatFromHeader(text[i], header.Format)
		}
	}
	header.Samples = parseSamplesFromHeader(text[len(text)-1])
	return header
}

// parseFileFormatFromHeader parses a line beginning with ##fileformat to determine the vcf version.
func parseFileFormatFromHeader(line string) string {
	return strings.Split(line, "=")[1]
}

// parseChromsFromHeader parses a line beginning with ##contig to a map keying the chrom name to the chromInfo.
func parseChromsFromHeader(line string, chroms map[string]chromInfo.ChromInfo, chromIdx int) {
	fields := getHeaderFields(line)
	var chrom chromInfo.ChromInfo
	chrom.Order = chromIdx
	var err error

	// Parse chromInfo
	for i := range fields {
		switch {
		case strings.HasPrefix(fields[i], "ID="):
			chrom.Name = fields[i][3:] // cut off "ID="

		case strings.HasPrefix(fields[i], "length="):
			chrom.Size, err = strconv.Atoi(fields[i][7:]) // cut off "length=" and parse
			if err != nil {
				log.Panicf("trouble parsing length field in vcf header.\n Tried to parse '%s'", fields[i][7:])
			}
		}
	}

	// Add chrom to map
	_, alreadySawChrom := chroms[chrom.Name]
	if alreadySawChrom {
		log.Fatalf("ERROR: contig names in header must be unique. Saw %s multiple times", chrom.Name)
	}
	chroms[chrom.Name] = chrom
}

// parseInfoFromHeader parses a line beginning with ##INFO to a map keying the ID to the info fields information.
func parseInfoFromHeader(line string, info map[string]InfoHeader) {
	var answer InfoHeader
	answer.Id, answer.Number, answer.DataType, answer.Description, answer.Source, answer.Version = parseHeaderFields(line)

	_, duplicateId := info[answer.Id]
	if duplicateId {
		log.Fatalf("duplicate ID in info header: '%s'", answer.Id)
	}
	answer.IsFormat = false
	info[answer.Id] = answer
}

// parseFilterFromHeader parses a line beginning with ##FILTER to a map keying the ID to the filter information.
func parseFilterFromHeader(line string, filter map[string]FilterHeader) {
	var answer FilterHeader
	answer.Id, _, _, answer.Description, _, _ = parseHeaderFields(line)

	_, duplicateId := filter[answer.Id]
	if duplicateId {
		log.Fatalf("duplicate ID in filter header: '%s'", answer.Id)
	}

	filter[answer.Id] = answer
}

// parseFormatFromHeader parses a line beginning with ##FORMAT to a map keying the ID to the format information.
func parseFormatFromHeader(line string, formatData map[string]FormatHeader) {
	var answer FormatHeader
	answer.Id, answer.Number, answer.DataType, answer.Description, _, _ = parseHeaderFields(line)

	_, duplicateId := formatData[answer.Id]
	if duplicateId {
		log.Fatalf("duplicate ID in answer header: '%s'", answer.Id)
	}

	answer.IsFormat = true
	formatData[answer.Id] = answer
}

// parseSamplesFromHeader reads samples from the column names in the final line of the header. Returns a
// map keying sample name to index in the Samples field of the vcf.
func parseSamplesFromHeader(line string) map[string]int {
	colNames := strings.Split(line, "\t")
	if len(colNames) < 10 { // no samples in file (may be sites only vcf)
		return nil
	}
	if colNames[0] != "#CHROM" {
		log.Fatalf("ERROR: malformed header. Expected final header line to begin with '#CHROM'. Actually began with '%s'.\n", colNames[0])
	}
	sampleMap := make(map[string]int)
	samples := colNames[9:]
	var nameExists bool
	for i := range samples {
		_, nameExists = sampleMap[samples[i]]
		if nameExists {
			log.Fatalf("ERROR: cannot have duplicate sample names. Sample %s was present more than once.\n", samples[i])
		}
		sampleMap[samples[i]] = i
	}
	return sampleMap
}

// getHeaderFields parses the comma delimited fields within the '<' '>' delimited portion of a header line.
// e.g. ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype"> returns []string{ID=GT Number=1 Type=String Description="Genotype"}
func getHeaderFields(line string) []string {
	start := strings.IndexRune(line, '<')
	if line[len(line)-1] != '>' || start == -1 {
		log.Fatalf("ERROR: malformed header line (must have a field delimited by '<' and '>')\n%s\n", line)
	}
	return strings.Split(line[start+1:len(line)-1], ",")
}

// parseHeaderFields returns the Id, Number, Type, Description, Source, and Version present
// in a header line (if field is applicable, returns zero value otherwise).
func parseHeaderFields(line string) (Id string, Number string, Type InfoType, Description string, Source string, Version string) {
	fields := getHeaderFields(line)
	for i := range fields {
		switch {
		case strings.HasPrefix(fields[i], "ID="):
			Id = fields[i][3:]

		case strings.HasPrefix(fields[i], "Number="):
			Number = fields[i][7:]

		case strings.HasPrefix(fields[i], "Type="):
			switch fields[i][5:] {
			case "Integer":
				Type = Integer
			case "Float":
				Type = Float
			case "Flag":
				Type = Flag
			case "Character":
				Type = Character
			case "String":
				Type = String
			default:
				log.Panicf("Unrecognized type in vcf header: %s", fields[i][5:])
			}

		case strings.HasPrefix(fields[i], "Description="):
			Description = strings.Trim(fields[i][12:], "\"")

		case strings.HasPrefix(fields[i], "Source="):
			Source = strings.Trim(fields[i][7:], "\"")

		case strings.HasPrefix(fields[i], "Version="):
			Version = strings.Trim(fields[i][8:], "\"")
		}
	}
	return
}

func processHeader(header Header, line string) Header {
	if strings.HasPrefix(line, "#") {
		header.Text = append(header.Text, line)
	} else {
		log.Fatal("There was an error reading the header line")
	}
	return header
}

//If you have multiple samples to add to the header, use strings.Join(samples[], "\t") as an argument that combines multiple samples by tabs
//Name input is strictly used to push in a name for the sample column.

func NewHeader(name string) Header {
	var header Header
	header.Text = append(header.Text, "##fileformat=VCFv4.2")
	header.Text = append(header.Text, fmt.Sprintf("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"))
	return header
}

//The other option is add everything as one string, so we don';t have to keep appending, parsing will in general be the same, but append is more consistent with how we read in line, by line.
/*
var text string =
		"##fileformat=VCFv4.2\n"+
		"##fileDate="+t.Format("20060102")+"\n"+
		"##source=github.com/vertgenlab/gonomics\n"+
		"##phasing=none\n"+
		"##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant: DEL, INS, DUP, INV, CNV, BND\">\n"+
		"##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant described in this record\">\n"+
		"##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n"+
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"+
		fmt.Sprintf("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s", name)

*/

func NewWriteHeader(file io.Writer, header Header) {
	var err error
	for h := 0; h < len(header.Text); h++ {
		_, err = fmt.Fprintf(file, "%s\n", header.Text[h])
		common.ExitIfError(err)
	}
}

func WriteMultiSamplesHeader(file io.Writer, header Header, listNames []string) {
	var err error
	for h := 0; h < len(header.Text); h++ {
		if strings.Contains(header.Text[h], "#CHROM\t") {
			name := strings.Join(listNames, "\t")
			_, err = fmt.Fprintf(file, fmt.Sprintf("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n", name))
			common.ExitIfError(err)
		} else {
			_, err = fmt.Fprintf(file, "%s\n", header.Text[h])
			common.ExitIfError(err)
		}
	}
}

//Uses Vcf header to create 2 hash maps 1) is the sample index that maps the which allele each sample has in Vcf 2) hash reference chromsome names to an index (used to build uint64 containing chromID and position)
func HeaderToMaps(header Header) *SampleHash {
	var name string
	var index, hapIdx int16
	var hash *SampleHash = &SampleHash{Fa: make(map[string]int16), GIndex: make(map[string]int16)}
	for _, line := range header.Text {
		if strings.HasPrefix(line, "##contig") {
			name = strings.Split(strings.Split(line, "=")[2], ",")[0]
			_, ok := hash.Fa[name]
			if !ok {
				hash.Fa[name] = index
				index++
			}
		} else if strings.HasPrefix(line, "#CHROM") {
			words := strings.Split(line, "\t")[9:]
			for hapIdx = 0; hapIdx < int16(len(words)); hapIdx++ {
				hash.GIndex[words[hapIdx]] = hapIdx
			}
		}
	}
	return hash
}

//HeaderGetSampleList returns an ordered list of the samples present in the header of a Vcf file. Useful when adding or removing samples from a VCF.
func HeaderGetSampleList(header Header) []string {
	var answer []string
	for _, line := range header.Text {
		if strings.HasPrefix(line, "#CHROM") {
			return strings.Split(line, "\t")[9:]
		}
	}
	log.Fatalf("No Sample info in VCF line, cannot parse sample names.")
	return answer
}

//HeaderUpdateSampleList can be provided with a new list of samples to update the sample list in a Header.
func HeaderUpdateSampleList(header Header, newSamples []string) Header {
	var line string
	for i := 0; i < len(header.Text); i++ {
		if strings.HasPrefix(header.Text[i], "#CHROM") {
			line = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
			for j := 0; j < len(newSamples); j++ {
				line += "\t" + newSamples[j]
			}
			header.Text[i] = line
		}
	}
	return header
}

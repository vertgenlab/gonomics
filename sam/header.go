package sam

import (
	"fmt"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"strings"
)

// Header encodes the header of a sam file as both the raw header (Text),
// and semi-parsed fields (Metadata and Chroms).
type Header struct {
	Text     []string
	Metadata Metadata
	Chroms   []chromInfo.ChromInfo // tags SQ - SN and SQ - LN
}

// Metadata stores semi-parsed header data with several explicitly
// parsed fields (e.g. Version) an the rest of the header encoded in
// the AllTags field as a TagMap.
type Metadata struct {
	Version   string      // tag HD - VN
	SortOrder []SortOrder // tag HD - SO (for len == 1) or HD - SS (for len > 1)
	Grouping  Grouping    // tag HD - GO
	AllTags   TagMap      // map of all tags present in file. Contains duplicates of parsed fields
	// TODO: Additional chromInfo tags
	// TODO: ReadGroup Parsing
	// TODO: Program parsing
	Comments []string // tag CO
}

// Tag is a 2 rune identifier of data encoded in a sam header.
// Tags occur in sets where each header line begins with a tag
// which is further split into subtags. e.g. @SQ SN:ref LN:45
type Tag [2]rune

// TagMap organizes all tags into a map where the line tag (e.g. SQ)
// stores a slice where each element in the slice corresponds to one
// line that has the keyed line tag. For instance, the line tag 'SQ'
// occurs once per chromosome. Using the key 'SQ' in the tag map
// gives each SQ line as a slice indexed by order of occurrence.
//
// Each element in the resulting struct is itself a map keyed by Tags
// present in the line. For example, the data stored in the 'SN' tag
// of the 5th 'SQ' line would be retrieved by TapMap[SQ][4][SN].
//
// For convenience, the TagMap also stores tags that are further parsed
// in other fields. e.g. the version number can be retrieved by calling
// either Metadata.Version or by Metadata.AllTags[HD][0][VN]
//
// Note that comment lines (line tag 'CO') are not stored in TagMap.
type TagMap map[Tag][]map[Tag]string

type SortOrder string

const (
	Unknown         SortOrder = "unknown"
	Unsorted        SortOrder = "unsorted"
	QueryName       SortOrder = "queryname"
	Coordinate      SortOrder = "coordinate"
	Lexicographical SortOrder = "lexicographical"
	Natural         SortOrder = "natural"
	Umi             SortOrder = "umi"
)

var sortOrderMap = map[string]SortOrder{
	"unknown":         Unknown,
	"unsorted":        Unsorted,
	"queryname":       QueryName,
	"coordinate":      Coordinate,
	"lexicographical": Lexicographical,
	"natural":         Natural,
	"umi":             Umi,
}

type Grouping string

const (
	None      Grouping = "none"
	Query     Grouping = "query"
	Reference Grouping = "reference"
)

var groupingMap = map[string]Grouping{
	"none":      None,
	"query":     Query,
	"reference": Reference,
}

func processHeaderLine(header Header, line string) {
	var chrCount int = 0

	header.Text = append(header.Text, line)
	if strings.HasPrefix(line, "@SQ") && strings.Contains(line, "SN:") && strings.Contains(line, "LN:") {
		curr := chromInfo.ChromInfo{Name: "", Size: 0, Order: chrCount}
		chrCount++
		words := strings.Fields(line)
		for i := 1; i < len(words); i++ {
			elements := strings.Split(words[i], ":")
			switch elements[0] {
			case "SN":
				curr.Name = elements[1]
			case "LN":
				curr.Size = common.StringToInt(elements[1])
			}
		}
		if curr.Name == "" || curr.Size == 0 {
			//	log.Fatal(fmt.Errorf("Thought I would get a name and non-zero size on this line: %s\n", line))
		}
		header.Chroms = append(header.Chroms, curr)
	}
}

func ReadHeader(er *fileio.EasyReader) Header {
	var line string
	var err error
	var nextBytes []byte
	var header Header

	for nextBytes, err = er.Peek(1); err == nil && nextBytes[0] == '@'; nextBytes, err = er.Peek(1) {
		line, _ = fileio.EasyNextLine(er)
		processHeaderLine(header, line)
	}
	return header
}

func WriteHeaderToFileHandle(file *fileio.EasyWriter, header Header) error {
	var err error

	for i := range header.Text {
		_, err = fmt.Fprintf(file, "%s\n", header.Text[i])
		common.ExitIfError(err)
	}
	return nil
}

func ChromInfoSamHeader(chromSize []chromInfo.ChromInfo) Header {
	var header Header
	header.Text = append(header.Text, "@HD\tVN:1.6\tSO:unsorted")
	var words string

	for i := 0; i < len(chromSize); i++ {
		words = fmt.Sprintf("@SQ\tSN:%s\tLN:%d", chromSize[i].Name, chromSize[i].Size)
		header.Text = append(header.Text, words)
	}
	return header
}

func ChromInfoMapSamHeader(chromSize map[string]chromInfo.ChromInfo) *Header {
	var header Header
	header.Text = append(header.Text, "@HD\tVN:1.6\tSO:unsorted")
	var words string
	var i int
	for i = 0; i < len(chromSize); {
		for j := range chromSize {
			if i == chromSize[j].Order {
				words = fmt.Sprintf("@SQ\tSN:%s\tLN:%d", chromSize[j].Name, chromSize[j].Size)
				header.Text = append(header.Text, words)
				i++
			}
		}
	}
	return &header
}

func FastaHeader(ref []fasta.Fasta) *Header {
	var header Header
	header.Text = append(header.Text, "@HD\tVN:1.6\tSO:unsorted")
	var words string

	for i := 0; i < len(ref); i++ {
		words = fmt.Sprintf("@SQ\tSN:%s\tLN:%d", ref[i].Name, len(ref[i].Seq))
		header.Text = append(header.Text, words)
		header.Chroms = append(header.Chroms, chromInfo.ChromInfo{Name: ref[i].Name, Size: len(ref[i].Seq)})
	}
	return &header
}

package sam

import (
	"bytes"
	"github.com/vertgenlab/gonomics/chromInfo"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
	"strconv"
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

// Tag is a 2 byte identifier of data encoded in a sam header.
// Tags occur in sets where each header line begins with a tag
// which is further split into subtags. e.g. @SQ SN:ref LN:45
type Tag [2]byte

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

// Sort order defines whether the file is sorted and if so, how it was sorted
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

// sortOrderMap provides easy lookup of string to SortOrder
var sortOrderMap = map[string]SortOrder{
	"unknown":         Unknown,
	"unsorted":        Unsorted,
	"queryname":       QueryName,
	"coordinate":      Coordinate,
	"lexicographical": Lexicographical,
	"natural":         Natural,
	"umi":             Umi,
}

// Grouping defines how the reads are grouped, if they are not sorted
type Grouping string

const (
	None      Grouping = "none"
	Query     Grouping = "query"
	Reference Grouping = "reference"
)

// groupingMap provides easy lookup of string to Grouping
var groupingMap = map[string]Grouping{
	"none":      None,
	"query":     Query,
	"reference": Reference,
}

// ReadHeaderBytes processes the contiguous header from a ByteReader
// and advances the Reader past the header lines.
func ReadHeaderBytes(br *fileio.ByteReader) Header {
	var answer Header
	var buff *bytes.Buffer
	var done bool
	for peek, err := br.Peek(1); err == nil && peek[0] == '@' && !done; peek, err = br.Peek(1) {
		buff, done = fileio.ReadLine(br)
		answer.Text = append(answer.Text, buff.String())
	}

	answer.Metadata.AllTags, answer.Metadata.Comments = parseTagsAndComments(answer.Text)
	answer.Chroms = getChromInfo(answer.Metadata.AllTags)
	answer.Metadata.Version = getVersion(answer.Metadata.AllTags)
	answer.Metadata.SortOrder = getSortOrder(answer.Metadata.AllTags)
	answer.Metadata.Grouping = getGrouping(answer.Metadata.AllTags)
	return answer
}

// parseTagsAndComments parses header text into a TagMap and a slice of comment lines
func parseTagsAndComments(text []string) (tags TagMap, comments []string) {
	var currTag Tag
	tags = make(TagMap)
	for _, line := range text {
		words := strings.Split(line, "\t")

		if words[0][0] != '@' || len(words[0]) != 3 {
			log.Fatalf("Error: malformed header line: %s", line)
		}

		if words[0] == "@CO" {
			comments = append(comments, line)
			continue
		}

		copy(currTag[:], words[0][1:]) // copy Tag

		tags[currTag] = append(tags[currTag], parseSubTags(words[1:]))
	}

	return
}

// parseSubTags parses a single line of tab delimited tagsets
func parseSubTags(tagsets []string) map[Tag]string {
	var currTag Tag
	var alreadyUsed bool
	answer := make(map[Tag]string)
	for _, tagset := range tagsets {
		if tagset[2] != ':' {
			log.Printf("Warning: ignoring malformed tag in header: %s\n", tagset)
			continue
		}
		copy(currTag[:], tagset[0:2]) // copy Tag

		if _, alreadyUsed = answer[currTag]; alreadyUsed {
			log.Fatalf("Error: same tag used multiple times per line: %s", tagset[0:2])
		}

		answer[currTag] = tagset[3:]
	}

	return answer
}

// getChromInfo further parses tags stored in a TagMap to extract ChromInfo
func getChromInfo(tags TagMap) []chromInfo.ChromInfo {
	chroms := tags[[2]byte{'S', 'Q'}]
	answer := make([]chromInfo.ChromInfo, len(chroms))
	for idx, chrom := range chroms {
		size, err := strconv.Atoi(chrom[[2]byte{'L', 'N'}])
		if err != nil {
			log.Fatalf("Error: could not convert chromosome size: \"%s\" to an integer in sam header", chrom[[2]byte{'L', 'N'}])
		}
		answer[idx] = chromInfo.ChromInfo{
			Name:  chrom[[2]byte{'S', 'N'}],
			Size:  size,
			Order: idx,
		}
	}
	return answer
}

// getVersion pulls the version number from a TagMap
func getVersion(tags TagMap) string {
	return tags[[2]byte{'H', 'D'}][0][[2]byte{'V', 'N'}]
}

// getSortOrder pulls the sort order from a TagMap
func getSortOrder(tags TagMap) []SortOrder {
	var answer []SortOrder
	order := tags[[2]byte{'H', 'D'}][0][[2]byte{'S', 'O'}]
	if !strings.Contains(order, ":") { // if only one sort order
		return append(answer, sortOrderMap[order])
	} else { // if multiple colon delimited sort orders
		words := strings.Split(order, ":")
		for _, word := range words {
			if _, ok := sortOrderMap[word]; ok {
				answer = append(answer, sortOrderMap[word])
			} else {
				log.Fatalf("Error: unrecognized sort order in sam header: %s", word)
			}
		}
	}
	return answer
}

// getGrouping pulls the grouping from a TagMap
func getGrouping(tags TagMap) Grouping {
	return groupingMap[tags[[2]byte{'H', 'D'}][0][[2]byte{'G', 'O'}]]
}

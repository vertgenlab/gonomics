package sam

import (
	"log"
	"strconv"
	"strings"

	"github.com/vertgenlab/gonomics/chromInfo"
)

// Metadata stores semi-parsed header data with several explicitly
// parsed fields (e.g. Version) and the rest of the header encoded in
// the AllTags field as a HeaderTagMap.
type Metadata struct {
	Version   string      // tag HD - VN
	SortOrder []SortOrder // tag HD - SO (for len == 1) or HD - SS (for len > 1)
	Grouping  Grouping    // tag HD - GO
	// TODO: Additional tag parsing (ReadGroup would be good)
	AllTags  HeaderTagMap // map of all tags present in file. Contains duplicates of parsed fields
	Comments []string     // tag CO
}

// Tag is a 2 byte identifier of data encoded in a sam header.
// Tags occur in sets where each header line begins with a tag
// which is further split into subtags. e.g. @SQ SN:ref LN:45.
type Tag [2]byte

// HeaderTagMap organizes all tags into a map where the line tag (e.g. SQ)
// stores a slice where each element in the slice corresponds to one
// line that has the keyed line tag. For instance, the line tag 'SQ'
// occurs once per chromosome. Using the key 'SQ' in the tag map
// gives each SQ line as a slice indexed by order of occurrence.
//
// Each element in the resulting slice is itself a map keyed by Tags
// present in the line. For example, the data stored in the 'SN' tag
// of the 5th 'SQ' line would be retrieved by TapMap[SQ][4][SN].
//
// For convenience, the HeaderTagMap also stores tags that are further parsed
// in other fields. e.g. the version number can be retrieved by calling
// either Metadata.Version or by Metadata.AllTags[HD][0][VN]
//
// Note that comment lines (line tag 'CO') are not stored in HeaderTagMap.
type HeaderTagMap map[Tag][]map[Tag]string

// Sort order defines whether the file is sorted and if so, how it was sorted.
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

// sortOrderMap provides easy lookup of string to SortOrder.
var sortOrderMap = map[string]SortOrder{
	"unknown":         Unknown,
	"unsorted":        Unsorted,
	"queryname":       QueryName,
	"coordinate":      Coordinate,
	"lexicographical": Lexicographical,
	"natural":         Natural,
	"umi":             Umi,
}

// Grouping defines how the reads are grouped, if they are not sorted.
type Grouping string

const (
	None      Grouping = "none"
	Query     Grouping = "query"
	Reference Grouping = "reference"
)

// groupingMap provides easy lookup of string to Grouping.
var groupingMap = map[string]Grouping{
	"none":      None,
	"query":     Query,
	"reference": Reference,
}

// ParseHeaderText parses the raw text of in the header struct and
// fills the appropriate fields in the rest of the header struct.
func ParseHeaderText(h Header) Header {
	if h.Text != nil {
		h.Metadata.AllTags, h.Metadata.Comments = parseTagsAndComments(h.Text)
		h.Chroms = getChromInfo(h.Metadata.AllTags)
		if _, ok := h.Metadata.AllTags[[2]byte{'H', 'D'}]; ok {
			h.Metadata.Version = getVersion(h.Metadata.AllTags)
			h.Metadata.SortOrder = getSortOrder(h.Metadata.AllTags)
			h.Metadata.Grouping = getGrouping(h.Metadata.AllTags)
		}
	}
	return h
}

// parseTagsAndComments parses header text into a HeaderTagMap and a slice of comment lines.
func parseTagsAndComments(text []string) (tags HeaderTagMap, comments []string) {
	var currTag Tag
	tags = make(HeaderTagMap)
	for _, line := range text {
		words := strings.Split(line, "\t")

		if words[0][:3] == "@CO" {
			comments = append(comments, line)
			continue
		}

		if words[0][0] != '@' || len(words[0]) != 3 {
			log.Fatalf("malformed sam header line: %s", line)
		}

		copy(currTag[:], words[0][1:]) // copy Tag

		tags[currTag] = append(tags[currTag], parseSubTags(words[1:]))
	}

	return
}

// parseSubTags parses a single line of tab delimited tagsets.
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

// getChromInfo further parses tags stored in a HeaderTagMap to extract ChromInfo.
func getChromInfo(tags HeaderTagMap) []chromInfo.ChromInfo {
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

// getVersion pulls the version number from a HeaderTagMap.
func getVersion(tags HeaderTagMap) string {
	return tags[[2]byte{'H', 'D'}][0][[2]byte{'V', 'N'}]
}

// getSortOrder pulls the sort order from a HeaderTagMap.
func getSortOrder(tags HeaderTagMap) []SortOrder {
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

// getGrouping pulls the grouping from a HeaderTagMap.
func getGrouping(tags HeaderTagMap) Grouping {
	return groupingMap[tags[[2]byte{'H', 'D'}][0][[2]byte{'G', 'O'}]]
}

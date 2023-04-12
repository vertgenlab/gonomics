// Package popgen contains tools for population genetic analysis, specifically for selection and population structure.
package popgen

import (
	"strings"

	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
)

// The Group struct and associated functions provides a system to partition a list of samples or species for subsequent analysis of subsections of the overall data.
type Group struct {
	Name    string
	Members []string
}

// ReadGroups parses a Group format file into a slice of Group structs.
func ReadGroups(filename string) []*Group {
	var line string
	var doneReading bool = false
	var index int64 = -1
	answer := make([]*Group, 0)
	//answer[0] = &Group{Name: "", Members: make([]string, 0)}
	//answer[1] = &Group{Name: "", Members: make([]string, 0)}

	groupFile := fileio.EasyOpen(filename)

	for line, doneReading = fileio.EasyNextRealLine(groupFile); !doneReading; line, doneReading = fileio.EasyNextRealLine(groupFile) {
		if strings.HasPrefix(line, ">") {
			tmp := Group{Name: line[1:], Members: make([]string, 0)}
			answer = append(answer, &tmp)
			index++
		} else {
			answer[index].Members = append(answer[index].Members, line)
		}
	}
	groupFile.Close()
	return answer
}

// GroupListsAreEqual returns true if all the groups in a list of groups are equal to another list of groups, false otherwise.
func GroupListsAreEqual(a []*Group, b []*Group) bool {
	if GroupListsCompare(a, b) != 0 {
		return false
	}
	return true
}

// GroupListsCompare compares two slices of Groups a and b for sorting and equality testing.
func GroupListsCompare(a []*Group, b []*Group) int {
	var res int
	stop := min(len(a), len(b))
	for i := 0; i < stop; i++ {
		res = GroupCompare(a[i], b[i])
		if res != 0 {
			return res
		}
	}
	if len(a) < len(b) {
		return -1
	} else if len(a) > len(b) {
		return 1
	}
	return 0
}

// GroupCompare compares two Groups a and b for sorting or equality testing.
func GroupCompare(a *Group, b *Group) int {
	if strings.Compare(a.Name, b.Name) != 0 {
		//DEBUG:log.Printf("It was the names. a: %s. b: %s.", a.Name, b.Name)
		return strings.Compare(a.Name, b.Name)
	}
	var res int
	stop := min(len(a.Members), len(b.Members))
	for i := 0; i < stop; i++ {
		res = strings.Compare(a.Members[i], b.Members[i])
		if res != 0 {
			//DEBUG:log.Printf("It was the members. a: %s. b: %s.", a.Members[i], b.Members[i])
			return res
		}
	}
	if len(a.Members) < len(b.Members) {
		return -1
	} else if len(a.Members) > len(b.Members) {
		return 1
	}
	return 0
}

// min is a local implmentation of a minimum int function to avoid a numbers import.
func min(a int, b int) int {
	if a < b {
		return a
	}
	return b
}

// GroupsContains returns true if any groups within a slice of groups g contains a string s, false otherwise.
func GroupsContains(g []*Group, s string) bool {
	for i := range g {
		if GroupContains(g[i], s) {
			return true
		}
	}
	return false
}

// GroupContains returns true if an input string s is contained within the members of group g, false otherwise.
func GroupContains(g *Group, s string) bool {
	for i := range g.Members {
		if strings.Compare(g.Members[i], s) == 0 {
			return true
		}
	}
	return false
}

// FindMissingGroupMembers returns a string of all of the entries in a Group slice that are not contained in the names of a multiFa alignment.
func FindMissingGroupMembers(aln []fasta.Fasta, g []*Group) string {
	var answer string = "Missing: "
	var missing bool = false
	for i := range g {
		answer = answer + g[i].Name + ": "
		for j := range g[i].Members {
			missing = true
			for k := 0; k < len(aln); k++ {
				if aln[k].Name == g[i].Members[j] {
					missing = false
				}
			}
			if missing {
				answer = answer + g[i].Members[j] + ", "
			}
		}
	}
	return answer
}

// FilterMultByGroup takes in a multiFa alignment returns a multiFa containing only the entries that are contained in an input slice of Group structs.
func FilterMultByGroup(aln []fasta.Fasta, g []*Group) []fasta.Fasta {
	var answer []fasta.Fasta
	for i := range aln {
		for j := range g {
			for k := 0; k < len(g[j].Members); k++ {
				if aln[i].Name == g[j].Members[k] {
					answer = append(answer, aln[i])
				}
			}
		}
	}
	return answer
}

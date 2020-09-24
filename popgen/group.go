//Package popgen contains tools for population genetic analysis, specifically for selection and population structure.
package popgen

import (
	"github.com/vertgenlab/gonomics/fasta"
	"github.com/vertgenlab/gonomics/fileio"
	"strings"
)

//The Group struct and associated functions provides a system to partition a list of samples or species for subsequent analysis of subsections of the overall data.
type Group struct {
	Name    string
	Members []string
}

//ReadGroups parses a Group format file into a slice of Group structs.
func ReadGroups(filename string) []*Group {
	var line string
	var doneReading bool = false
	var index int64 = -1
	answer := make([]*Group, 0)
	//answer[0] = &Group{Name: "", Members: make([]string, 0)}
	//answer[1] = &Group{Name: "", Members: make([]string, 0)}

	groupFile := fileio.EasyOpen(filename)
	defer groupFile.Close()

	for line, doneReading = fileio.EasyNextRealLine(groupFile); !doneReading; line, doneReading = fileio.EasyNextRealLine(groupFile) {
		if strings.HasPrefix(line, ">") {
			tmp := Group{Name: line[1:len(line)], Members: make([]string, 0)}
			answer = append(answer, &tmp)
			index++
		} else {
			answer[index].Members = append(answer[index].Members, line)
		}
	}
	return answer
}

//FindMissingGroupMembers returns a string of all of the entries in a Group slice that are not contained in the names of a multiFa alignment.
func FindMissingGroupMembers(aln []*fasta.Fasta, g []*Group) string {
	var answer string = "Missing: "
	var missing bool = false
	for i := 0; i < len(g); i++ {
		answer = answer + g[i].Name + ": "
		for j := 0; j < len(g[i].Members); j++ {
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

//FilterMultByGroup takes in a multiFa alignment returns a multiFa containing only the entries that are contained in an input slice of Group structs.
func FilterMultByGroup(aln []*fasta.Fasta, g []*Group) []*fasta.Fasta {
	var answer []*fasta.Fasta
	for i := 0; i < len(aln); i++ {
		for j := 0; j < len(g); j++ {
			for k := 0; k < len(g[j].Members); k++ {
				if aln[i].Name == g[j].Members[k] {
					answer = append(answer, aln[i])
				}
			}
		}
	}
	return answer
}

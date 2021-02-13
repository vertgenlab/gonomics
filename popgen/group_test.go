package popgen

import (
	"testing"
	"strings"
)

var GroupList []*Group = []*Group{&Group{Name: "apples", Members: strings.Split("Fuji:Honeycrisp:RedDelicious", ":")}, &Group{Name: "bananas", Members: strings.Split("Cavendish:Plantain", ":")}}

var GroupTests =[]struct {
	inputFile string
	expected []*Group
}{
	{"testdata/test.group", GroupList},
}

func TestGroupRead(t *testing.T) {
	for _, v := range GroupTests {
		records := ReadGroups(v.inputFile)
		if !GroupListsAreEqual(records, v.expected) {
			t.Errorf("Error in reading Group files.")
		}
	}
}

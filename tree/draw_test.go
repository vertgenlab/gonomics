package tree

import (
	"image/png"
	"os"
	"testing"
)

var treeDrawTests = []struct {
	treeText    string
	tmpFilename string
}{
	{"(human,chimp)ancestor;", "testdata/tmp.1.png"},
	{"((human,chimp)humanChimpAncestor,rhesus);", "testdata/tmp.2.png"},
	{"(((human,chimp),(mouse,rat)),dog);", "testdata/tmp.3.png"},
	{"((platypus,(opossum,((((rat,mouse),rabbit),human),dog))),(lizard,bird));", "testdata/tmp.4.png"},
	{"(human:0.5,chimp:0.2);", "testdata/tmp.5.png"},
	{"((human:0.5,chimp:0.2):0.3,rhesus:0.3);", "testdata/tmp.6.png"},
	{"(((human:0.5,chimp:0.2):0.3,(mouse:0.1,rat:0.6):0.2):0.1,dog:0.7);", "testdata/tmp.7.png"},
}

// TODO: This is not a great test because there is no real check, other
// than the code running without an error.  Maybe there is a better way?
func TestTreeDraw(t *testing.T) {
	for _, test := range treeDrawTests {
		tree, err := ParseNewick(test.treeText)
		if err != nil {
			t.Error(err)
		}
		img, err := Draw(tree, 1000, 200)
		if err != nil {
			t.Error(err)
		}
		imgOutFile, err := os.Create(test.tmpFilename)
		if err != nil {
			t.Error(err)
		}
		err = png.Encode(imgOutFile, img)
		if err != nil {
			t.Error(err)
		}
		imgOutFile.Close()
		os.Remove(test.tmpFilename)
	}
}

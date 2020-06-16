package interval

import (
	"github.com/vertgenlab/gonomics/bed"
	"reflect"
	"testing"
)

func TestQuery(t *testing.T) {
	var testIntervals []Interval = []Interval{
		&bed.Bed{Chrom: "1", ChromStart: 1, ChromEnd: 8},
		&bed.Bed{Chrom: "1", ChromStart: 3, ChromEnd: 5},
		&bed.Bed{Chrom: "1", ChromStart: 4, ChromEnd: 4},
		&bed.Bed{Chrom: "1", ChromStart: 5, ChromEnd: 9},
		&bed.Bed{Chrom: "1", ChromStart: 6, ChromEnd: 11},
		&bed.Bed{Chrom: "1", ChromStart: 7, ChromEnd: 7},
		&bed.Bed{Chrom: "1", ChromStart: 8, ChromEnd: 10},
		&bed.Bed{Chrom: "1", ChromStart: 9, ChromEnd: 10},
	}
	tree := BuildTree(testIntervals)

	q := &bed.Bed{Chrom: "1", ChromStart: 4, ChromEnd: 4}

	//TODO: build in more types of relationship tests
	answer := Query(tree, q, "e")

	if reflect.DeepEqual(answer[0].(*bed.Bed), *q) {
		t.Errorf("ERROR: Problem with querying tree")
	}
}

func TestBuildTree(t *testing.T) {
	var testIntervals []Interval = []Interval{
		&bed.Bed{Chrom: "1", ChromStart: 1, ChromEnd: 8},
		&bed.Bed{Chrom: "1", ChromStart: 3, ChromEnd: 5},
		&bed.Bed{Chrom: "1", ChromStart: 4, ChromEnd: 4},
		&bed.Bed{Chrom: "1", ChromStart: 5, ChromEnd: 9},
		&bed.Bed{Chrom: "1", ChromStart: 6, ChromEnd: 11},
		&bed.Bed{Chrom: "1", ChromStart: 7, ChromEnd: 7},
		&bed.Bed{Chrom: "1", ChromStart: 8, ChromEnd: 10},
		&bed.Bed{Chrom: "1", ChromStart: 9, ChromEnd: 10},
	}
	tree := BuildTree(testIntervals)
	if tree.xMid != 5 {
		t.Errorf("ERROR: Problem building tree")
	}
	if tree.lChild.xMid != 3 {
		t.Errorf("ERROR: Problem building tree")
	}
	if tree.rChild.xMid != 7 {
		t.Errorf("ERROR: Problem building tree")
	}
	if tree.lChild.lChild.xMid != 1 {
		t.Errorf("ERROR: Problem building tree")
	}
	if tree.lChild.rChild.xMid != 4 {
		t.Errorf("ERROR: Problem building tree")
	}
	if tree.rChild.lChild.xMid != 6 {
		t.Errorf("ERROR: Problem building tree")
	}
	if tree.rChild.rChild.xMid != 8 {
		t.Errorf("ERROR: Problem building tree")
	}
}

func TestBuildFCIndex(t *testing.T) {
	var testIntervals []Interval = []Interval{
		&bed.Bed{Chrom: "1", ChromStart: 1, ChromEnd: 8},
		&bed.Bed{Chrom: "1", ChromStart: 3, ChromEnd: 5},
		&bed.Bed{Chrom: "1", ChromStart: 4, ChromEnd: 4},
		&bed.Bed{Chrom: "1", ChromStart: 5, ChromEnd: 9},
		&bed.Bed{Chrom: "1", ChromStart: 6, ChromEnd: 11},
		&bed.Bed{Chrom: "1", ChromStart: 7, ChromEnd: 7},
		&bed.Bed{Chrom: "1", ChromStart: 8, ChromEnd: 10},
		&bed.Bed{Chrom: "1", ChromStart: 9, ChromEnd: 10},
	}

	pLeft := []Interval{
		&bed.Bed{Chrom: "1", ChromStart: 1, ChromEnd: 8},
		&bed.Bed{Chrom: "1", ChromStart: 3, ChromEnd: 5},
		&bed.Bed{Chrom: "1", ChromStart: 4, ChromEnd: 4},
		&bed.Bed{Chrom: "1", ChromStart: 5, ChromEnd: 9},
	}

	sortIntervals(testIntervals, yLess)
	sortIntervals(pLeft, yLess)

	answer := createFCIndex(testIntervals, pLeft)
	if answer[0] != 0 || answer[1] != 1 || answer[2] != 2 ||
		answer[3] != 2 || answer[4] != 3 || answer[5] != -1 ||
		answer[6] != -1 || answer[7] != -1 {
		t.Errorf("ERROR: Problem creating FC index")
	}
}

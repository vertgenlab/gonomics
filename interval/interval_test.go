package interval

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/common"
	"math/rand"
	"reflect"
	"testing"
	"time"
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

const (
	numIntervals = 1000
	rangeLow     = 0
	rangeHigh    = 100
)

func BenchmarkQueryO(b *testing.B) {
	rand.Seed(time.Now().UnixNano())
	intervals := generateIntervals(numIntervals, rangeLow, rangeHigh)
	tree := BuildTree(intervals)

	for n := 0; n < b.N; n++ {
		Query(tree, randInterval(rangeLow, rangeHigh), "o")
	}
}
func BenchmarkQueryOi(b *testing.B) {
	rand.Seed(time.Now().UnixNano())
	intervals := generateIntervals(numIntervals, rangeLow, rangeHigh)
	tree := BuildTree(intervals)

	for n := 0; n < b.N; n++ {
		Query(tree, randInterval(rangeLow, rangeHigh), "oi")
	}
}
func BenchmarkQueryD(b *testing.B) {
	rand.Seed(time.Now().UnixNano())
	intervals := generateIntervals(numIntervals, rangeLow, rangeHigh)
	tree := BuildTree(intervals)

	for n := 0; n < b.N; n++ {
		Query(tree, randInterval(rangeLow, rangeHigh), "d")
	}
}
func BenchmarkQueryDi(b *testing.B) {
	rand.Seed(time.Now().UnixNano())
	intervals := generateIntervals(numIntervals, rangeLow, rangeHigh)
	tree := BuildTree(intervals)

	for n := 0; n < b.N; n++ {
		Query(tree, randInterval(rangeLow, rangeHigh), "di")
	}
}
func BenchmarkQueryS(b *testing.B) {
	rand.Seed(time.Now().UnixNano())
	intervals := generateIntervals(numIntervals, rangeLow, rangeHigh)
	tree := BuildTree(intervals)

	for n := 0; n < b.N; n++ {
		Query(tree, randInterval(rangeLow, rangeHigh), "s")
	}
}
func BenchmarkQuerySi(b *testing.B) {
	rand.Seed(time.Now().UnixNano())
	intervals := generateIntervals(numIntervals, rangeLow, rangeHigh)
	tree := BuildTree(intervals)

	for n := 0; n < b.N; n++ {
		Query(tree, randInterval(rangeLow, rangeHigh), "si")
	}
}
func BenchmarkQueryF(b *testing.B) {
	rand.Seed(time.Now().UnixNano())
	intervals := generateIntervals(numIntervals, rangeLow, rangeHigh)
	tree := BuildTree(intervals)

	for n := 0; n < b.N; n++ {
		Query(tree, randInterval(rangeLow, rangeHigh), "f")
	}
}
func BenchmarkQueryFi(b *testing.B) {
	rand.Seed(time.Now().UnixNano())
	intervals := generateIntervals(numIntervals, rangeLow, rangeHigh)
	tree := BuildTree(intervals)

	for n := 0; n < b.N; n++ {
		Query(tree, randInterval(rangeLow, rangeHigh), "fi")
	}
}
func BenchmarkQueryM(b *testing.B) {
	rand.Seed(time.Now().UnixNano())
	intervals := generateIntervals(numIntervals, rangeLow, rangeHigh)
	tree := BuildTree(intervals)

	for n := 0; n < b.N; n++ {
		Query(tree, randInterval(rangeLow, rangeHigh), "m")
	}
}
func BenchmarkQueryMi(b *testing.B) {
	rand.Seed(time.Now().UnixNano())
	intervals := generateIntervals(numIntervals, rangeLow, rangeHigh)
	tree := BuildTree(intervals)

	for n := 0; n < b.N; n++ {
		Query(tree, randInterval(rangeLow, rangeHigh), "mi")
	}
}
func BenchmarkQueryE(b *testing.B) {
	rand.Seed(time.Now().UnixNano())
	intervals := generateIntervals(numIntervals, rangeLow, rangeHigh)
	tree := BuildTree(intervals)

	for n := 0; n < b.N; n++ {
		Query(tree, randInterval(rangeLow, rangeHigh), "e")
	}
}
func BenchmarkQueryLT(b *testing.B) {
	rand.Seed(time.Now().UnixNano())
	intervals := generateIntervals(numIntervals, rangeLow, rangeHigh)
	tree := BuildTree(intervals)

	for n := 0; n < b.N; n++ {
		Query(tree, randInterval(rangeLow, rangeHigh), "lt")
	}
}
func BenchmarkQueryGT(b *testing.B) {
	rand.Seed(time.Now().UnixNano())
	intervals := generateIntervals(numIntervals, rangeLow, rangeHigh)
	tree := BuildTree(intervals)

	for n := 0; n < b.N; n++ {
		Query(tree, randInterval(rangeLow, rangeHigh), "gt")
	}
}
func BenchmarkQueryAny(b *testing.B) {
	rand.Seed(time.Now().UnixNano())
	intervals := generateIntervals(numIntervals, rangeLow, rangeHigh)
	tree := BuildTree(intervals)

	for n := 0; n < b.N; n++ {
		Query(tree, randInterval(rangeLow, rangeHigh), "any")
	}
}
func BenchmarkQueryWithin(b *testing.B) {
	rand.Seed(time.Now().UnixNano())
	intervals := generateIntervals(numIntervals, rangeLow, rangeHigh)
	tree := BuildTree(intervals)

	for n := 0; n < b.N; n++ {
		Query(tree, randInterval(rangeLow, rangeHigh), "within")
	}
}
func BenchmarkQueryStart(b *testing.B) {
	rand.Seed(time.Now().UnixNano())
	intervals := generateIntervals(numIntervals, rangeLow, rangeHigh)
	tree := BuildTree(intervals)

	for n := 0; n < b.N; n++ {
		Query(tree, randInterval(rangeLow, rangeHigh), "start")
	}
}
func BenchmarkQueryEnd(b *testing.B) {
	rand.Seed(time.Now().UnixNano())
	intervals := generateIntervals(numIntervals, rangeLow, rangeHigh)
	tree := BuildTree(intervals)

	for n := 0; n < b.N; n++ {
		Query(tree, randInterval(rangeLow, rangeHigh), "end")
	}
}

type testInterval struct {
	chr   string
	start int
	end   int
}

func (t *testInterval) GetChrom() string {
	return t.chr
}
func (t *testInterval) GetChromStart() int {
	return t.start
}
func (t *testInterval) GetChromEnd() int {
	return t.end
}

func generateIntervals(num int, rangeLow int, rangeHigh int) []Interval {
	var answer []Interval

	for i := 0; i < num; i++ {
		answer = append(answer, randInterval(rangeLow, rangeHigh))
	}

	return answer
}

func randInterval(rangeLow int, rangeHigh int) Interval {
	var curr testInterval
	var currI Interval
	curr.chr = "1"
	curr.start = common.RandIntInRange(rangeLow, rangeHigh)
	curr.end = common.RandIntInRange(curr.start, rangeHigh)
	currI = &curr
	return currI
}

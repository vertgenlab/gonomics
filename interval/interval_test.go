package interval

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/numbers"
	"io"
	"math/rand"
	"testing"
	"time"
)

func TestQuery(t *testing.T) {
	type testQuery struct {
		query        *bed.Bed
		relationship string
		hasOverlap   bool
	}

	var targetIntervals []Interval = []Interval{
		&bed.Bed{Chrom: "1", ChromStart: 1, ChromEnd: 8},
		&bed.Bed{Chrom: "1", ChromStart: 3, ChromEnd: 5},
		&bed.Bed{Chrom: "1", ChromStart: 4, ChromEnd: 5},
		&bed.Bed{Chrom: "1", ChromStart: 5, ChromEnd: 9},
		&bed.Bed{Chrom: "1", ChromStart: 6, ChromEnd: 11},
		&bed.Bed{Chrom: "1", ChromStart: 7, ChromEnd: 8},
		&bed.Bed{Chrom: "1", ChromStart: 8, ChromEnd: 10},
		&bed.Bed{Chrom: "1", ChromStart: 9, ChromEnd: 10},
		&bed.Bed{Chrom: "chr2", ChromStart: 1, ChromEnd: 3},
		&bed.Bed{Chrom: "chr2", ChromStart: 4, ChromEnd: 10},
		&bed.Bed{Chrom: "chr3", ChromStart: 10, ChromEnd: 100},
	}

	var queries []testQuery = []testQuery{
		{
			query:        &bed.Bed{Chrom: "1", ChromStart: 4, ChromEnd: 5},
			relationship: "e",
			hasOverlap:   true,
		},
		{
			query:        &bed.Bed{Chrom: "chr2", ChromStart: 3, ChromEnd: 4},
			relationship: "any",
			hasOverlap:   false,
		},
		{
			query:        &bed.Bed{Chrom: "chr3", ChromStart: 98, ChromEnd: 101},
			relationship: "o",
			hasOverlap:   true,
		},
		{
			query:        &bed.Bed{Chrom: "chr3", ChromStart: 9, ChromEnd: 12},
			relationship: "oi",
			hasOverlap:   true,
		},
		{
			query:        &bed.Bed{Chrom: "chr3", ChromStart: 9, ChromEnd: 101},
			relationship: "d",
			hasOverlap:   true,
		},
		{
			query:        &bed.Bed{Chrom: "chr3", ChromStart: 11, ChromEnd: 99},
			relationship: "di",
			hasOverlap:   true,
		},
		{
			query:        &bed.Bed{Chrom: "chr3", ChromStart: 99, ChromEnd: 101},
			relationship: "m",
			hasOverlap:   true,
		},
		{
			query:        &bed.Bed{Chrom: "chr3", ChromStart: 9, ChromEnd: 11},
			relationship: "mi",
			hasOverlap:   true,
		},
		{
			query:        &bed.Bed{Chrom: "chr3", ChromStart: 10, ChromEnd: 101},
			relationship: "s",
			hasOverlap:   true,
		},
		{
			query:        &bed.Bed{Chrom: "chr3", ChromStart: 10, ChromEnd: 99},
			relationship: "si",
			hasOverlap:   true,
		},
		{
			query:        &bed.Bed{Chrom: "chr3", ChromStart: 9, ChromEnd: 100},
			relationship: "f",
			hasOverlap:   true,
		},
		{
			query:        &bed.Bed{Chrom: "chr3", ChromStart: 11, ChromEnd: 100},
			relationship: "fi",
			hasOverlap:   true,
		},
	}

	tree := BuildTree(targetIntervals)
	var answer []Interval

	//TODO: build in more types of relationship tests
	for _, q := range queries {
		answer = Query(tree, q.query, q.relationship)
		if q.hasOverlap && len(answer) == 0 {
			t.Errorf("Error: did not find interval during query of tree: %s %d %d", q.query.Chrom, q.query.ChromStart, q.query.ChromEnd)
		} else if !q.hasOverlap && len(answer) != 0 {
			t.Errorf("Error: found something when you should not have during query of tree.")
		}
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
	if tree["1"].xMid != 5 {
		t.Errorf("ERROR: Problem building tree")
	}
	if tree["1"].lChild.xMid != 3 {
		t.Errorf("ERROR: Problem building tree")
	}
	if tree["1"].rChild.xMid != 7 {
		t.Errorf("ERROR: Problem building tree")
	}
	if tree["1"].lChild.lChild.xMid != 1 {
		t.Errorf("ERROR: Problem building tree")
	}
	if tree["1"].lChild.rChild.xMid != 4 {
		t.Errorf("ERROR: Problem building tree")
	}
	if tree["1"].rChild.lChild.xMid != 6 {
		t.Errorf("ERROR: Problem building tree")
	}
	if tree["1"].rChild.rChild.xMid != 8 {
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
	numIntervals = 10000
	rangeLow     = 0
	rangeHigh    = 248956422
)

func BenchmarkBuildTree(b *testing.B) {
	rand.Seed(time.Now().UnixNano())
	intervals := generateIntervals(numIntervals, rangeLow, rangeHigh)
	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		BuildTree(intervals)
	}
}
func BenchmarkQueryO(b *testing.B) {
	rand.Seed(time.Now().UnixNano())
	intervals := generateIntervals(numIntervals, rangeLow, rangeHigh)
	tree := BuildTree(intervals)
	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		Query(tree, randInterval(rangeLow, rangeHigh), "o")
	}
}
func BenchmarkQueryOi(b *testing.B) {
	rand.Seed(time.Now().UnixNano())
	intervals := generateIntervals(numIntervals, rangeLow, rangeHigh)
	tree := BuildTree(intervals)
	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		Query(tree, randInterval(rangeLow, rangeHigh), "oi")
	}
}
func BenchmarkQueryD(b *testing.B) {
	rand.Seed(time.Now().UnixNano())
	intervals := generateIntervals(numIntervals, rangeLow, rangeHigh)
	tree := BuildTree(intervals)
	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		Query(tree, randInterval(rangeLow, rangeHigh), "d")
	}
}
func BenchmarkQueryDi(b *testing.B) {
	rand.Seed(time.Now().UnixNano())
	intervals := generateIntervals(numIntervals, rangeLow, rangeHigh)
	tree := BuildTree(intervals)
	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		Query(tree, randInterval(rangeLow, rangeHigh), "di")
	}
}
func BenchmarkQueryS(b *testing.B) {
	rand.Seed(time.Now().UnixNano())
	intervals := generateIntervals(numIntervals, rangeLow, rangeHigh)
	tree := BuildTree(intervals)
	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		Query(tree, randInterval(rangeLow, rangeHigh), "s")
	}
}
func BenchmarkQuerySi(b *testing.B) {
	rand.Seed(time.Now().UnixNano())
	intervals := generateIntervals(numIntervals, rangeLow, rangeHigh)
	tree := BuildTree(intervals)
	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		Query(tree, randInterval(rangeLow, rangeHigh), "si")
	}
}
func BenchmarkQueryF(b *testing.B) {
	rand.Seed(time.Now().UnixNano())
	intervals := generateIntervals(numIntervals, rangeLow, rangeHigh)
	tree := BuildTree(intervals)
	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		Query(tree, randInterval(rangeLow, rangeHigh), "f")
	}
}
func BenchmarkQueryFi(b *testing.B) {
	rand.Seed(time.Now().UnixNano())
	intervals := generateIntervals(numIntervals, rangeLow, rangeHigh)
	tree := BuildTree(intervals)
	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		Query(tree, randInterval(rangeLow, rangeHigh), "fi")
	}
}
func BenchmarkQueryM(b *testing.B) {
	rand.Seed(time.Now().UnixNano())
	intervals := generateIntervals(numIntervals, rangeLow, rangeHigh)
	tree := BuildTree(intervals)
	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		Query(tree, randInterval(rangeLow, rangeHigh), "m")
	}
}
func BenchmarkQueryMi(b *testing.B) {
	rand.Seed(time.Now().UnixNano())
	intervals := generateIntervals(numIntervals, rangeLow, rangeHigh)
	tree := BuildTree(intervals)
	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		Query(tree, randInterval(rangeLow, rangeHigh), "mi")
	}
}
func BenchmarkQueryE(b *testing.B) {
	rand.Seed(time.Now().UnixNano())
	intervals := generateIntervals(numIntervals, rangeLow, rangeHigh)
	tree := BuildTree(intervals)
	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		Query(tree, randInterval(rangeLow, rangeHigh), "e")
	}
}
func BenchmarkQueryLT(b *testing.B) {
	rand.Seed(time.Now().UnixNano())
	intervals := generateIntervals(numIntervals, rangeLow, rangeHigh)
	tree := BuildTree(intervals)
	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		Query(tree, randInterval(rangeLow, rangeHigh), "lt")
	}
}
func BenchmarkQueryGT(b *testing.B) {
	rand.Seed(time.Now().UnixNano())
	intervals := generateIntervals(numIntervals, rangeLow, rangeHigh)
	tree := BuildTree(intervals)
	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		Query(tree, randInterval(rangeLow, rangeHigh), "gt")
	}
}
func BenchmarkQueryAny(b *testing.B) {
	rand.Seed(time.Now().UnixNano())
	intervals := generateIntervals(numIntervals, rangeLow, rangeHigh)
	tree := BuildTree(intervals)
	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		Query(tree, randInterval(rangeLow, rangeHigh), "any")
	}
}
func BenchmarkQueryWithin(b *testing.B) {
	rand.Seed(time.Now().UnixNano())
	intervals := generateIntervals(numIntervals, rangeLow, rangeHigh)
	tree := BuildTree(intervals)
	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		Query(tree, randInterval(rangeLow, rangeHigh), "within")
	}
}
func BenchmarkQueryStart(b *testing.B) {
	rand.Seed(time.Now().UnixNano())
	intervals := generateIntervals(numIntervals, rangeLow, rangeHigh)
	tree := BuildTree(intervals)
	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		Query(tree, randInterval(rangeLow, rangeHigh), "start")
	}
}
func BenchmarkQueryEnd(b *testing.B) {
	rand.Seed(time.Now().UnixNano())
	intervals := generateIntervals(numIntervals, rangeLow, rangeHigh)
	tree := BuildTree(intervals)
	b.ResetTimer()
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
func (t *testInterval) WriteToFileHandle(file io.Writer) {
	return
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
	curr.start = numbers.RandIntInRange(rangeLow, rangeHigh)
	curr.end = numbers.RandIntInRange(curr.start, curr.start+1000)
	currI = &curr
	return currI
}

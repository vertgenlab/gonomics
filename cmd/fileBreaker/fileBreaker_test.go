package fileBreaker

import "testing"

var chunk1 = "bananas"
var chunk2 = "apples"
var chunk3 = "oranges"

func TestFileBreaker(t *testing.T) {
	FileBreaker("testdata/input.txt", "5", "testdata/fruits.txt")
}

package common

import (
	"log"
)

func ExitIfError(err error) {
	if err != nil {
		log.Fatal(err)
	}
}

func DigitsBaseTen(x int64) int {
	var count int = 1
	if x < 0 {
		x = -1 * x
		count++
	}
	for x >= 10 {
		x = x / 10
		count++
	}
	return count
}

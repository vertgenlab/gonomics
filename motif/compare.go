package motif

import (
	"fmt"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/numbers"
	"github.com/vertgenlab/gonomics/numbers/parse"
	"strings"
)

// AlmostEqualTest determines if floating-point numbers within two files are equal within a specified epsilon level.
func AlmostEqualTest(alpha, beta string, epsilon float64) bool {
	query, answer := fileio.Read(alpha), fileio.Read(beta)
	indexes := []int{7, 8} // Assuming these are the indexes you want to compare

	if len(query) != len(answer) {
		fmt.Errorf("Error: %s and %s do not have the same number of bed records...", alpha, beta)
		return false
	}

	for i := 0; i < len(query); i++ {
		queryFields := strings.Split(query[i], "\t")
		answerFields := strings.Split(answer[i], "\t")

		// Check if both lines have the same number of fields
		if len(queryFields) != len(answerFields) {
			fmt.Errorf("Error: Line %d of files %s and %s have different number of fields", i, alpha, beta)
			return false
		}

		// Compare the specific fields for near equality
		for _, index := range indexes {
			if index >= len(queryFields) || index >= len(answerFields) {
				fmt.Errorf("Error: Index out of range for line %d", i)
				return false
			}

			queryValue := parse.StringToFloat64(queryFields[index])
			answerValue := parse.StringToFloat64(answerFields[index])

			// Compare the parsed values for near equality
			if !numbers.AlmostEqual(queryValue, answerValue, epsilon) {
				fmt.Errorf("Error: Values on line %d at index %d are not almost equal: %v, %v", i, index, queryValue, answerValue)
				return false
			}
		}
	}
	return true
}

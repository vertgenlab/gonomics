package numbers

import (
	"fmt"
	"math"
	"testing"
)

var LogPowIntTests = []struct {
	x      float64
	y      int
	answer float64
}{
	{0.0, 0, 0.0},
	{9.0, 0, 0.0},
	{9.0, 10, 21.9722458},
	{14.0, 12608, 33273.2348},
	{0.05, 56401, -168962},
}

func TestLogPowInt(t *testing.T) {
	for _, test := range LogPowIntTests {
		calculated := LogPowInt(test.x, test.y)
		if math.Abs(calculated-test.answer)/test.answer > 0.000001 {
			fmt.Println(math.Abs(calculated-test.answer) / test.answer)
			t.Errorf("LogPowerInt for x: %f y:%d returned %f. Expected %f.", test.x, test.y, calculated, test.answer)
		}
	}
}

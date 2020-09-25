package numbers

import (
	//DEBUG: "fmt"
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
			t.Errorf("LogPowerInt for x: %f y:%d returned %f. Expected %f.", test.x, test.y, calculated, test.answer)
		}
	}
}

var MultiplyLogTests = []struct {
	x      float64
	y      float64
	answer float64
}{
	{math.Inf(-1), 3, math.Inf(-1)},
	{4, math.Inf(1), math.Inf(1)},
	{4, 12, 16},
}

func TestMultiplyLog(t *testing.T) {
	for _, test := range MultiplyLogTests {
		calculated := MultiplyLog(test.x, test.y)
		if math.Abs(calculated-test.answer)/test.answer > 0.000001 {
			t.Errorf("TestMultiplyLog for x: %f y:%f returned %f. Expected %f.", test.x, test.y, calculated, test.answer)
		}
	}
}

var DivideLogTests = []struct {
	x      float64
	y      float64
	answer float64
}{
	{math.Inf(-1), 400000.0, math.Inf(-1)},
	{16, 4, 12},
}

func TestDivideLog(t *testing.T) {
	for _, test := range DivideLogTests {
		calculated := DivideLog(test.x, test.y)
		if math.Abs(calculated-test.answer)/test.answer > 0.000001 {
			t.Errorf("TestDivideLog for x: %f y:%f returned %f. Expected %f.", test.x, test.y, calculated, test.answer)
		}
	}
}

var AddLogTests = []struct {
	x      float64
	y      float64
	answer float64
}{
	{math.Inf(-1), 4100, 4100},
	{4100, math.Inf(-1), 4100},
	{4100, 3, 4100},
	{2, 3, 3.313262},
	{20, 15, 20.0067153}, //from wolframAlpha
}

func TestAddLog(t *testing.T) {
	for _, test := range AddLogTests {
		calculated := AddLog(test.x, test.y)
		if math.Abs(calculated-test.answer)/test.answer > 0.000001 {
			t.Errorf("TestAddLog for x: %f y:%f returned %f. Expected %f.", test.x, test.y, calculated, test.answer)
		}
	}
}

var SubtractLogTests = []struct {
	x      float64
	y      float64
	answer float64
}{
	{20, 20, math.Inf(-1)},
	{400, math.Inf(-1), 400},
	{410000, 8, 410000},
	{420, 418, 419.8545865},
}

func TestSubtractLog(t *testing.T) {
	for _, test := range SubtractLogTests {
		calculated := SubtractLog(test.x, test.y)
		if math.Abs(calculated-test.answer)/test.answer > 0.000001 {
			t.Errorf("TestSubtractLog for x: %f y:%f returned %f. Expected %f.", test.x, test.y, calculated, test.answer)
		}
	}
}

var MidpointLogTests = []struct {
	x      float64
	y      float64
	answer float64
}{
	{16, 13, 15.35544},
	{2, 3, 2.620115},
}

func TestMidpointLog(t *testing.T) {
	for _, test := range MidpointLogTests {
		calculated := MidpointLog(test.x, test.y)
		if math.Abs(calculated-test.answer)/test.answer > 0.000001 {
			t.Errorf("TestMidpointLog for x: %f y:%f returned %f. Expected %f.", test.x, test.y, calculated, test.answer)
		}
	}
}

var AverageLogTests = []struct {
	x      []float64
	answer float64
}{
	{[]float64{10.0, 10.0, 10.0, 10.0}, 10.0},
	{[]float64{10.0, 8.0, 14.0, 4.0}, 12.6343},
}

//Note: AverageLog is significantly less accurate than other functions, due to the inaccuracy of AddLog, which compounds here.
//TODO: Maybe there's a smarter way to do this for more accuracy. Obviously averages in normal space are better, this should be used for rough estimates if the numbers can only be expressed in logSpace.
func TestAverageLog(t *testing.T) {
	for _, test := range AverageLogTests {
		calculated := AverageLog(test.x)
		if math.Abs(calculated-test.answer)/test.answer > 0.01 {
			t.Errorf("TestAverageLog for x: %f returned %f. Expected %f.", test.x, calculated, test.answer)
		}
	}
}

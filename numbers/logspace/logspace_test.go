package logspace

import (
	"math"
	"testing"
)

var PowTests = []struct {
	x      float64
	y      float64
	answer float64
}{
	{math.Log(0.0), 0, 0.0},
	{math.Log(9.0), 0, 0.0},
	{math.Log(9.0), 10, 21.9722458},
	{math.Log(14.0), 12608, 33273.2348},
	{math.Log(0.05), 56401, -168962},
}

func TestPow(t *testing.T) {
	for _, test := range PowTests {
		calculated := Pow(test.x, test.y)
		if math.Abs(calculated-test.answer)/test.answer > 0.000001 {
			t.Errorf("LogPow for x: %f y: %f returned %f. Expected %f.", test.x, test.y, calculated, test.answer)
		}
	}
}

var MultiplyTests = []struct {
	x      float64
	y      float64
	answer float64
}{
	{math.Inf(-1), 3, math.Inf(-1)},
	{4, math.Inf(1), math.Inf(1)},
	{4, 12, 16},
}

func TestMultiply(t *testing.T) {
	for _, test := range MultiplyTests {
		calculated := Multiply(test.x, test.y)
		if math.Abs(calculated-test.answer)/test.answer > 0.000001 {
			t.Errorf("TestMultiply for x: %f y:%f returned %f. Expected %f.", test.x, test.y, calculated, test.answer)
		}
	}
}

var DivideTests = []struct {
	x      float64
	y      float64
	answer float64
}{
	{math.Inf(-1), 400000.0, math.Inf(-1)},
	{16, 4, 12},
}

func TestDivide(t *testing.T) {
	for _, test := range DivideTests {
		calculated := Divide(test.x, test.y)
		if math.Abs(calculated-test.answer)/test.answer > 0.000001 {
			t.Errorf("TestDivide for x: %f y:%f returned %f. Expected %f.", test.x, test.y, calculated, test.answer)
		}
	}
}

var AddTests = []struct {
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

func TestAdd(t *testing.T) {
	for _, test := range AddTests {
		calculated := Add(test.x, test.y)
		if math.Abs(calculated-test.answer)/test.answer > 0.000001 {
			t.Errorf("TestAdd for x: %f y:%f returned %f. Expected %f.", test.x, test.y, calculated, test.answer)
		}
	}
}

var SubtractTests = []struct {
	x      float64
	y      float64
	answer float64
}{
	{20, 20, math.Inf(-1)},
	{400, math.Inf(-1), 400},
	{410000, 8, 410000},
	{420, 418, 419.8545865},
}

func TestSubtract(t *testing.T) {
	for _, test := range SubtractTests {
		calculated := Subtract(test.x, test.y)
		if math.Abs(calculated-test.answer)/test.answer > 0.000001 {
			t.Errorf("TestSubtractLog for x: %f y:%f returned %f. Expected %f.", test.x, test.y, calculated, test.answer)
		}
	}
}

var AverageTests = []struct {
	x      float64
	y      float64
	answer float64
}{
	{16, 13, 15.35544},
	{2, 3, 2.620115},
}

func TestAverage(t *testing.T) {
	for _, test := range AverageTests {
		calculated := Average(test.x, test.y)
		if math.Abs(calculated-test.answer)/test.answer > 0.000001 {
			t.Errorf("TestAverage for x: %f y:%f returned %f. Expected %f.", test.x, test.y, calculated, test.answer)
		}
	}
}

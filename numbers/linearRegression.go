package numbers

import (
	"fmt"
	"log"
	"math"
)

type matrix [][]float64

func LinearRegressionOld(data [][]float64, learningRate float64, maxInteration int, tol float64) {
	if learningRate <= 0 || learningRate >= 1 {
		log.Fatal()
	}
	var m, c, prevErr, err float64
	var pred [][]float64
	m = 3
	c = 0
	for i := 0; i < maxInteration; i++ {
		err, pred = MSE(m, c, data)
		if math.Abs(prevErr-err) < 1e-11 {
			fmt.Println("iteration: ", i)
			break
		}
		slopeGradient := o2(pred, data)
		interceptGradent := o1(pred, data)
		m = m - learningRate*slopeGradient
		c = c - learningRate*interceptGradent
		fmt.Println("cost: ", err)
		fmt.Printf("y = %f x + %f\n", m, c)
		prevErr = err
	}
}

func GradientDescent(start, tol, learnRate float64, gradient func(float64) float64) []float64 {
	var steps []float64
	var diff float64
	val := start
	steps = append(steps, val)
	for i := 0; i < 1000; i++ {
		diff = learnRate * gradient(val)
		if math.Abs(diff) < tol {
			break
		}
		val -= diff
		steps = append(steps, val)
	}
	return steps
}

func MSE(slope float64, intercept float64, data [][]float64) (float64, [][]float64) {
	var p, err, sum float64
	var pred [][]float64
	for i := range data {
		p = slope*data[i][0] + intercept
		pred = append(pred, []float64{data[i][0], p})
		err = p - data[i][1]
		sum += math.Pow(err, 2)
	}
	return sum / float64(len(data)), pred
}

func o1(pred, data [][]float64) float64 {
	var sum float64
	for i := range pred {
		sum += (data[i][1] - pred[i][1]) * -1
	}
	return sum / float64(len(data))
}

func o2(pred, data [][]float64) float64 {
	var sum float64
	for i := range pred {
		sum += ((data[i][1] - pred[i][1]) * -1) * data[i][0]
	}
	return 2 * sum / float64(len(data))
}

func SimulateLinearData(n int, slope, intercept, maxErr float64) [][]float64 {
	var data [][]float64
	var ideal, rand float64
	for i := 0; i < n; i++ {
		ideal = float64(i)*slope + intercept
		rand = RandFloat64InRange(-maxErr, maxErr)
		data = append(data, []float64{float64(i), ideal + rand})
	}
	return data
}

package numbers

import (
	"fmt"
	"math"
	"math/rand"
	"testing"
	"time"
)

func TestRunMixtureModel(t *testing.T) {
	rand.Seed(time.Now().UnixNano())
	data := []float64{15, 23, 23, 25, 25, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 33, 35, 35, 35, 35, 37, 37, 37, 37, 37, 37, 37, 37, 37, 40}
	//for i := 0; i < 5; i++ {
	//	data = append(data, data...)
	//}

	//numComponents := 2
	//countPerComponent := 100
	//means := []float64{120, 480, 240, 360, 320}
	//std := []float64{20, 30, 10, 40, 5}
	//data := make([]float64, 0, 200)
	//for i := 0; i < numComponents; i++ {
	//	data = append(data, generateData(countPerComponent, means[i], std[i])...)
	//}

	maxIterations := 200
	mm := new(MixtureModel)
	for i := 0; i < 100; i++ {
		RunMixtureModel(data, 2, maxIterations, mm)
		fmt.Println(mm.Means, math.Abs(mm.Means[0]-mm.Means[1]))
		//fmt.Println(mm.Variances)
		//fmt.Println()
	}
}

func generateData(num int, mean, std float64) []float64 {
	ans := make([]float64, num)
	var noise float64
	for i := range ans {
		noise = 2 * std * rand.Float64()
		noise -= 10
		ans[i] = ((rand.NormFloat64() * std) + mean) + noise
	}
	return ans
}

package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

//Tests were made with the following cmds.
//~/go/bin/multiFaAcceleration -windowSize 50 -rawInitialVelBranchLengths test.RawInitial.bed -rawVelBranchLengths test.RawVel.bed chr1 test.fa test.vel.expected.bed test.accel.expected.bed test.initialVel.expected.bed
//~/go/bin/multiFaAcceleration -windowSize 50 -searchSpaceBed test.searchspace.bed chr1 test.fa test.vel.searchspace.expected.bed test.accel.searchspace.expected.bed test.initialVel.searchspace.expected.bed
//~/go/bin/multiFaAcceleration -windowSize 50 -searchSpaceBed test.searchspace.bed -useSnpDistance chr1 test.fa test.vel.snpDistance.expected.bed test.accel.snpDistance.expected.bed test.initialVel.snpDistance.expected.bed

var multiFaAccelerationTests = []struct {
	InFile                     string
	ChromName                  string
	VelOut                     string
	AccelOut                   string
	InitialVelOut              string
	SearchSpaceBed             string
	SearchSpaceProportion      float64
	WindowSize                 int
	Verbose                    bool
	VelExpected                string
	AccelExpected              string
	InitialVelExpected         string
	UseSnpDistance             bool
	Epsilon                    float64
	AllowNegative              bool
	ZeroDistanceWeightConstant float64
	RawVelOut                  string
	RawInitialOut              string
	RawVelExpected             string
	RawInitialExpected         string
	CavalliSforzaQ             bool
}{
	{"testdata/test.fa",
		"chr1",
		"testdata/test.vel.bed",
		"testdata/test.accel.bed",
		"testdata/test.initialVel.bed",
		"",
		0.5,
		50,
		false,
		"testdata/test.vel.expected.bed",
		"testdata/test.accel.expected.bed",
		"testdata/test.initialVel.expected.bed",
		false,
		1e-8,
		false,
		1000,
		"testdata/test.RawVel.bed",
		"testdata/test.RawInitial.bed",
		"testdata/expected.RawVel.bed",
		"testdata/expected.RawInitial.bed",
		false,
	},
	{"testdata/test.fa",
		"chr1",
		"testdata/test.vel.bed",
		"testdata/test.accel.bed",
		"testdata/test.initialVel.bed",
		"testdata/test.searchspace.bed",
		0.5,
		50,
		false,
		"testdata/test.vel.searchspace.expected.bed",
		"testdata/test.accel.searchspace.expected.bed",
		"testdata/test.initialVel.searchspace.expected.bed",
		false,
		1e-8,
		false,
		1000,
		"",
		"",
		"",
		"",
		false,
	},
	{"testdata/test.fa",
		"chr1",
		"testdata/test.vel.bed",
		"testdata/test.accel.bed",
		"testdata/test.initialVel.bed",
		"testdata/test.searchspace.bed",
		0.5,
		50,
		false,
		"testdata/test.vel.snpDistance.expected.bed",
		"testdata/test.accel.snpDistance.expected.bed",
		"testdata/test.initialVel.snpDistance.expected.bed",
		true,
		1e-8,
		false,
		1000,
		"",
		"",
		"",
		"",
		false,
	},
}

func TestMultiFaAcceleration(t *testing.T) {
	var err error
	for _, v := range multiFaAccelerationTests {
		s := Settings{
			InFile:                     v.InFile,
			ChromName:                  v.ChromName,
			VelOut:                     v.VelOut,
			AccelOut:                   v.AccelOut,
			InitialVelOut:              v.InitialVelOut,
			SearchSpaceBed:             v.SearchSpaceBed,
			SearchSpaceProportion:      v.SearchSpaceProportion,
			WindowSize:                 v.WindowSize,
			Verbose:                    v.Verbose,
			UseSnpDistance:             v.UseSnpDistance,
			Epsilon:                    v.Epsilon,
			AllowNegative:              v.AllowNegative,
			ZeroDistanceWeightConstant: v.ZeroDistanceWeightConstant,
			RawVelOut:                  v.RawVelOut,
			RawInitialOut:              v.RawInitialOut,
			CavalliSforzaQ:             v.CavalliSforzaQ,
		}
		multiFaAcceleration(s)
		if !fileio.AreEqual(v.VelOut, v.VelExpected) {
			t.Errorf("Error in multiFaAcceleration, velOut did not match expected.")
		} else {
			err = os.Remove(v.VelOut)
			exception.PanicOnErr(err)
		}
		if !fileio.AreEqual(v.AccelOut, v.AccelExpected) {
			t.Errorf("Error in multiFaAcceleration, accelOut did not match expected.")
		} else {
			err = os.Remove(v.AccelOut)
			exception.PanicOnErr(err)
		}
		if !fileio.AreEqual(v.InitialVelOut, v.InitialVelExpected) {
			t.Errorf("Error in multiFaAcceleration, initialVelOut did not match expected.")
		} else {
			err = os.Remove(v.InitialVelOut)
			exception.PanicOnErr(err)
		}
		if v.RawVelOut != "" {
			if !fileio.AreEqual(v.RawVelOut, v.RawVelExpected) {
				t.Errorf("Error in multiFaAcceleration, RawVel did not match expected.")
			} else {
				err = os.Remove(v.RawVelOut)
				exception.PanicOnErr(err)
			}
		}
		if v.RawInitialOut != "" {
			if !fileio.AreEqual(v.RawInitialOut, v.RawInitialExpected) {
				t.Errorf("Error in multiFaAcceleration, RawInitial did not match expected.")
			} else {
				err = os.Remove(v.RawInitialOut)
				exception.PanicOnErr(err)
			}
		}
	}
}

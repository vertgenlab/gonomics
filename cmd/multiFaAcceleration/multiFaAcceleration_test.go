package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"os"
	"testing"
)

var multiFaAccelerationTests = []struct {
	InFile                string
	ChromName             string
	VelOut                string
	AccelOut              string
	InitialVelOut         string
	SearchSpaceBed        string
	SearchSpaceProportion float64
	WindowSize            int
	Verbose               bool
	VelExpected           string
	AccelExpected         string
	InitialVelExpected    string
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
	},
}

func TestMultiFaAcceleration(t *testing.T) {
	var err error
	for _, v := range multiFaAccelerationTests {
		s := Settings{
			InFile:                v.InFile,
			ChromName:             v.ChromName,
			VelOut:                v.VelOut,
			AccelOut:              v.AccelOut,
			InitialVelOut:         v.InitialVelOut,
			SearchSpaceBed:        v.SearchSpaceBed,
			SearchSpaceProportion: v.SearchSpaceProportion,
			WindowSize:            v.WindowSize,
			Verbose:               v.Verbose,
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
	}
}

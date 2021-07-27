package main

import (
	"testing"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/exception"
	"os"
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
	VelExpected string
	AccelExpected string
	InitialVelExpected string
}{
}

func TestMultiFaAcceleration(t *testing.T) {
	var err error
	for _, v := range multiFaAccelerationTests {
		s := Settings{
			InFile: v.InFile,
			ChromName: v.ChromName,
			VelOut: v.VelOut,
			AccelOut: v.AccelOut,
			InitialVelOut: v.InitialVelOut,
			SearchSpaceBed: v.SearchSpaceBed,
			SearchSpaceProportion: v.SearchSpaceProportion,
			WindowSize: v.WindowSize,
			Verbose: v.Verbose,
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

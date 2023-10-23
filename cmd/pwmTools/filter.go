package main

import (
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/motif"
)

type FilterSettings struct {
	InFile     string
	OutFile    string
	MatrixType string
	MinLength  int
	MaxLength  int
}

func pwmFilter(s FilterSettings) {
	var err error
	var pass bool
	records := motif.ReadJaspar(s.InFile, s.MatrixType)
	out := fileio.EasyCreate(s.OutFile)

	for i := range records {
		pass = true
		if len(records[i].Mat[0]) < s.MinLength {
			pass = false
		} else if len(records[i].Mat[0]) > s.MaxLength {
			pass = false
		}
		if pass {
			motif.WritePositionMatrixJaspar(out, records[i])
		}
	}

	err = out.Close()
	exception.PanicOnErr(err)
}

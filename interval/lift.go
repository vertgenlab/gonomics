package interval

import (
	"github.com/vertgenlab/gonomics/fileio"
)

type Lift interface {
	GetChrom() string
	GetChromStart() int
	GetChromEnd() int
	WriteToFileHandle(*fileio.EasyWriter)
	UpdateLift(string, int, int)
}

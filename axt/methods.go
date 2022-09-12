package axt

import (
	"io"
)

func (a Axt) GetChrom() string {
	return a.RName
}

func (a Axt) GetChromStart() int {
	return int(a.RStart - 1)
}

func (a Axt) GetChromEnd() int {
	return int(a.REnd)
}

func (a *Axt) UpdateLift(c string, start int, end int) {
	a.RName = c
	a.RStart = start
	a.REnd = end
}

func (a Axt) WriteToFileHandle(file io.Writer) {
	WriteToFileHandle(file, a, 0)
}

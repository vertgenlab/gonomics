package gtf

import (
	"io"
)

func (g *Gene) GetChrom() string {
	return g.Transcripts[0].Chr
}

func (g *Gene) GetChromStart() int {
	return g.Transcripts[0].Start - 1
}

func (g *Gene) GetChromEnd() int {
	return g.Transcripts[0].End
}

func (g *Gene) WriteToFileHandle(file io.Writer) {
	WriteToFileHandle(file, g)
}

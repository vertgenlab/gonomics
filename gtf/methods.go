package gtf

import (
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fileio"
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

func (g *Gene) WriteToFileHandle(file *fileio.EasyWriter) {
	err := WriteToFileHandle(file, g)
	common.ExitIfError(err)
}

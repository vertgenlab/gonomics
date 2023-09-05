package gtf

import (
	"io"
)

// GetChrom returns the name of the chromosome where the gene is located.
func (g *Gene) GetChrom() string {
	return g.Transcripts[0].Chr
}

// GetChromStart returns the genomic coordinate where the gene starts.
func (g *Gene) GetChromStart() int {
	return g.Transcripts[0].Start - 1
}

// GetChromEnd returns the genomics coordinate where the gene ends.
func (g *Gene) GetChromEnd() int {
	return g.Transcripts[0].End
}

// WriteToFileHandle writes the gene to an io.Writer
func (g *Gene) WriteToFileHandle(file io.Writer) {
	WriteToFileHandle(file, g)
}

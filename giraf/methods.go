package giraf

import (
	"fmt"
	"io"
)

func (g *Giraf) GetChrom() string {
	return fmt.Sprintf("%d", g.Path.Nodes[0])
}

func (g *Giraf) GetChromStart() int {
	return g.Path.TStart
}

func (g *Giraf) GetChromEnd() int {
	// BEWARE this returns the end position on the final node
	// There may be a gotcha if you assume GetChromStart() and End()
	// refer to positions on the same chromosome
	// Not sure how best to handle this???
	return g.Path.TEnd
}

func (g *Giraf) WriteToFileHandle(file io.Writer) {
	WriteGirafToFileHandle(file, g)
}

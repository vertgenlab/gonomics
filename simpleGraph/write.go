package simpleGraph

import (
	"bytes"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/giraf"
	"io"
	"sync"
)

func SimpleWriteGirafPair(filename string, input <-chan GirafGsw, wg *sync.WaitGroup) {
	file := fileio.EasyCreate(filename)
	var buf *bytes.Buffer

	var simplePool = sync.Pool{
		New: func() interface{} {
			return &bytes.Buffer{}
		},
	}
	var err error
	for gp := range input {
		buf = simplePool.Get().(*bytes.Buffer)
		_, err = buf.WriteString(giraf.ToString(&gp.ReadOne))
		common.ExitIfError(err)
		err = buf.WriteByte('\n')
		common.ExitIfError(err)
		_, err = buf.WriteString(giraf.ToString(&gp.ReadTwo))
		common.ExitIfError(err)
		err = buf.WriteByte('\n')
		common.ExitIfError(err)

		io.Copy(file, buf)
		buf.Reset()
		simplePool.Put(buf)
	}
	file.Close()
	wg.Done()
}

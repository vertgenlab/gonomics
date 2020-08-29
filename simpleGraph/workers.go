package simpleGraph

import (
	"bytes"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/giraf"
	"io"
	"sync"
)

/*
// this is just for speed testing to see how much of the speed slowdown is due to alignment time
func gswNothingWorker(gg *SimpleGraph, seedHash map[uint64][]uint64, seedLen int, stepSize int, incomingFastqs <-chan *fastq.FastqBig, outgoingSams chan<- *sam.SamAln, wg *sync.WaitGroup) {
	var curr *sam.SamAln
	for read := range incomingFastqs {
		curr = &sam.SamAln{QName: read.Name, Flag: 4, RName: "*", Pos: 0, MapQ: 255, Cigar: []*cigar.Cigar{&cigar.Cigar{Op: '*'}}, RNext: "*", PNext: 0, TLen: 0, Seq: read.Seq, Qual: string(read.Qual), Extra: "BZ:i:0"}
		outgoingSams <- curr
	}
	wg.Done()
}

func gswWorkerMemPool(gg *SimpleGraph, seedHash map[uint64][]uint64, seedLen int, stepSize int, scoreMatrix [][]int64, incomingFastqs <-chan *fastq.FastqBig, outgoingSams chan<- *sam.SamAln, wg *sync.WaitGroup) {
	m, trace := swMatrixSetup(10000)
	memChunk := make([]SeedDev, 100000)
	for i := 0; i < len(memChunk)-1; i++ {
		memChunk[i].Next = &memChunk[i+1]
	}
	memStart := &(memChunk[0])
	for read := range incomingFastqs {
		outgoingSams <- GraphSmithWatermanMemPool(gg, read, seedHash, seedLen, stepSize, scoreMatrix, m, trace, &memStart)
	}
	wg.Done()
}*/

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

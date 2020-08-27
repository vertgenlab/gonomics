package simpleGraph

import (
	"bytes"

	"github.com/vertgenlab/gonomics/cigar"
	"github.com/vertgenlab/gonomics/common"
	"github.com/vertgenlab/gonomics/dna"
	"github.com/vertgenlab/gonomics/fastq"
	"github.com/vertgenlab/gonomics/fileio"
	"github.com/vertgenlab/gonomics/giraf"
	"io"
	"strconv"
	"sync"
)

/*
func WriteSimpleGiraf(filename string, input <-chan GirafGsw, wg *sync.WaitGroup) {
	file := fileio.EasyCreate(filename)
	var buf *bytes.Buffer

	var simplePool = sync.Pool{
		New: func() interface{} {
			return &bytes.Buffer{}
		},
	}

	for gp := range input {
		buf = simplePool.Get().(*bytes.Buffer)
		buf = GirafStringBuilder(gp, buf)
		io.Copy(file, buf)
		buf.Reset()
		simplePool.Put(buf)
	}
	file.Close()
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
	for gp := range input {
		buf = simplePool.Get().(*bytes.Buffer)
		buf = GirafStringBuilder(gp.ReadOne, buf)
		buf = GirafStringBuilder(gp.ReadTwo, buf)
		io.Copy(file, buf)
		buf.Reset()
		simplePool.Put(buf)
	}
	file.Close()
	wg.Done()
}

func GirafStringBuilder(g giraf.Giraf, buf *bytes.Buffer) *bytes.Buffer {
	var err error
	buf = buildGirafString(buf, g, err)
	err = buf.WriteByte('\n')
	common.ExitIfError(err)
	return buf
}

func buildGirafString(buf *bytes.Buffer, g giraf.Giraf, err error) *bytes.Buffer {
	_, err = buf.WriteString(g.QName)
	common.ExitIfError(err)
	err = buf.WriteByte('\t')
	common.ExitIfError(err)
	_, err = buf.WriteString(strconv.Itoa(g.QStart))
	common.ExitIfError(err)
	err = buf.WriteByte('\t')
	common.ExitIfError(err)
	_, err = buf.WriteString(strconv.Itoa(g.QEnd))
	common.ExitIfError(err)
	err = buf.WriteByte('\t')
	common.ExitIfError(err)
	_, err = buf.WriteString(strconv.Itoa(int(g.Flag)))
	common.ExitIfError(err)
	err = buf.WriteByte('\t')
	common.ExitIfError(err)
	_, err = buf.WriteRune(common.StrandToRune(g.PosStrand))
	common.ExitIfError(err)
	err = buf.WriteByte('\t')
	common.ExitIfError(err)
	_, err = buf.WriteString(giraf.PathToString(g.Path))
	common.ExitIfError(err)
	err = buf.WriteByte('\t')
	common.ExitIfError(err)
	_, err = buf.WriteString(cigar.ByteCigarToString(g.ByteCigar))
	common.ExitIfError(err)
	err = buf.WriteByte('\t')
	common.ExitIfError(err)
	_, err = buf.WriteString(strconv.Itoa(g.AlnScore))
	common.ExitIfError(err)
	err = buf.WriteByte('\t')
	common.ExitIfError(err)
	_, err = buf.WriteString(strconv.Itoa(int(g.MapQ)))
	common.ExitIfError(err)
	err = buf.WriteByte('\t')
	common.ExitIfError(err)
	_, err = buf.WriteString(dna.ByteDnaBasesToString(g.Seq))
	common.ExitIfError(err)
	err = buf.WriteByte('\t')
	common.ExitIfError(err)
	_, err = buf.WriteString(fastq.Uint8QualToString(g.Qual))
	//common.ExitIfError(err)
	//err = buf.WriteByte('\t')
	common.ExitIfError(err)
	_, err = buf.WriteString(giraf.NotesToString(g.Notes))
	common.ExitIfError(err)
	return buf
}

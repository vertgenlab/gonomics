package bigWig

import (
	"encoding/binary"
	"github.com/vertgenlab/gonomics/exception"
	"github.com/vertgenlab/gonomics/fileio"
	"log"
)

// Read parses a BigWig struct from an input filename.
func Read(filename string) BigWig {
	var err error
	var answer = BigWig{}
	file := fileio.EasyOpen(filename)
	answer.BbiHeader = readBbiHeader(file)
	for currZoomLevel := 0; currZoomLevel < int(answer.BbiHeader.ZoomLevels); currZoomLevel++ {
		answer.ZoomHeaders = append(answer.ZoomHeaders, readZoomHeader(file))
	}
	answer.TotalSummaryBlock = readTotalSummary(file)
	answer.ChromTreeHeader = readChromTreeHeader(file)
	answer.ChromTreeNodes = readChromTreeNodes(file, answer.ChromTreeHeader)
	err = binary.Read(file, binary.LittleEndian, &answer.DataCount)
	exception.PanicOnErr(err)

	return answer
}

// We currently support only reading bbi headers in LittleEndian.
func readBbiHeader(file *fileio.EasyReader) BbiHeader {
	var header = BbiHeader{}
	err := binary.Read(file, binary.LittleEndian, &header.Magic)
	exception.PanicOnErr(err)
	if header.Magic == bigWigMagicBigEndian {
		log.Fatalf("Error: bigWig file appears to be in big endian. Current functionality only supports little endian bigWig files.\n")
	}
	if header.Magic != bigWigMagic {
		log.Fatalf("Error: bigWig magic was not as expected. Found: %v. Expected: %v.\n", header.Magic, bigWigMagic)
	}
	err = binary.Read(file, binary.LittleEndian, &header.Version)
	exception.PanicOnErr(err)
	err = binary.Read(file, binary.LittleEndian, &header.ZoomLevels)
	exception.PanicOnErr(err)
	err = binary.Read(file, binary.LittleEndian, &header.ChromosomeTreeOffset)
	exception.PanicOnErr(err)
	err = binary.Read(file, binary.LittleEndian, &header.FullDataOffset)
	exception.PanicOnErr(err)
	err = binary.Read(file, binary.LittleEndian, &header.FullIndexOffset)
	exception.PanicOnErr(err)
	err = binary.Read(file, binary.LittleEndian, &header.FieldCount)
	exception.PanicOnErr(err)
	if header.FieldCount != 0 {
		log.Fatalf("Error: bigWig header fieldCount field must be zero. Found: %v.\n", header.FieldCount)
	}
	err = binary.Read(file, binary.LittleEndian, &header.DefinedFieldCount)
	exception.PanicOnErr(err)
	if header.DefinedFieldCount != 0 {
		log.Fatalf("Error: bigWig header definedFieldCount field must be zero. Found: %v.\n", header.DefinedFieldCount)
	}
	err = binary.Read(file, binary.LittleEndian, &header.AutoSqlOffset)
	exception.PanicOnErr(err)
	if header.AutoSqlOffset != 0 {
		log.Fatalf("Error: bigWig header AutoSeqlOffset field must be zero. Found: %v.\n", header.AutoSqlOffset)
	}
	err = binary.Read(file, binary.LittleEndian, &header.TotalSummaryOffset)
	exception.PanicOnErr(err)
	err = binary.Read(file, binary.LittleEndian, &header.UncompressBufferSize)
	exception.PanicOnErr(err)
	err = binary.Read(file, binary.LittleEndian, &header.ExtensionOffset)
	exception.PanicOnErr(err)
	return header
}

func readZoomHeader(file *fileio.EasyReader) ZoomHeader {
	var answer = ZoomHeader{}
	var err error
	err = binary.Read(file, binary.LittleEndian, &answer.ReductionLevel)
	exception.PanicOnErr(err)
	err = binary.Read(file, binary.LittleEndian, &answer.Reserved)
	exception.PanicOnErr(err)
	err = binary.Read(file, binary.LittleEndian, &answer.DataOffset)
	exception.PanicOnErr(err)
	err = binary.Read(file, binary.LittleEndian, &answer.IndexOffset)
	exception.PanicOnErr(err)
	return answer
}

func readTotalSummary(file *fileio.EasyReader) TotalSummaryBlock {
	var answer = TotalSummaryBlock{}
	var err error
	err = binary.Read(file, binary.LittleEndian, &answer.BasesCovered)
	exception.PanicOnErr(err)
	err = binary.Read(file, binary.LittleEndian, &answer.MinVal)
	exception.PanicOnErr(err)
	err = binary.Read(file, binary.LittleEndian, &answer.MaxVal)
	exception.PanicOnErr(err)
	err = binary.Read(file, binary.LittleEndian, &answer.SumData)
	exception.PanicOnErr(err)
	err = binary.Read(file, binary.LittleEndian, &answer.SumSquares)
	exception.PanicOnErr(err)
	return answer
}

func readChromTreeHeader(file *fileio.EasyReader) ChromTreeHeader {
	var answer = ChromTreeHeader{}
	var err error
	err = binary.Read(file, binary.LittleEndian, &answer.Magic)
	exception.PanicOnErr(err)
	if answer.Magic == chromTreeMagicBigEndian {
		log.Fatalf("Error: bigWig package found a big endian chromosome tree header. Currently our package only supports little endian bigWig files.\n")
	}
	if answer.Magic != chromTreeMagic {
		log.Fatalf("Error: Expected to find chromosome tree magic number (2026540177) in file. Found: %v.\n", answer.Magic)
	}
	err = binary.Read(file, binary.LittleEndian, &answer.BlockSize)
	exception.PanicOnErr(err)
	err = binary.Read(file, binary.LittleEndian, &answer.KeySize)
	exception.PanicOnErr(err)
	err = binary.Read(file, binary.LittleEndian, &answer.ValSize)
	exception.PanicOnErr(err)
	err = binary.Read(file, binary.LittleEndian, &answer.ItemCount)
	exception.PanicOnErr(err)
	err = binary.Read(file, binary.LittleEndian, &answer.Reserved)
	exception.PanicOnErr(err)
	if answer.Reserved != 0 {
		log.Fatalf("Error: expected chromosome tree header reserved field to be equal to 0. Found: %v.\n", answer.Reserved)
	}
	return answer
}

func readChromTreeNodes(file *fileio.EasyReader, header ChromTreeHeader) []ChromTreeNode {
	var answer = make([]ChromTreeNode, 0)
	var currNode ChromTreeNode
	var currItem ChromTreeItem
	var currByteCounter uint32
	var err error
	var countIdx uint16
	numItems := 0

	for numItems < int(header.ItemCount) {
		// first we make a Node and read the node header
		currNode = ChromTreeNode{}
		err = binary.Read(file, binary.LittleEndian, &currNode.IsLeaf)
		exception.PanicOnErr(err)
		err = binary.Read(file, binary.LittleEndian, &currNode.Reserved)
		exception.PanicOnErr(err)
		err = binary.Read(file, binary.LittleEndian, &currNode.Count)
		exception.PanicOnErr(err)

		for countIdx = 0; countIdx < currNode.Count; countIdx++ {
			currItem = ChromTreeItem{}
			// first we parse the key (same for leaf and non-leaf nodes, and comes first)
			currItem.Key = make([]byte, header.KeySize)
			for currByteCounter = 0; currByteCounter < header.KeySize; currByteCounter++ {
				err = binary.Read(file, binary.LittleEndian, &currItem.Key[currByteCounter])
				exception.PanicOnErr(err)
			}
			if currNode.IsLeaf {
				err = binary.Read(file, binary.LittleEndian, &currItem.ChromId)
				exception.PanicOnErr(err)
				err = binary.Read(file, binary.LittleEndian, &currItem.ChromSize)
				exception.PanicOnErr(err)
			} else {
				err = binary.Read(file, binary.LittleEndian, &currItem.ChildOffset)
				exception.PanicOnErr(err)
			}
			currNode.Items = append(currNode.Items, currItem)
			numItems++
		}
		answer = append(answer, currNode)
	}

	return answer
}

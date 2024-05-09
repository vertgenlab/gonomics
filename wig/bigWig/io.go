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
	decodeBinaryField(file, &header.Magic)
	if header.Magic == bigWigMagicBigEndian {
		log.Fatalf("Error: bigWig file appears to be in big endian. Current functionality only supports little endian bigWig files.\n")
	}
	if header.Magic != bigWigMagic {
		log.Fatalf("Error: bigWig magic was not as expected. Found: %v. Expected: %v.\n", header.Magic, bigWigMagic)
	}
	decodeBinaryField(file, &header.Version)
	decodeBinaryField(file, &header.ZoomLevels)
	decodeBinaryField(file, &header.ChromosomeTreeOffset)
	decodeBinaryField(file, &header.FullDataOffset)
	decodeBinaryField(file, &header.FullIndexOffset)
	decodeBinaryField(file, &header.FieldCount)
	if header.FieldCount != 0 {
		log.Fatalf("Error: bigWig header fieldCount field must be zero. Found: %v.\n", header.FieldCount)
	}
	decodeBinaryField(file, &header.DefinedFieldCount)
	if header.DefinedFieldCount != 0 {
		log.Fatalf("Error: bigWig header definedFieldCount field must be zero. Found: %v.\n", header.DefinedFieldCount)
	}
	decodeBinaryField(file, &header.AutoSqlOffset)
	if header.AutoSqlOffset != 0 {
		log.Fatalf("Error: bigWig header AutoSeqlOffset field must be zero. Found: %v.\n", header.AutoSqlOffset)
	}
	decodeBinaryField(file, &header.TotalSummaryOffset)
	decodeBinaryField(file, &header.UncompressBufferSize)
	decodeBinaryField(file, &header.ExtensionOffset)
	return header
}

func readZoomHeader(file *fileio.EasyReader) ZoomHeader {
	var answer = ZoomHeader{}
	decodeBinaryField(file, &answer.ReductionLevel)
	decodeBinaryField(file, &answer.Reserved)
	decodeBinaryField(file, &answer.DataOffset)
	decodeBinaryField(file, &answer.IndexOffset)
	return answer
}

func readTotalSummary(file *fileio.EasyReader) TotalSummaryBlock {
	var answer = TotalSummaryBlock{}
	decodeBinaryField(file, &answer.BasesCovered)
	decodeBinaryField(file, &answer.MinVal)
	decodeBinaryField(file, &answer.MaxVal)
	decodeBinaryField(file, &answer.SumData)
	decodeBinaryField(file, &answer.SumSquares)
	return answer
}

func readChromTreeHeader(file *fileio.EasyReader) ChromTreeHeader {
	var answer = ChromTreeHeader{}
	decodeBinaryField(file, &answer.Magic)
	if answer.Magic == chromTreeMagicBigEndian {
		log.Fatalf("Error: bigWig package found a big endian chromosome tree header. Currently our package only supports little endian bigWig files.\n")
	}
	if answer.Magic != chromTreeMagic {
		log.Fatalf("Error: Expected to find chromosome tree magic number (2026540177) in file. Found: %v.\n", answer.Magic)
	}
	decodeBinaryField(file, &answer.BlockSize)
	decodeBinaryField(file, &answer.KeySize)
	decodeBinaryField(file, &answer.ValSize)
	decodeBinaryField(file, &answer.ItemCount)
	decodeBinaryField(file, &answer.Reserved)
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
	var countIdx uint16
	numItems := 0

	for numItems < int(header.ItemCount) {
		// first we make a Node and read the node header
		currNode = ChromTreeNode{}
		decodeBinaryField(file, &currNode.IsLeaf)
		decodeBinaryField(file, &currNode.Reserved)
		decodeBinaryField(file, &currNode.Count)

		for countIdx = 0; countIdx < currNode.Count; countIdx++ {
			currItem = ChromTreeItem{}
			// first we parse the key (same for leaf and non-leaf nodes, and comes first)
			currItem.Key = make([]byte, header.KeySize)
			for currByteCounter = 0; currByteCounter < header.KeySize; currByteCounter++ {
				decodeBinaryField(file, &currItem.Key[currByteCounter])
			}
			if currNode.IsLeaf {
				decodeBinaryField(file, &currItem.ChromId)
				decodeBinaryField(file, &currItem.ChromSize)
			} else {
				decodeBinaryField(file, &currItem.ChildOffset)
			}
			currNode.Items = append(currNode.Items, currItem)
			numItems++
		}
		answer = append(answer, currNode)
	}

	return answer
}

// decodeBinaryField reads binary data from an input *fileio.EasyReader to a variable.
func decodeBinaryField(file *fileio.EasyReader, data any) {
	err := binary.Read(file, binary.LittleEndian, data)
	exception.PanicOnErr(err)
}

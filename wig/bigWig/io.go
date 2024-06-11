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
	fileio.DecodeBinaryField(file, &header.Magic)
	if header.Magic == bigWigMagicBigEndian {
		log.Fatalf("Error: bigWig file appears to be in big endian. Current functionality only supports little endian bigWig files.\n")
	}
	if header.Magic != bigWigMagic {
		log.Fatalf("Error: bigWig magic was not as expected. Found: %v. Expected: %v.\n", header.Magic, bigWigMagic)
	}
	fileio.DecodeBinaryField(file, &header.Version)
	fileio.DecodeBinaryField(file, &header.ZoomLevels)
	fileio.DecodeBinaryField(file, &header.ChromosomeTreeOffset)
	fileio.DecodeBinaryField(file, &header.FullDataOffset)
	fileio.DecodeBinaryField(file, &header.FullIndexOffset)
	fileio.DecodeBinaryField(file, &header.FieldCount)
	if header.FieldCount != 0 {
		log.Fatalf("Error: bigWig header fieldCount field must be zero. Found: %v.\n", header.FieldCount)
	}
	fileio.DecodeBinaryField(file, &header.DefinedFieldCount)
	if header.DefinedFieldCount != 0 {
		log.Fatalf("Error: bigWig header definedFieldCount field must be zero. Found: %v.\n", header.DefinedFieldCount)
	}
	fileio.DecodeBinaryField(file, &header.AutoSqlOffset)
	if header.AutoSqlOffset != 0 {
		log.Fatalf("Error: bigWig header AutoSeqlOffset field must be zero. Found: %v.\n", header.AutoSqlOffset)
	}
	fileio.DecodeBinaryField(file, &header.TotalSummaryOffset)
	fileio.DecodeBinaryField(file, &header.UncompressBufferSize)
	fileio.DecodeBinaryField(file, &header.ExtensionOffset)
	return header
}

func readZoomHeader(file *fileio.EasyReader) ZoomHeader {
	var answer = ZoomHeader{}
	fileio.DecodeBinaryField(file, &answer.ReductionLevel)
	fileio.DecodeBinaryField(file, &answer.Reserved)
	fileio.DecodeBinaryField(file, &answer.DataOffset)
	fileio.DecodeBinaryField(file, &answer.IndexOffset)
	return answer
}

func readTotalSummary(file *fileio.EasyReader) TotalSummaryBlock {
	var answer = TotalSummaryBlock{}
	fileio.DecodeBinaryField(file, &answer.BasesCovered)
	fileio.DecodeBinaryField(file, &answer.MinVal)
	fileio.DecodeBinaryField(file, &answer.MaxVal)
	fileio.DecodeBinaryField(file, &answer.SumData)
	fileio.DecodeBinaryField(file, &answer.SumSquares)
	return answer
}

func readChromTreeHeader(file *fileio.EasyReader) ChromTreeHeader {
	var answer = ChromTreeHeader{}
	fileio.DecodeBinaryField(file, &answer.Magic)
	if answer.Magic == chromTreeMagicBigEndian {
		log.Fatalf("Error: bigWig package found a big endian chromosome tree header. Currently our package only supports little endian bigWig files.\n")
	}
	if answer.Magic != chromTreeMagic {
		log.Fatalf("Error: Expected to find chromosome tree magic number (2026540177) in file. Found: %v.\n", answer.Magic)
	}
	fileio.DecodeBinaryField(file, &answer.BlockSize)
	fileio.DecodeBinaryField(file, &answer.KeySize)
	fileio.DecodeBinaryField(file, &answer.ValSize)
	fileio.DecodeBinaryField(file, &answer.ItemCount)
	fileio.DecodeBinaryField(file, &answer.Reserved)
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
		fileio.DecodeBinaryField(file, &currNode.IsLeaf)
		fileio.DecodeBinaryField(file, &currNode.Reserved)
		fileio.DecodeBinaryField(file, &currNode.Count)

		for countIdx = 0; countIdx < currNode.Count; countIdx++ {
			currItem = ChromTreeItem{}
			// first we parse the key (same for leaf and non-leaf nodes, and comes first)
			currItem.Key = make([]byte, header.KeySize)
			for currByteCounter = 0; currByteCounter < header.KeySize; currByteCounter++ {
				fileio.DecodeBinaryField(file, &currItem.Key[currByteCounter])
			}
			if currNode.IsLeaf {
				fileio.DecodeBinaryField(file, &currItem.ChromId)
				fileio.DecodeBinaryField(file, &currItem.ChromSize)
			} else {
				fileio.DecodeBinaryField(file, &currItem.ChildOffset)
			}
			currNode.Items = append(currNode.Items, currItem)
			numItems++
		}
		answer = append(answer, currNode)
	}
	return answer
}

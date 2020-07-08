package bam

/*
import(
	"encoding/binary"
	"github.com/vertgenlab/gonomics/common"
	"log"
)
func bamLine(reader *BamReader) {
	var blockSize uint32 = 0
  	var flagNc    uint32 = 0
  	var binMqNl   uint32 = 0
  	for {
  		//reader.PipeIo.Data = make([]byte, blockSize)
  		reader.PipeIo.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &blockSize)
  		reader.PipeIo.Data = make([]byte, blockSize)

  		var i, j int = 0, 0
  		for i = 0; i < int(blockSize); {
  			bai := BinAln{}
  			reader.PipeIo.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &bai.RName)
  			common.ExitIfError(reader.PipeIo.Debug)

  			//position
  			reader.PipeIo.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &bai.Pos)
  			common.ExitIfError(reader.PipeIo.Debug)

  			reader.PipeIo.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &binMqNl)
  			common.ExitIfError(reader.PipeIo.Debug)

  			bai.Bai      = uint16((binMqNl >> 16) & 0xffff)
    		bai.MapQ     = uint8 ((binMqNl >>  8) & 0xff)
    		bai.RNLength = uint8 ((binMqNl >>  0) & 0xff)

    		reader.PipeIo.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &flagNc)
  			common.ExitIfError(reader.PipeIo.Debug)

  			bai.Flag = BamFlag(flagNc >> 16)
  			bai.NCigarOp = uint16(flagNc & 0xffff)

  			reader.PipeIo.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &bai.LSeq)
  			common.ExitIfError(reader.PipeIo.Debug)

  			reader.PipeIo.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &bai.NextRefID)
  			common.ExitIfError(reader.PipeIo.Debug)

  			reader.PipeIo.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &bai.NextPos)
  			common.ExitIfError(reader.PipeIo.Debug)



  			reader.PipeIo.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &bai.TLength)
  			common.ExitIfError(reader.PipeIo.Debug)



  			j, reader.PipeIo.Debug = reader.gunzip.Read(reader.PipeIo.Data[i:])
  			common.ExitIfError(reader.PipeIo.Debug)



  			log.Printf("%s\n", string(reader.PipeIo.Data))
  			i += j

  		}


  		//bai := BinAln{}
  		//ref name
  		//reader.PipeIo.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &bai.RName)
  		//common.ExitIfError(reader.PipeIo.Debug)

  		//position
  		//reader.PipeIo.Debug = binary.Read(reader.gunzip, binary.LittleEndian, &bai.Pos)
  		//common.ExitIfError(reader.PipeIo.Debug)




 	}
}*/

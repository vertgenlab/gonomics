package bam

//func bamBlockToSam()
/*
func bamTryAgain(reader *BamReader) {
	var blockSize int32
  	var flagNc    uint32
  	var stats   uint32
  	var err error
  	for {
  		block := BamData{}
  		buf := bytes.NewBuffer([]byte{})
  		// read block size
  		if err = binary.Read(reader.gunzip, binary.LittleEndian, &blockSize); err != nil {
  			if err == io.EOF {
  				return
      		}
      		common.ExitIfError(err)
  		}
  		//read block data
  		if err = binary.Read(reader.gunzip, binary.LittleEndian, &block.RName); err != nil {
  			common.ExitIfError(err)
  		}
  		if err = binary.Read(reader.gunzip, binary.LittleEndian, &block.Pos); err != nil {
      		common.ExitIfError(err)
      	}
      	if err = binary.Read(reader.gunzip, binary.LittleEndian, &stats); err != nil {
      		common.ExitIfError(err)
      	}
      	block.Bai      = uint16((stats >> 16) & 0xffff)
    	block.MapQ     = uint8 ((stats >>  8) & 0xff)
    	block.RNLength = uint8 ((stats >>  0) & 0xff)
    	if err = binary.Read(reader.gunzip, binary.LittleEndian, &flagNc); err != nil {
    		common.ExitIfError(err)
    	}
    	// get Flag and NCigarOp from FlagNc
    	block.Flag     = uint16(flagNc >> 16)
    	block.NCigarOp = uint16(flagNc & 0xffff)
    	if err = binary.Read(reader.gunzip, binary.LittleEndian, &block.LSeq); err != nil {
    		common.ExitIfError(err)
    	}
    	if err = binary.Read(reader.gunzip, binary.LittleEndian, &block.NextRefID); err != nil {
    		common.ExitIfError(err)
    	}
    	if err = binary.Read(reader.gunzip, binary.LittleEndian, &block.NextPos); err != nil {
    		common.ExitIfError(err)
    	}
    	if err = binary.Read(reader.gunzip, binary.LittleEndian, &block.TLength); err != nil {
    	 	common.ExitIfError(err)
    	}
    	// parse the read name
    	var b byte
    	for {
    		if err = binary.Read(reader.gunzip, binary.LittleEndian, &b); err != nil {
    			common.ExitIfError(err)

    		}
    		if b == 0 {
        		block.QName = buf.String()
        		break
      		}
      		buf.WriteByte(b)
    	}
    	var i int
    	// parse cigar block
    	block.Cigar = make(BamCigar, block.NCigarOp)
    	for i := 0; i < int(block.NCigarOp); i++ {
    		if err = binary.Read(reader.gunzip, binary.LittleEndian, &block.Cigar[i]); err != nil {
    			common.ExitIfError(err)
    		}
    	}
    	 // parse seq
    	block.Seq = make([]byte, (block.LSeq+1)/2)
    	for i = 0; i < int((block.LSeq+1)/2); i++ {
    		if err = binary.Read(reader.gunzip, binary.LittleEndian, &block.Seq[i]); err != nil {
    			common.ExitIfError(err)
    		}
    	}
    	block.Qual = make([]byte, block.LSeq)
    	for i = 0; i < int(block.LSeq); i++ {
    		if err = binary.Read(reader.gunzip, binary.LittleEndian, &block.Qual[i]); err != nil {
    			common.ExitIfError(err)
    		}
    	}
    	// read auxiliary data
    	position := 8*4 + int(block.RNLength) + 4*int(block.NCigarOp) + int((block.LSeq + 1)/2) + int(block.LSeq)
    	//for i := 0; position + i < int(blockSize); {

    	//}
    	if _, err = io.CopyN(ioutil.Discard, reader.gunzip, int64(blockSize) - int64(position)); err != nil {
    		common.ExitIfError(err)
    	}
    	log.Printf("Block Value: %v\n", block)
  	}
}*/

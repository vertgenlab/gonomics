package main

// TODO: Disabling tests until final program is complete.
/*
var StrawToBedPeTests = []struct {
	FileList string
	OutFile  string
	BinSize  int
	Expected string
}{
	{FileList: "testdata/fileList.txt",
		OutFile:  "testdata/out.bedpe",
		BinSize:  5000,
		Expected: "testdata/expected.out.bedpe"},
	{FileList: "testdata/fileList.txt",
		OutFile:  "testdata/out.interChrom.bedpe",
		BinSize:  5000,
		Expected: "testdata/expected.out.interChrom.bedpe"},
}

func TestStrawToBedpe(t *testing.T) {
	var err error
	for _, v := range StrawToBedPeTests {
		strawToBedpe(v.FileList, v.OutFile, v.BinSize)
		if !bedpe.AllAreEqual(bedpe.Read(v.Expected), bedpe.Read(v.OutFile)) {
			t.Errorf("outFile: %s did not match expected file: %s.", v.OutFile, v.Expected)
		} else {
			err = os.Remove(v.OutFile)
			exception.PanicOnErr(err)
		}
	}
}
*/

package fastq

import (
	"github.com/vertgenlab/gonomics/dna"
	"reflect"
	"testing"
)

func TestTrimmingBarcodes(t *testing.T) {
	var nameOneFwd string = "ST-E00545:471:HVJ32CCXY:2:1101:6705:1538_BX:NTCGGACGTGGGAAGG"
	var bxOne []dna.Base = dna.StringToBases("NTCGGACGTGGGAAGG")
	var umiOne []dna.Base = dna.StringToBases("AGGTGG")
	var seqOneFwd []dna.Base = dna.StringToBases("TGGCATGTTGGGCGGGTCAGACGTTGGCCATTGCTCCTCGAATGGTACTTGTAGTTGGACAATAAGCTCAGTTCCCCTTAAAAGTCGCAGTTCCCAAAACCCGGGATTCAAACGCGGGCGTGGTTTTG")
	var qualOneFwd []byte = ToQual([]byte("JAJAJFJ<F<-7FFF--7A-7-7-7F7-F-7-----7A-7--<AJ7-<F7A7F--<7AFJJ---F-7-----7-<<---7F<FAA--AF-77FA<FFFFF7-A<7-A-A-FJ7A--7-)7AA)7-FFF"))

	var nameTwoFwd string = "ST-E00545:471:HVJ32CCXY:2:1101:7740:1538_BX:NCCATGGCATCATGGT"
	var bxTwo []dna.Base = dna.StringToBases("NCCATGGCATCATGGT")
	var umiTwo []dna.Base = dna.StringToBases("AGATAT")
	var seqTwoFwd []dna.Base = dna.StringToBases("AGAGTTGTAGATGTAGAGTTGTAGATGCAGAGATGTAGAGTTGTAGATGCAGAGTTGTAGATGCAGAGTTGTAGAGATGTAGATGTAGATGTAGAGTTCTAGATGTAGATGCAGAGTTGTAGATGTAG")
	var qualTwoFwd []byte = ToQual([]byte("JJAFFJJJJJJJJJAFJFFFJJJFFJ-FA7JA-<FAJJF7--7FJJFJF777FA-<<FFJFJFFFFJJ-<JJJJ7AFFJ<AJJ<JFFFJ7FAJJAF---7FJF777-7A<A-JFAJ<---AFAFA7JA"))

	var nameOneRev string = "ST-E00545:471:HVJ32CCXY:2:1101:6705:1538_BX:NTCGGACGTGGGAAGG"
	var seqOneRev []dna.Base = dna.StringToBases("NACCTCAATAATGTAATCGACTTTCAGATTACGAGATAGAATGGTTAAGGCAATCTGCAGAGCTTTGTGGGACGGAAGAGGTATGTTACTTCAAAAGACCATCACAGGTAATCTGCAGCATCCCAACGAGCACAGAACCCGAAGGGCACA")
	var qualOneRev []byte = ToQual([]byte("#-AAAFFAF-F7--AFAA-77---FF--7---7<<-<---<---<------A---<7-<7AF<A-F--A--7----7--7AF-7--7-7-7-7-7A----7---7<-)---7-77-<)))--)7--7)--)))-77-7-))))7)))7-)"))

	var nameTwoRev string = "ST-E00545:471:HVJ32CCXY:2:1101:7740:1538_BX:NCCATGGCATCATGGT"
	var seqTwoRev []dna.Base = dna.StringToBases("NTACAACTCTACATCTACAACTCTACATCTCTGCATCTACAACTCTACATCTACAACTCTACATCTACAACTCTCCATCTACAACTCTCTACATCTACATCTATACATCAACAACTACATACAAATCTAACACTACAACTCTCCACAACT")
	var qualTwoRev []byte = ToQual([]byte("#AAAFFJJJJJJJJJJJJJJFJAFJFFJJJJJJFJFFAFAJJJJJF-AFJJJJJJFJFAJ7FJFAJ<7F-7FJF---77A-AF----A7-7---7-7FF-7-7---7<--<-------7<--7---<------7---7-----)7)----"))

	var linkedReadOne LinkedRead = LinkedRead{}
	linkedReadOne.Fwd = &Fastq{Name: nameOneFwd, Seq: seqOneFwd, Qual: qualOneFwd}
	linkedReadOne.Rev = &Fastq{Name: nameOneRev, Seq: seqOneRev, Qual: qualOneRev}
	linkedReadOne.Bx = bxOne
	linkedReadOne.Umi = umiOne
	var linkedReadTwo LinkedRead = LinkedRead{}
	linkedReadTwo.Fwd = &Fastq{Name: nameTwoFwd, Seq: seqTwoFwd, Qual: qualTwoFwd}
	linkedReadTwo.Rev = &Fastq{Name: nameTwoRev, Seq: seqTwoRev, Qual: qualTwoRev}
	linkedReadTwo.Bx = bxTwo
	linkedReadTwo.Umi = umiTwo

	var readPairNumber int = 0
	tenX := make(chan *LinkedRead)
	go ReadToChanLinked("testdata/barcode10x_R1.fastq", "testdata/barcode10x_R2.fastq", tenX)

	for fq := range tenX {
		if readPairNumber == 0 && !reflect.DeepEqual(fq, &linkedReadOne) {
			t.Errorf("Error: linked reads were not equal\n")
		}
		if readPairNumber == 1 && !reflect.DeepEqual(fq, &linkedReadTwo) {
			t.Errorf("Error: linked reads were not equal\n")
		}
		readPairNumber++
	}
}

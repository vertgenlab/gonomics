package bed

//moved to convert package

/*
import (
	"github.com/vertgenlab/gonomics/fasta"
	"log"
)

func toFasta(b *Bed, ref []*fasta.Fasta) *fasta.Fasta {
	for i := 0; i < len(ref); i++ {
		if b.Chrom == ref[i].Name {
			return fasta.Extract(ref[i], b.ChromStart, b.ChromEnd, b.Name)
		}
	}
	log.Fatalf("Chrom not found in fasta")
	return nil
}

func SliceToFasta(b []*Bed, ref []*fasta.Fasta) []*fasta.Fasta {
	outlist := make([]*fasta.Fasta, len(b))
	for i := 0; i < len(b); i++ {
		outlist[i] = toFasta(b[i], ref)
	}
	return outlist
}
*/

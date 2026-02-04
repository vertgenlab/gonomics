package reconstruct

import (
	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/fasta"
	"log"
)

func PostProb_intergenic_regions_extract(regions []bed.Bed, inputFa []fasta.Fasta, specifyChrom string) []fasta.Fasta{
	extracted := fasta.Fasta{}

	/// the input fasta file is X fastas, where X is the number of species in the alignment
	/// I need to initialise a X-length []fasta.Fasta{} output file
	/// and then iterate over the regions in the bed file
	/// extracting that region in each of the X sequences in the input file, 
	//  and append them onto each of the X-sequences in the output file 
	
	///////////////////////////////////// get only bed regions for specified chrom (if available), else return all regions
	filteredRegions := []bed.Bed{}
	if specifyChrom != "" {
		for _, r := range regions {
			if r.Chrom != specifyChrom {
				continue
			} else {
				filteredRegions = append(filteredRegions, r)
			}
		}
	} else {
		filteredRegions = regions
	}
	
	if len(filteredRegions) == 0 {
		log.Fatalf("Error: expecting more than 1 region in bed file")
	}

	///////////////////////////////////// get alignment positions
	var alnChromStart int
	var alnChromEnd int

	alnRegions := make([]bed.Bed, len(filteredRegions))

	//// first region (idx 0)
	alnStart := fasta.RefPosToAlnPos(inputFa[0], filteredRegions[0].ChromStart)
	alnEnd := fasta.RefPosToAlnPos(inputFa[0], filteredRegions[0].ChromEnd)

	alnRegions[0].ChromStart = alnStart
	alnRegions[0].ChromEnd = alnEnd
	alnRegions[0].Chrom = filteredRegions[0].Chrom

	refStart := filteredRegions[0].ChromEnd
	alnStart = alnEnd

	for idx, r := range filteredRegions[1:len(filteredRegions)-1] {
		alnChromStart = fasta.RefPosToAlnPosCounter(inputFa[0], r.ChromStart, refStart, alnStart)
		alnChromEnd = fasta.RefPosToAlnPosCounter(inputFa[0], r.ChromEnd, r.ChromStart, alnChromStart)
		alnRegions[idx+1].Chrom = r.Chrom
		alnRegions[idx+1].ChromStart = alnChromStart
		alnRegions[idx+1].ChromEnd = alnChromEnd

		refStart = r.ChromEnd
		alnStart = alnChromEnd
	}

	//// last region
	region := filteredRegions[len(filteredRegions)-1]
	alnChromStart = fasta.RefPosToAlnPosCounter(inputFa[0], region.ChromStart, refStart, alnStart)

	alnChromEnd = fasta.RefPosToAlnPosCounter(inputFa[0], region.ChromEnd, region.ChromStart, alnChromStart)

	alnRegions[len(filteredRegions)-1].Chrom = region.Chrom
	alnRegions[len(filteredRegions)-1].ChromStart = alnChromStart
	alnRegions[len(filteredRegions)-1].ChromEnd = alnChromEnd


	///////////////////////////////////// extract from multifa based on alignment bed
	ans := make([]fasta.Fasta, len(inputFa))
	for idx, inSeq := range inputFa {
		ans[idx].Name = inSeq.Name

		for _, r := range alnRegions {
			extracted = fasta.Extract(inSeq, r.ChromStart, r.ChromEnd, inSeq.Name)
			ans[idx].Seq = append(ans[idx].Seq, extracted.Seq...)
		}
	}
	return ans
}
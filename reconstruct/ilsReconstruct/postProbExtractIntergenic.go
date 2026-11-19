package ilsReconstruct

import (
	"log"

	"github.com/vertgenlab/gonomics/bed"
	"github.com/vertgenlab/gonomics/fasta"
)

// concatenateMultifa concatenates concatenates record-wise each multifa in a list of multifas
func concatenateMultifa(inputMultifa [][]fasta.Fasta) []fasta.Fasta {
	ans := make([]fasta.Fasta, len(inputMultifa))
	for i_idx, multi := range inputMultifa {
		ans[i_idx].Name = multi[0].Name
		for _, seq := range multi {
			ans[i_idx].Seq = append(ans[i_idx].Seq, seq.Seq...)
		}
	}
	return ans
}

// PostProbExtractIntergenic takes as input regions, sorted bed file,
// and inputFa, a multifasta, where all sequences are aligned to the first sequence (index 0),
// then extracts all regions in 'regions' from inputFa and output them as a list of
func PostProbExtractIntergenic(regions []bed.Bed, inputFa []fasta.Fasta, specifyChrom string) ([][]fasta.Fasta, []fasta.Fasta) {
	extracted := fasta.Fasta{}

	/// the input fasta file is X fastas, where X is the number of species in the alignment
	/// I need to initialise a X-length []fasta.Fasta{} output file
	/// and then iterate over the regions in the bed file
	/// extracting that region in each of the X sequences in the input file,
	//  and append them onto each of the X-sequences in the output file

	// Retain only BED regions matching the specified chromosome.
	// If no chromosome is specified, use all regions.
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

	//// Map first region directly from reference to alignment coordinates.
	if len(inputFa) == 0 {
		log.Fatalf("Must have at least 1 FASTA record.")
	}
	alnStart := fasta.RefPosToAlnPos(inputFa[0], filteredRegions[0].ChromStart)
	alnEnd := fasta.RefPosToAlnPos(inputFa[0], filteredRegions[0].ChromEnd)

	alnRegions[0].ChromStart = alnStart
	alnRegions[0].ChromEnd = alnEnd
	alnRegions[0].Chrom = filteredRegions[0].Chrom

	refStart := filteredRegions[0].ChromEnd
	alnStart = alnEnd

	//// Map intermediate (between 1st and last) regions using the previous reference and alignment
	//// endpoints to avoid repeatedly scanning from the start of the alignment
	for idx, r := range filteredRegions[1 : len(filteredRegions)-1] {
		// Convert BED coordinates (relative to the ungapped reference) to multiple sequence
		// alignment coordinates, accounting for gaps in the reference sequence
		alnChromStart = fasta.RefPosToAlnPosCounterBed(inputFa[0], r.ChromStart, refStart, alnStart)
		alnChromEnd = fasta.RefPosToAlnPosCounterBed(inputFa[0], r.ChromEnd, r.ChromStart, alnChromStart)
		alnRegions[idx+1].Chrom = r.Chrom
		alnRegions[idx+1].ChromStart = alnChromStart
		alnRegions[idx+1].ChromEnd = alnChromEnd

		refStart = r.ChromEnd
		alnStart = alnChromEnd
	}

	//// last region
	ans := make([][]fasta.Fasta, len(inputFa))
	if len(filteredRegions) > 1 {
		region := filteredRegions[len(filteredRegions)-1]
		alnChromStart = fasta.RefPosToAlnPosCounterBed(inputFa[0], region.ChromStart, refStart, alnStart)

		alnChromEnd = fasta.RefPosToAlnPosCounterBed(inputFa[0], region.ChromEnd, region.ChromStart, alnChromStart)

		alnRegions[len(filteredRegions)-1].Chrom = region.Chrom
		alnRegions[len(filteredRegions)-1].ChromStart = alnChromStart
		alnRegions[len(filteredRegions)-1].ChromEnd = alnChromEnd

		// Extract each region from every species using the same alignment coords
		for idx, rec := range inputFa {
			for _, r := range alnRegions {
				extracted = fasta.Extract(rec, r.ChromStart, r.ChromEnd, rec.Name)
				ans[idx] = append(ans[idx], extracted)
			}
		}
	} else { // filteredRegions == 1
		// If only one region selected, extract directly for each species
		for idx, rec := range inputFa {
			firstRegion := alnRegions[0]
			ans[idx] = append(ans[idx], fasta.Extract(rec, firstRegion.ChromStart, firstRegion.ChromEnd, rec.Name))
		}
	}
	return ans, concatenateMultifa(ans)
}

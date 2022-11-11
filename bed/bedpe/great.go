package bedpe

//input: bedpe, our tssBed
//within program: fillSpaceNoHiddenValues
//calculate hidden values for all bedpes
//fillSpaceWithHiddenValues
//output: phantomGenes and properGenes beds, spaceFilled

//func Fill3dSpace(contacts []BedPe, tss []bed.Bed, sizes map[string]chromInfo.ChromInfo) []bed.Bed {
//	var answer []bed.Bed
//	var closest2dGenesIntervals []interval.Interval
//	closest2dGene := bed.FillSpaceNoHiddenValue(tss, sizes)
//
//	for i := range closest2dGene {
//		closest2dGenesIntervals = append(closest2dGenesIntervals, closest2dGene[i])
//	}
//	closest2dGeneTree := interval.BuildTree(closest2dGenesIntervals)
//
//	return answer
//}

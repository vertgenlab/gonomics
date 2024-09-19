[1mdiff --git a/simulate/bed.go b/simulate/bed.go[m
[1mindex 22f7a3cf..d2dcf036 100644[m
[1m--- a/simulate/bed.go[m
[1m+++ b/simulate/bed.go[m
[36m@@ -42,14 +42,14 @@[m [mfunc GenerateBedRegion(searchSpace []bed.Bed, randPos int, regionLength int) (be[m
 [m
 		// Decrement randomly generated overall position (corresponds to start of generated region) until it fits within a region[m
 		// must have randPos < chromWindows (randPos is 0-indexed, chromWindows is not), at most randPos + 1 = chromWindows[m
[31m-		if randPos > chromWindows {[m
[32m+[m		[32mif randPos - chromWindows > -1 {[m
 			randPos -= chromWindows[m
 		} else {[m
 			log.Print("generated\n\n\n")[m
 			return bed.Bed{[m
 				Chrom: searchSpace[j].Chrom, [m
[31m-				ChromStart: searchSpace[j].ChromStart + randPos - 1, [m
[31m-				ChromEnd: searchSpace[j].ChromStart + randPos - 1 + regionLength, [m
[32m+[m				[32mChromStart: searchSpace[j].ChromStart + randPos,[m[41m [m
[32m+[m				[32mChromEnd: searchSpace[j].ChromStart + randPos + regionLength,[m[41m [m
 				Name: searchSpace[j].Name,[m
 				FieldsInitialized: 4}, true[m
 		}[m

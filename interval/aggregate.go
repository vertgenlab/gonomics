package interval

import "github.com/vertgenlab/gonomics/fileio"

type AggregateInterval struct {
	chr 		string
	start 		int
	end 		int
	components 	[]Interval
}

func (a *AggregateInterval) GetChrom() string {
	return a.chr
}
func (a *AggregateInterval) GetChromStart() int {
	return a.start
}
func (a *AggregateInterval) GetChromEnd() int {
	return a.end
}
func (a *AggregateInterval) SetExclude() {
	a.chr = "EXCLUDE"
}
func (a *AggregateInterval) WriteToFileHandle(file *fileio.EasyWriter) {
	for _, val := range a.components {
		val.WriteToFileHandle(file)
	}
}

func MergeIntervals(intervals []Interval) []Interval {
	m := splitIntervalsByChr(intervals)
	answer := make([]Interval, 0, len(intervals))
	for _, input := range m {
		sortIntervals(input, xLess)
		var curr *AggregateInterval
		for i := 0; i < len(input); i++ {

			if i + 1 < len(input) - 1 && input[i+1].GetChromStart() <= input[i].GetChromEnd() {
				curr = &AggregateInterval{
					chr: 	input[i].GetChrom(),
					start: 	input[i].GetChromStart(),
					end: 	input[i].GetChromEnd(),
					components: []Interval{input[i]},
				}
				for i++; curr.GetChrom() == input[i].GetChrom() &&
					input[i].GetChromStart() <= curr.GetChromEnd(); i++ {
					if curr.GetChromEnd() < input[i].GetChromEnd() {
						curr.end = input[i].GetChromEnd()
					}
					curr.components = append(curr.components, input[i])
					if i + 1 > len(input) - 1 {
						answer = append(answer, curr)
						return answer
					}
				}
				answer = append(answer, curr)
			}
			answer = append(answer, input[i])
		}
	}
	return answer
}
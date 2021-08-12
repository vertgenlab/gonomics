#!/bin/csh -ef

~/go/bin/multiFaAcceleration -windowSize 50 chr1 test.fa test.vel.expected.bed test.accel.expected.bed test.initialVel.expected.bed
~/go/bin/multiFaAcceleration -windowSize 50 -searchSpaceBed test.searchspace.bed chr1 test.fa test.vel.searchspace.expected.bed test.accel.searchspace.expected.bed test.initialVel.searchspace.expected.bed
~/go/bin/multiFaAcceleration -windowSize 50 -searchSpaceBed test.searchspace.bed -useSnpDistance chr1 test.fa test.vel.snpDistance.expected.bed test.accel.snpDistance.expected.bed test.initialVel.snpDistance.expected.bed

echo DONE

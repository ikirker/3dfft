#!/bin/bash
# This script gathers the data produced by the various runs and
#  turns it into a collected results file and a few specific
#  comparison data sets.
#
# ian.kirker@gmail.com
#

# The C code that produces the CSV output.
<< EOF 
printf("fft-results:%d,%d,%s,%s,%s,%g,%g,%g\n",
			size,
			extent,
			((decompType==1)?"slab":"rod"),
			((use2DFFT==1)?"2DFFT":"1DFFT"),
			FFT_NAME,
			
			/* Communication/memory reorg time */ 
			(phaseTime[2] - phaseTime[1]) + 
			 (phaseTime[4] - phaseTime[3]),
			
			/* FFT time */
			(phaseTime[1] - phaseTime[0]) + 
			 (phaseTime[3] - phaseTime[2])+
			 (phaseTime[5] - phaseTime[4]),
			 
			 /* Total time */
			phaseTime[5] - phaseTime[0]
			);
EOF

# So, size, extent, decomp, use2DFFT, lib name, reorg time, fft time, total time
hostname >> all_data.csv
grep -h "fft-results:" $@ | sed 's/fft-results://' >> all_data.csv

# Example record
# 32,64,slab,1DFFT,fftw3,1.5,2.5,6

# cat all_data.csv | dbInsert.pl


# Low tech - replaced by db pl scripts.
#awk -F"," '{ OFS=","; if ( ($4 != "2DFFT") && ($5 == "fftw3") ){ print $1,$2,$3,$6,$7,$8 } }' all_data.csv > fftw3-decomp-comparison.csv
#awk -F"," '{ OFS=","; if ( ($4 != "2DFFT") && ($5 == "fftw2") ){ print $1,$2,$3,$6,$7,$8 } }' all_data.csv > fftw2-decomp-comparison.csv
#awk -F"," '{ OFS=","; if ( ($4 != "2DFFT") && ($5 == "amcl" ) ){ print $1,$2,$3,$6,$7,$8 } }' all_data.csv > amcl-decomp-comparison.csv
#awk -F"," '{ OFS=","; if ( ($4 != "2DFFT") && ($5 == "mkl"  ) ){ print $1,$2,$3,$6,$7,$8 } }' all_data.csv > mkl-decomp-comparison.csv
#awk -F"," '{ OFS=","; if ( ($4 != "2DFFT") && ($5 == "essl" ) ){ print $1,$2,$3,$6,$7,$8 } }' all_data.csv > essl-decomp-comparison.csv
#awk -F"," '{ OFS=","; if ( ($3 == "rod"  )                    ){ print $1,$2,$5,$6,$7,$8 } }' all_data.csv > 1D-over-libraries-comparison.csv
#awk -F"," '{ OFS=","; if ( ($3 == "slab" )                    ){ print $1,$2,$4,$5,$6,$7,$8 } }' all_data.csv > 2D-over-libraries-comparison.csv


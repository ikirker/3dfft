#!/bin/bash
# This script gathers together the information required to produce the
#   script files wanted, asks the user whether to produce them, then
#   produces them using a template file supplied. It produces two
#   temporary files, header and perexec, and deletes them at the end.
# It relies on sed, cat and bc.
#
# ian.kirker@gmail.com

#####
# First, get all the info we need.
#####

#Choose the job script template:

## Gives numbered aliases for each template in the templates dir.
let loopIndex=1
for template in ${0%batchmaker.sh}/templates/*
do
	# Loop over the template files, assigning each to a variable
	#  of the form $template1,$template2 etc.
	varName="template$loopIndex"
	eval $"template$loopIndex"=$template
	echo "${loopIndex}) $template" 
	let loopIndex=loopIndex+1
done	

echo -n "Select one of the above templates for the batch files. > "
read batchSelection
## Read in a number, append it to "$template" and evaluate the variable produced.
eval template=\$"template$batchSelection"

echo "$template selected."


# Lower limit of processor count
echo -n "Enter the minimum number of processors to use. [32] > "
read minNumProcs

if [ "$minNumProcs" == "" ]
then let minNumProcs=32
fi

# Upper limit of processor count
echo -n "Enter the maximum number of processors to use. [1024] > "
read maxNumProcs

if [ "$maxNumProcs"  == "" ]
then let maxNumProcs=1024
fi

# Lower limit of extents
echo -n "Enter the smallest extent to use. [32] > "
read minExtent

if [ "$minExtent"  == "" ]
then let minExtent=32
fi

# Upper limit of extents
echo -n "Enter the largest extent to use. [1024] > "
read maxExtent

if [ "$maxExtent"  == "" ]
then let maxExtent=1024
fi

# Approximate RAM per processor
echo "Enter the approximate amount of RAM per processor in megabytes [2048] "
echo " ( HECToR - 3072MB ) "
echo " ( HPCx - 2048MB ) "
echo " ( Ness - 2048MB ) "
echo " ( BlueGene - 256MB ) "
echo " ( Eddie - 2048MB ) "
echo -n " > "
read ram

if [ "$ram"  == "" ]
then let ram=2048
fi

######
# I think that's everything, so check details, compile, and then start making batch scripts
######

# Output details for verification

echo
echo "You have provided the following details:"
echo "----------------------------------------"
echo "Batch Template: $template"
echo "Processors:     $minNumProcs - $maxNumProcs"
echo "Extent range:   $minExtent - $maxExtent"
echo

# Function to check whether 1st argument number is a power of two
IsPowerOfTwo()
{
	#Fixed point arithmetic
	a=$1
	let e=a*10000
	while (( $e>10000 ))
	do
		let e=e/2
	done
	if (( e==10000 ))
	then 
		return 1
	else
		return 0
	fi
}

GetDecompDim()
{ # This returns the smaller of the two numbers of processors
  #  that makes a 2D decomposition.
bc <<HEREDOC
n=$1
for( i=(sqrt(n)); i>0; i-- ) {
  if ((i % 2) == 0) if ((n % i) == 0) {
    o=i
    break
  }
}
o=n/o
o
HEREDOC
}
	

MakeFilteredOutput()
{ # Uses temporary files created by template, replacing variables
let halfcpus=$2/2
cat "$1" | sed \
-e "s:%^HALFCPUS%^:${halfcpus}:g" \
-e "s:%^CPUS%^:$2:g" \
-e "s:%^EXTENT%^:$3:g" \
-e "s:%^EXEC%^:$4:g" \
-e "s:%^DECOMP%^:$5:g" \
-e "s:%^COUNT%^:$6:g" >> $7
}

#Function for generating job scripts with one header and multiple execs.
MakeJobs()
{
	template=$1
	version=$2
	minNumProcs=$3
	maxNumProcs=$4
	minExtent=$5
	maxExtent=$6
	ram=$7
	
	let p=minNumProcs
	let x=minExtent
	let adding=$minExtent/4
	
	let jobfileCount=1
	
	# Creates temporary files named header and perexec from the template file
	source $template
	
	while (( $p <= $maxNumProcs ))
	do
		let execCount=1
		MakeFilteredOutput "header" $p $x $version 1 0 job-$version-$p.nys		
	
		while (( $x <= $maxExtent ))
		do 
			#Check we're not using too much memory
			memPerProc=`echo "(5 * ( $x ^ 3 ) * 2 * 8)/(2 * $p * (1024^2) )"|bc`
			
			if (( ram>memPerProc )); then
				IsPowerOfTwo $x
				if [[ "$?" == "1" ]]
				then
					if (( x>=p )) # Check we can do a slab decomposition
					then
						# Slab decomp
						MakeFilteredOutput "perexec" $p $x $version 1 $execCount job-$version-$p.nys
						let execCount=$execCount+1
						
						#Slab decomp w/ 2D FFT on slabs
						MakeFilteredOutput "perexec" $p $x $version 3 $execCount job-$version-$p.nys
						let execCount=$execCount+1
				
						#Auto decomp
						if [ "$version" == "fft-fftw2" ]
						then
							MakeFilteredOutput "perexec" $p $x $version 0 $execCount job-$version-$p.nys
							let execCount=$execCount+1
						fi
						if [ "$version" == "fft-essl" ]
						then
							MakeFilteredOutput "perexec" $p $x $version 0 $execCount job-$version-$p.nys
							let execCount=$execCount+1
						fi
					fi
			
					# Update increment so that 3 data points
					#  are obtained between each power of two.
					let adding=x/4 
				fi
				
				# Rod decomp is impossible for p < 4
				if (( p>=4 ))
				then
					# Check x is a multiple of the rod decomposition dims and
					#  is more than the dims multiplied
					dim=$(GetDecompDim $p)
					let odim=p/dim
					# This is for debugging the decomposition algorithm
					#echo "dim: $dim, odim: $odim, p: $p"
					if (( ( ( x/dim ) * dim == x) && (p==(dim*odim) ) ))
					then
						#Rod decomp
						MakeFilteredOutput "perexec" $p $x $version 2 $execCount job-$version-$p.nys
						let execCount=$execCount+1
					fi
				fi
			fi #RAM check
			let x=x+adding
		done # Loop over x
		let p=p*2
		let x=$minExtent
	done # Loop over p
	
	rm header perexec # Cleanup.
}


# Generate set of batch scripts for a run.
echo "Would you like to generate batch scripts now? (y/N)"
read yn

if [ "$yn" == "y" ]
then
	let p=minNumProcs
	let x=minExtent
	let adding=$minExtent/4
	
	# If no arguments were provided, make scripts for fftw3 - default version, otherwise, make one for each.
	if [ "$1" == "" ]
	then
		echo "No versions specified, using default: fft-fftw3"
		echo "Executables to use must be specified as $0 exec1 exec2 exec3 ..."
		versions="fft-fftw3"
	else
		versions="$1 $2 $3 $4" # There shouldn't be more than 4 versions present.
	fi
	
	for version in ${versions}
	do
		MakeJobs $template $version $minNumProcs $maxNumProcs $minExtent $maxExtent $ram
	done
fi

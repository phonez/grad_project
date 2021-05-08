#!/bin/bash

# extract_scores.bash score.sc

awk 'BEGIN{OFS=","}{
	
		if( NR != 1){
			if( $2 == "total_score" ){
				for(x=1; x<= NF; ++x){
					if($x == "description"){
						file=x
					}
					if($x == "total_score"){
						total=x
					}
					if($x == "rmsd"){
						rmsd=x
					}
				} 
			}else {
				print $file, $total, $rmsd
			}
		}
	}' $1


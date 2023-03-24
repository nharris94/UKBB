samtools depth 6018697_23143_0_0.cram   | awk '{sum+=$3} END { print "Average = ",sum/NR}' 

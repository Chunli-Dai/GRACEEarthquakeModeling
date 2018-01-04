#! /bin/bash 

if [ $# -lt 3 ]
then
  echo "usage: $0 filelist renameodir renameofile"
  exit
fi

filelist=$1
renameodir=$2
renameofile=$3

rm -f $renameofile

for line in ` cat $filelist `
do
  #GSM-2_2012336-2012366_0029_EIGEN_G---_0005
  #2003_016_0.txt
# echo $line
  year1=${line:6:4}
  doy1=${line:10:3}
  year2=${line:14:4}
  doy2=${line:18:3}
# echo $year1 $doy1 $year2 $doy2
# day=($doy1+$doy2)
  oneyr=365
  remainder=`expr $year1 % 4`   # remainder=$(($year1%4))
  if [ $remainder -eq 0 ]; then
#   oneyr=366;                  # Should use this line, to be consistent, use the next
    doy2=`echo "$doy2+1" | bc`  # Mistake in Guo Lao-shi's code
  fi
  if [ $doy2 -lt $doy1 ]; then
  doy2=`echo "$doy2+$oneyr" | bc`
  fi
# day=($doy1+$doy2)/2
#  day=`expr $doy1 + $doy2`
  day=`echo "($doy1 + $doy2 )/2." | bc`
# echo $day

  if [ $day -lt 10 ]; then
    doy=00"$day"
  elif [ $day -lt 100 ]; then
    doy=0"$day"
  else
    doy="$day"
  fi
  ofile=$year1'_'$doy'_0.txt'
# echo $ofile >> $renameofile
  echo $year1'_'$doy >> $renameofile
  mv $line $renameodir/$ofile
done

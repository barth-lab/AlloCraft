# This script reads xtc trajectories and transforms them to dcd
# Read the directory, number of runs, and stride to use for the trajs
set dirname [lindex $argv 0]
set nruns [lindex $argv 1]
set stride [lindex $argv 2]
set pdbName [lindex $argv 3]
set xtcName [lindex $argv 4]

for {set i 1} {$i <= $nruns} {incr i 1} {
    mol new $pdbName
    mol rename top "Run$i"
    mol addfile $dirname/run$i/$xtcName step $stride waitfor all
    mol modcolor 0 top "ColorID" $i
#    animate delete  beg 0 end 0 skip 0 top
#    mol addrep top
#    mol representation licorice
#    mol selection resname AMAN BMAN BGLC
#    mol material Diffuse
#    mol modrep 1 top

    set sel [atomselect top "all"]
    animate write dcd $dirname/run$i/traj.dcd beg 1 end -1 sel $sel
    mol off top
}
exit

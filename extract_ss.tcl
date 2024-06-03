# This script takes hte input pdb file and outputs the secondary structure to a file
set protpath [lindex $argv 0]
set outpath [lindex $argv 1]
mol new $protpath
set outfile [open $outpath w]
set lookup {H G I}
set frame_num [molinfo top get numframes]
set full [atomselect top "name CA"]
set len [llength [$full get resid]]
$full delete

for {set i 0} {$i < $frame_num} {incr i} {
    animate goto $i
    set sel [atomselect top "name CA"]
    mol ssrecalc top
    set struc_string [$sel get structure]
    set helix 0
   # puts $outfile "$i\t$struc_string"
    puts $outfile $struc_string
    $sel delete
}
close $outfile
exit

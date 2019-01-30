display projection orthographic     
display depthcue off     
axes location off
color Display Background white
#set target [lindex $argv 0]
set name [lindex $argv 0]
set drugname [lindex $argv 1]
set drugu [string toupper $drugname]
set aaregion_str $::env(aaregion_str)

#mol new er-$target.pdb
mol new $name.pdb
mol modselect 0 0 "noh and not ($aaregion_str)"
mol modstyle 0 0 Licorice 0.1 35.0 35.0
mol modmaterial 0 0 Transparent
mol addrep 0
#David Koes suggests having helix 12 in orange.
mol modselect 1 0 "noh and (resid 228 to 243)"
mol modstyle 1 0 Licorice 0.4 35.0 35.0
mol modcolor 1 0 ColorID 3
mol addrep 0
mol modselect 2 0 "noh and ($aaregion_str) and not (resname $drugu)"
mol modstyle 2 0 Licorice 0.4 35.0 35.0
mol addrep 0
mol modselect 3 0 "name CA and not ($aaregion_str)"
mol modstyle 3 0 VDW 0.25 35.0
mol modcolor 3 0 ColorID 7
mol addrep 0
mol modselect 4 0 "noh and resname $drugu"
mol modstyle 4 0 Licorice 0.4 35.0 35.0
mol modcolor 4 0 ColorID 11


set caatoms [atomselect 0 "name CA and not ( $aaregion_str )"]
set caindices [$caatoms get index]
set caxyz [$caatoms get {x y z}]
set n [llength $caindices]
#draw go model bonds
draw color yellow
for {set i 0} {$i < $n} {incr i} {
	for {set j [expr $i + 2]} {$j < $n} {incr j} {
		set iatom [lindex $caindices $i]
		set jatom [lindex $caindices $j]
		set d [measure bond [list $iatom $jatom]]
		if { $d <= 8.0 } {
			puts "$iatom $jatom $d"
			set ixyz [lindex $caxyz $i]
			set jxyz [lindex $caxyz $j]
			puts "$ixyz $jxyz"
			draw cylinder $ixyz $jxyz radius 0.05 resolution 25 filled no
		}
	}
}
display resetview

#display resize 1338 1668
display resize 5352 6672

rotate z by -90
rotate y by 45
rotate x by -30 
#translate by 0.0 0.0 0.5
scale by 1.75

render TachyonInternal "go-model-$name.tga"
exit
quit


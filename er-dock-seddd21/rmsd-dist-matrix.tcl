proc read_equivalence_info {fname} {
        set equivalences [list]
        set current [list]
        set infile [open $fname "r"]
        while { [gets $infile line] >= 0 } {
                scan $line "%s %s" word1 word2
                #puts "$word1 $word2"
                if { [string compare $word1 "EQUIV"] == 0 } {
                        scan $word2 "%d" id
                        #puts "$junk $id"
                        if { $id > 1} {
                                lappend equivalences $current
                        }
                        set current [list]
                } else {
                        #scan $line "%s %s" name1 name2
                        set entry [list $word1 $word2]
                        lappend current $entry
                }
        }
        lappend equivalences $current
        #puts $equivalences
        return $equivalences
}

#atomlist is a list giving atom numbers by name: {C1 1} etc.  To be obtained from atomselect object
#equiv is a list giving equivalences {C1 C2}.  first atom is in original molecule, second is equivalent atom in molecule
proc calc_permuted_rmsd {coords ref alisttraj alistref equiv} {
        set rmsd 0.0
        set n [llength $equiv]
	#puts "calc_permuted_rmsd: $alisttraj $alistref"
        foreach pair $equiv {
                set origname [lindex $pair 0]
                set eqname [lindex $pair 1]
                #find atom numbers -- thanks to stackoverflow for this
                set iorig [lsearch -index 0 -all -inline -exact $alistref $origname]
		set iorig [lindex $iorig 0 1]
                set ieq [lsearch -index 0 -all -inline -exact $alisttraj $eqname] 
		set ieq [lindex $ieq 0 1]
                #puts "calc_permuted_rmsd: $origname $iorig $eqname $ieq"
                set origxyz [lindex $ref $iorig]
                set eqxyz [lindex $coords $ieq]
		#puts "calc_permuted_rmsd: $origxyz $eqxyz"
                set disp [vecsub $eqxyz $origxyz]
                set d2 [vecdot $disp $disp]
                set rmsd [expr $rmsd + $d2]
        }
        set rmsd [ expr sqrt( $rmsd / $n )]
}


set target [lindex $argv 0]
set drugname [lindex $argv 1]
#set ref [lindex $argv 2]
set eqfname [lindex $argv 2]
set pdb [lindex $argv 3]
#set start [lindex $argv 4]
#set end [lindex $argv 5]
#set nframes [lindex $argv 6]
#set nfm1 [expr $nfile - 1]
set outfname [lindex $argv 4]
set drugu [string toupper $drugname]
set aaregion_str $::env(aaregion_str)
set aaregion_str "($aaregion_str) or (noh and resname $drugu)" 

set equivalences [read_equivalence_info $eqfname]

mol new $pdb waitfor all
set nframes [molinfo top get numframes]
set nframesm1 [expr $nframes - 1]

#load each trajectory's ending into a separate mol, to account for differences in structure
set backboneref [atomselect 0 "protein and (name N or name CA or name C)" ] 
$backboneref num
for {set i 0} {$i < $nframes} {incr i} {
	puts "aligning $i"
	animate goto $i
	set backbone [atomselect 0 "protein and (name N or name CA or name C)"] 
	set allatoms [atomselect 0 "all"]
	set matrix [measure fit $backbone $backboneref]
        $allatoms move $matrix
}
#measure the RMSD for each frame.
set output [open $outfname "w"]
#measure pairwise rmsds.
for {set i 0} {$i < $nframes} {incr i} {
	set refcoords [[atomselect 0 "all" frame $i] get {x y z}];
	set drug1 [atomselect 0 "noh and resname $drugu"]
	set alistdrug1 [$drug1 get {name index}]
	for {set j $i} {$j < $nframes} {incr j} {
		puts "doing rmsd $i $j"
		set coords [[atomselect 0 "all" frame $j] get {x y z}];
		set minrmsd 1000.0
		set drug2 [atomselect 0 "noh and resname $drugu"]
	        set alistdrug2 [$drug2 get {name index}]
		foreach equiv $equivalences {
                	#puts "pre calc_permuted_rmsd: $alistdrug $alistref"
                	set rmsd [calc_permuted_rmsd $coords $refcoords $alistdrug1 $alistdrug2 $equiv]
                	if { $rmsd < $minrmsd } {
                        	set minrmsd $rmsd
                	}
        	}
		set aaregion1 [atomselect 0 "$aaregion_str" frame $i]
		set aaregion2 [atomselect 0 "$aaregion_str" frame $j]
		set aaregion_rmsd [measure rmsd $aaregion1 $aaregion2] 
		set id1 [expr $i + 1]
		set id2 [expr $j + 1]
		#set id1 [lindex $idslist $i]
		#set id2 [lindex $idslist $j]
		puts $output "$id1 $id2 $minrmsd $aaregion_rmsd"
	}
}
close $output
exit


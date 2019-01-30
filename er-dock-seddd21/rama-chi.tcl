proc getphipsi {mol res} {
        set prev [expr $res - 1]
        set succ [expr $res + 1]
        set prev_c_sel [atomselect $mol "resid $prev and name C"]
        set prev_c [lindex [$prev_c_sel get index] 0]
        $prev_c_sel delete
        set res_n_sel [atomselect $mol "resid $res and name N"]
        set res_n [lindex [$res_n_sel get index] 0]
        $res_n_sel delete
        set res_ca_sel [atomselect $mol "resid $res and name CA"]
        set res_ca [lindex [$res_ca_sel get index] 0]
	set rname [lindex [$res_ca_sel get resname] 0]
        $res_ca_sel delete
        set res_c_sel [atomselect $mol "resid $res and name C"]
        set res_c [lindex [$res_c_sel get index] 0]
        $res_c_sel delete
	set succ_n_sel [atomselect $mol "resid $succ and name N"]
        set succ_n [lindex [$succ_n_sel get index] 0]
        $succ_n_sel delete
	puts "$res $rname $prev_c $res_n $res_ca $res_c $succ_n"
	set phi [measure dihed [list $prev_c $res_n $res_ca $res_c] molid $mol]
	set psi [measure dihed [list $res_n $res_ca $res_c $succ_n] molid $mol]
	if { ( $rname == "GLY" ) || ( $rname == "ALA")  } {
		return [list $phi $psi "x" "x"]
	} else {
		puts "marker 1"
		set cdname ""
		if { ( $rname == "SER" ) } {
			set cgname "OG"
		} elseif { ( $rname == "THR" ) } {
			set cgname "OG1"
		} elseif { ( $rname == "VAL" ) || ( $rname == "ILE" ) } {
			set cgname "CG1"
		} else {
			set cgname "CG"
		}
		puts "marker 2"
		set res_cb_sel [atomselect $mol "resid $res and name CB"]
		set n [$res_cb_sel num]
		set res_cb [lindex [$res_cb_sel get index] 0]
		$res_cb_sel delete
		set res_cg_sel [atomselect $mol "resid $res and name $cgname"]
                set res_cg [lindex [$res_cg_sel get index] 0]
                $res_cg_sel delete
		puts "$res $rname $res_ca $res_cb $res_cg"
		if { ($res_cb == "") || ($res_cg == "") } {
			return [list $phi $psi "x" "x"]
		}
		set chi1 [measure dihed [list $res_n $res_ca $res_cb $res_cg] molid $mol]
		if { ( $cgname == "CG" ) } {
			puts "marker 3"
			if { ( $rname == "MET" ) } {
				set cdname "SD"
			} elseif { ( $rname == "ASN" ) || ( $rname == "ASP" ) } {
				set cdname "OD1" 
			} elseif { ( $rname == "HID" ) || ( $rname == "HIE" ) || ( $rname == "HIP" ) || ( $rname == "HIS" ) } {
				set cdname "ND1" 
			} elseif { ( $rname == "LEU" ) || ( $rname == "ILE" ) || \
				( $rname == "PHE" ) || ( $rname == "TYR" ) || ( $rname == "TRP" ) } {
				set cdname "CD1"
			} else {
				set cdname "CD"
			}
			puts "marker 4"
			set res_cd_sel [atomselect $mol "resid $res and name $cdname"] 
			set res_cd [lindex [$res_cd_sel get index] 0]
			$res_cd_sel delete
			puts "$res $rname $res_ca $res_cb $res_cg $res_cd"
			set chi2 [measure dihed [list $res_ca $res_cb $res_cg $res_cd] molid $mol]
			#set chi2 "x"
		} else {
			set chi2 "x"
		}
		return [list $phi $psi $chi1 $chi2]
	} 
}



#proc dihrmsd {dih1 dih2} {
#	set n [llength $dih1]
	#set n2 [llength $dih2]
	#if
	


set aaregion_str $::env(aaregion_str)
set psf [lindex $argv 0]
set dcd [lindex $argv 1]
#set drugname [lindex $argv 3]
set ref [lindex $argv 2]
set outfname [lindex $argv 3]
#set drugu [string toupper $drugname]


mol new $ref
mol new $psf
#mol addfile $initpdb
mol addfile $dcd waitfor all
set nframes [molinfo top get numframes]
set nframesm1 [expr $nframes - 1]

set sel1 [atomselect top "$aaregion_str"]
set reslist [$sel1 get resid]
set reslist [lsort -integer -unique $reslist]
#puts $reslist
set output [open $outfname "w"]
animate goto $nframesm1
foreach res $reslist {
	set sel [atomselect top "resid $res and name CA"]
	set rname [lindex [$sel get resname] 0]
	$sel delete
	if { $res != 243 } {
		set phipsi1 [getphipsi 0 $res]
		set phipsi2 [getphipsi 1 $res]
		puts "$res $rname $phipsi1 $phipsi2"
		puts "marker 5"
		puts $output "$res $rname $phipsi1 $phipsi2"
	}
}
close $output
quit

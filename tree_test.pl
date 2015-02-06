# by Viacheslav Shalisko 2015

use strict;
use Bio::TreeIO;
use Bio::Tree::Draw::Cladogram;

my $output_file = "prueba.csv";
my $logfile = "pruebas_log.log";
# ArcGIS python path
my $python_path = "C:\\Python27\\ArcGIS10.1\\python.exe";


# redirecting warnings to log file
open ( LOG, ">$logfile" ) || die "$0: Can't open log file!";
open ( STDERR, ">>$logfile" );


my $treeio = Bio::TreeIO->new(-format => 'nexus', -file => 'RAxML_bipartitions.result_nexus');
my $binary_tree = $treeio->next_tree;
my $original_size = $binary_tree->number_nodes;

# the analysis works with binary trees only, so forcing binary tree with arbitrary balanced branch splitting
$binary_tree->force_binary();

my $binary_size = $binary_tree->number_nodes;
my $root = $binary_tree->get_root_node;

# set tag ID with internal code nodes for all nodes
for my $n ( $binary_tree->get_nodes ) {
	if( !defined $n->id() ) {
		my $node_name = "N" . $n->internal_id();
		$n->id($node_name);
	}
	#$n->set_tag_value("myID",$n->internal_id());
}

# save copy of new tree (after binary)
my $treeout = Bio::TreeIO->new(-format => 'nexus', -file => '>modified_tree.nex');
$treeout->write_tree($binary_tree);

print "Modified tree:\n" . $binary_tree->as_text; 

# hash of relation between terminal nodes (leafs) and vector layers with area distributions
# [todo: better to move vector layer names directly to NEXUS format and get them from input file as characters]
 my %geographic_hash = (
 					"Moranopteris_inaccessans_VT349" => "Moranopteris_inaccessa_u04_b1km_ver2",
 					"Moranopteris_longisetosa" => "Moranopteris_longisetosa_u018_b1km",
 					"Moranopteris_basiattenuata" => "Moranopteris_bieattenuata_u005_b1km_f1",
					"Moranopteris_aphelolepis_557" => "Moranopteris_aphelolepis_u05_b1km",
					"Moranopteris_blepharidea" => "Moranopteris_blepharidea_u05_b1km",
					"Moranopteris_hyalina_1148" => "Moranopteris_hyalina_u015_b1km",
					"Moranopteris_hyalina_1426" => "Moranopteris_hyalina_u015_b1km",
					"Moranopteris_gradata" => "Moranopteris_gradata_u015_b1km",
					"Moranopteris_setosa" => "Moranopteris_setosa_u01_b1km",
					"Moranopteris_perpusilla" => "Moranopteris_perpusilla_u05_b1km",
					"Moranopteris_achilleifolia_" => "Moranopteris_achilleifolia_u03_b1km",
					"Moranopteris_serricula_8084" => "Moranopteris_serricula_u0025_b1km_f3",
					"Moranopteris_zurquina_" => "Moranopteris_zurquina_u005_b1km_f3",
					"Moranopteris_trichomanoides_114" => "Moranopteris_trichomanoides_u01_b1km",
					"Moranopteris_sherringii_VT22" => "Moranopteris_sherringii_u01_b1km_f2",
					"Moranopteris_williamsii" => "Moranopteris_williamsii_u05_b1km",
					"Moranopteris_cookii" => "Moranopteris_cookei_u01_b1km",
					"Moranopteris_caucana_182" => "Moranopteris_caucana_u01_b1km",
					"Moranopteris_caucana_13445" => "Moranopteris_caucana_u01_b1km",
					"Moranopteris_truncicola_862" => "Moranopteris_truncicola_u02_b1km",
					"Moranopteris_truncicola_626" => "Moranopteris_truncicola_u02_b1km",
					"Moranopteris_taenifolia_08116" => "Moranopteris_taenifolia_u015_b1km_f3",
					"Moranopteris_taenifolia_sn" => "Moranopteris_taenifolia_u015_b1km_f3",
					"Moranopteris_microlepis_2584" => "Moranopteris_microlepis_u04_b1km",
					"Moranopteris_microlepis_2003" => "Moranopteris_microlepis_u04_b1km",
					"Moranopteris_plicata" => "Moranopteris_plicata_u025_b1km",
					"Moranopteris_nana" => "Moranopteris_nana_u02_b1km",
					"Moranopteris_grisebachii" => "Moranopteris_grisebachii_u01_b1km_f3",
 						);

my %area_data_hash = ();

open( DATAFILE, ">$output_file" ) || die "\nUnable to open file $output_file\n";
print DATAFILE "node_id,root_distance,null_intersects,intersect_score,number_of_leaves,combinations,incomplete_combinations\n";

    # process all the children in tree
    foreach my $node ( $root->get_all_Descendents() ) {
        if( $node->is_Leaf ) { 
    		# print "leaf node ID: " . $node->id() . "\n";
        } else {
    		print "Node: " . $node->internal_id() . "\n";

    		my $node_intersections_rate = 0;
    		my $combinations_incomplete = 0;
    		my $combinations = 0;
    		my $null_intersects_counter = 0;
        	# as we deal with binary tree, each node can have only 2 direct descendents
        	my ($descendentA, $descendentB) = $node->each_Descendent;

        	# making array of all nodes in each subclade, including the first node and then filter only leaf nodes
        	my @descendentsA = $descendentA->get_all_Descendents;
        	push @descendentsA, $descendentA;
        	@descendentsA = &select_leafs_only(@descendentsA);
        	my @descendentsB = $descendentB->get_all_Descendents;
        	push @descendentsB, $descendentB;
        	@descendentsB = &select_leafs_only(@descendentsB);

        	my $node_leaf_counter = scalar @descendentsA + scalar @descendentsB;

        	foreach my $descendentA_node (@descendentsA) {
        		foreach my $descendentB_node (@descendentsB) {
        			$combinations++;
        			my $two_leaf_intersect = undef;
        			my $node_distance = &get_node_count_distance($descendentA_node, $descendentB_node, $binary_tree);

        			if( $geographic_hash{$descendentA_node->id()} && $geographic_hash{$descendentB_node->id()} ) {
        				my $area_A = $geographic_hash{$descendentA_node->id()};
        				my $area_B = $geographic_hash{$descendentB_node->id()};
        				

        				# check if we already have the intersect value in hash of geographic interscts
        				if( defined $area_data_hash{$area_A}{$area_B} ) {
        					# value is in hash
        					# do calculations using stored value
        					print "#";
        					$two_leaf_intersect = &calculate_area_part(
        								$area_data_hash{$area_A}{$area_B}->{A},
        								$area_data_hash{$area_A}{$area_B}->{B},
        								$area_data_hash{$area_A}{$area_B}->{Intersect},
        								);
        				} elsif( defined $area_data_hash{$area_B}{$area_A} ) {
        					print "#";
        					$two_leaf_intersect = &calculate_area_part(
        								$area_data_hash{$area_B}{$area_A}->{A},
        								$area_data_hash{$area_B}{$area_A}->{B},
        								$area_data_hash{$area_B}{$area_A}->{Intersect},
        								);        					
        				} else {
        					# value is not in hash, so we need to get it with external script
        					# output layer name include related node names and timestamp
        					my $area_output = "__N" . $node->internal_id() . "_A" . $descendentA_node->internal_id()
        									 . "_B" . $descendentB_node->internal_id() . "_T" . time;
        					# execute geographic intersect script
        					###################### external script intersect_areas.py call ##############################
        					my $external_output = `$python_path intersect_areas.py $area_A $area_B $area_output`;
        					if($external_output =~ /Exception/sm) {
        						# External script has generated some error so we don´t have result of this execution
        						print "!";
        						$combinations_incomplete++;
        					} else {
        						# External script has returned correct results
        						my ($out1, $out2, $out3) = split(/\s+/sm, $external_output);
         						$area_data_hash{$area_A}{$area_B} = {
        																A => $out1,
        																B => $out2,
        																Intersect => $out3,
        															    };
        						if($out3 > 0) {
        							print "+";
        						} else {
        							print "o";
        							$null_intersects_counter++;
        						}
        						$two_leaf_intersect = &calculate_area_part($out1, $out2, $out3);       				
        					}
        				}
        			} else {
        				# Skip leaf pair as we have no distribution data for one leaf or both
        				print ".";
        				$combinations_incomplete++;
        			}

        			my $two_leaf_intersect_component = 1 / ( 2 ** ( $node_distance - 1 ) ) * $two_leaf_intersect;
        			$node_intersections_rate = $node_intersections_rate + $two_leaf_intersect_component;
        		}
        	}

        	my @pair_node_root = ($node,$root);
	       	my $node_root_distance = $binary_tree->distance(-nodes => \@pair_node_root );

	       	if( $combinations_incomplete > 0) {
	       		print "\nD:". $node_root_distance . " I: incomplete\n";
	       		print DATAFILE $node->internal_id() ."," . 
	       				$node_root_distance . "," . 
	       				$null_intersects_counter .
	       				",incomplete," . 
	       				$node_leaf_counter . "," . 
	       				$combinations . "," . 
	       				$combinations_incomplete  . "\n";
	       	} else {
	       		print "\nD:". $node_root_distance . " I:" . $node_intersections_rate ."\n";
	       		print DATAFILE $node->internal_id() . "," . 
	       				$node_root_distance . "," . 
	       			    $null_intersects_counter . "," . 
	       				$node_intersections_rate . ",". 
	       				$node_leaf_counter . "," . 
						$combinations . "," . 
						$combinations_incomplete  . "\n";
	       	}
        	
        }
    }

close(DATAFILE);

print "Total node count in original tree: " . $original_size . "\n";
print "Total node count after convertion to binary tree: " . $binary_size . "\n";




sub get_node_count_distance {
	# get the number of nodes between any two nodes in tree

	my ($A,$B,$T) = @_; # node A, node B and tree objects of Bio::TreeIO
    my @A_linage = $T->get_lineage_nodes($A);
    my @B_linage = $T->get_lineage_nodes($B);

    my %node_count_hash = ();
    my @node_count_array = ();
    foreach my $N (@A_linage, @B_linage) { 
        $node_count_hash{$N}++; 
    }
    foreach my $node_counted (keys %node_count_hash) {
        if( $node_count_hash{$node_counted} == 1 ) {
        	push @node_count_array, $node_counted;
        }
    }
    # the node count array incude all unique nodes in the connection way between A & B, except their LCA, so we need to add 1
    return scalar @node_count_array + 1;
}

sub select_leafs_only {
	# return only leaf nodes form array of tree nodes

	my (@A) = @_;
	my @B = ();
	foreach my $N (@A) {
		if( $N->is_Leaf ) {
			push @B, $N;
		}
	}
	return @B;
}

sub calculate_area_part {

	my ($A,$B,$I) = @_;
	my $R = undef;
	if( $A != 0 and $A < $B ) {
		$R = $I / $A;
	} elsif( $B != 0 and $A >= $B ) {
		$R = $I / $B;
	} else {
		$R = 0;
	}
	return $R;

}

$SIG{__DIE__} = sub
{
    my @loc = caller(1);
    print STDOUT "Scripti execution interrumped, results are incomplete.\n";
    print STDOUT "Exception happened at line $loc[2] in $loc[1]:\n" . @_ . "\n";
    return 1;
};
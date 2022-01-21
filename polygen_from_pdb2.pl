use application "polytope";
use strict;
use warnings;

use File::Basename;
my $dirname = dirname(__FILE__);

my $name = "7N7F";
my $fullname = $name . ".pdb";
my $gz_fullname = $fullname . ".gz";

my $atomdata = load($dirname . "/" . $name . "-map.txt");

my $limit = 1000;
my @atom_coords_1;
my @atom_coords_2;
my @atom_coords_3;


my $atom_maps_1 = new Map<Int, Array<String>>();
my $atom_maps_2 = new Map<Int, Array<String>>();
my $atom_maps_3 = new Map<Int, Array<String>>();

sub seperate_amino {
  my $atom_seq = $_[0];
  my $origin = $_[1];
  
  my $new_atom_coord = new Vector<Rational>();
  my $amino_acid = "";
  my @atom_vertices;
  my @aminos = ();
  
  while (my ($k,$v)=each %$atom_seq) {
    if ($k == 0) {
      next;
    }
    
    if ($amino_acid ne $v->[1]) {
      $amino_acid = $v->[1];

      # finish with previous sequence
      if (@$new_atom_coord){
	print $new_atom_coord, "\n";
	push @atom_vertices, $new_atom_coord;

	my $atom_seq_pm = new Array<Vector<Rational>>(\@atom_vertices);

	print $atom_seq_pm, "\n";
	my $amino_hull = new Polytope<Rational>(POINTS=>$atom_seq_pm);
	my $translated_hull = translate($amino_hull, -$origin);
	#print $amino_hull;
	push @aminos, $translated_hull;
      }
      # start new amino sequence
      @atom_vertices = ();
      $new_atom_coord = new Vector<Rational>([1,$v->[4],$v->[5],$v->[6]]);
      push @atom_vertices, $new_atom_coord;
    } else {
      $new_atom_coord = new Vector<Rational>([1,$v->[4],$v->[5],$v->[6]]);
      push @atom_vertices, $new_atom_coord;
	
    }
  }

  return @aminos;
}

while (my ($k,$v)=each %$atomdata) {
   if ($k == 0){
      next;
    }

   my $chain = $v->[2];
   $chain =~ s/^\s+|\s+$//g;
   if ($chain eq "A"){
      my $new_atom_coord_1 = new Vector<Rational>([1,$v->[4],$v->[5],$v->[6]]);
      push @atom_coords_1,$new_atom_coord_1;
      $atom_maps_1->{$k} = $v;
   }
   if ($chain eq "B"){
      my $new_atom_coord_2 = new Vector<Rational>([1,$v->[4],$v->[5],$v->[6]]);
      push @atom_coords_2,$new_atom_coord_2;
      $atom_maps_2->{$k} = $v;
   }
   if ($chain eq "C"){
      my $new_atom_coord_3 = new Vector<Rational>([1,$v->[4],$v->[5],$v->[6]]);
      push @atom_coords_3,$new_atom_coord_3;
      $atom_maps_3->{$k} = $v;
   }

}

  
#Centering

my @atom_coords;
push @atom_coords,@atom_coords_1;
push @atom_coords,@atom_coords_2;
push @atom_coords,@atom_coords_3;
my $atom_coords_pm = new Array<Vector<Rational>>(\@atom_coords);
my $P = new Polytope<Rational>(POINTS=>$atom_coords_pm);
my $bary = new Vector<Rational>( @{barycenter($P->VERTICES)}[1..3]);


my $atom_coords_1_pm = new Array<Vector<Rational>>(\@atom_coords_1);
my $P1 = new Polytope<Rational>(POINTS=>$atom_coords_1_pm);
$P1 = translate($P1,-$bary);
#my $atom_coords_2_pm = new Array<Vector<Rational>>(\@atom_coords_2);
#my $P2 = new Polytope<Rational>(POINTS=>$atom_coords_2_pm);
#$P2 = translate($P2,-$bary);
#my $atom_coords_3_pm = new Array<Vector<Rational>>(\@atom_coords_3);
#my $P3 = new Polytope<Rational>(POINTS=>$atom_coords_3_pm);
#$P3 = translate($P3,-$bary);
#compose($P1->VISUAL(FacetColor=>"red"),$P2->VISUAL(FacetColor=>"blue"),$P3->VISUAL(FacetColor=>"green"));


my @chain_aminos = seperate_amino($atom_maps_1, $bary);

compose(map($_->VISUAL, @chain_aminos));


#save($P,$dirname . "/" . $name . "-nonscaled-" . $limit . ".poly");

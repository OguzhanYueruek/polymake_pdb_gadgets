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

while (my ($k,$v)=each %$atomdata) {
   if ($k == 0){
      next;
   }
   my $chain = $v->[2];
   $chain =~ s/^\s+|\s+$//g;
   if ($chain eq "A"){
      print "foo\n";
      my $new_atom_coord_1 = new Vector<Rational>([1,$v->[4],$v->[5],$v->[6]]);
      push @atom_coords_1,$new_atom_coord_1;

   }
   if ($chain eq "B"){
      my $new_atom_coord_2 = new Vector<Rational>([1,$v->[4],$v->[5],$v->[6]]);
      push @atom_coords_2,$new_atom_coord_2;

   }
   if ($chain eq "C"){
      my $new_atom_coord_3 = new Vector<Rational>([1,$v->[4],$v->[5],$v->[6]]);
      push @atom_coords_3,$new_atom_coord_3;

   }

}

#Centering
=pod
my @atom_coords = ();
push @atom_coords,@atom_coords_1;
push
my $atom_coords_pm = new Array<Vector<Rational>>(\@atom_coords);
my $P = new Polytope<Rational>(POINTS=>$atom_coords_pm);
my $bary = new Vector<Rational>( @{barycenter($P->VERTICES)}[1..3]);
=cut

my $atom_coords_1_pm = new Array<Vector<Rational>>(\@atom_coords_1);
my $P1 = new Polytope<Rational>(POINTS=>$atom_coords_1_pm);
$P1 = translate($P1,-$bary);
my $atom_coords_2_pm = new Array<Vector<Rational>>(\@atom_coords_2);
my $P2 = new Polytope<Rational>(POINTS=>$atom_coords_2_pm);
$P2 = translate($P2,-$bary);
my $atom_coords_3_pm = new Array<Vector<Rational>>(\@atom_coords_3);
my $P3 = new Polytope<Rational>(POINTS=>$atom_coords_3_pm);
$P3 = translate($P3,-$bary);
compose($P1->VISUAL(FacetColor=>"red"),$P2->VISUAL(FacetColor=>"blue"),$P3->VISUAL(FacetColor=>"green"));

#save($P,$dirname . "/" . $name . "-nonscaled-" . $limit . ".poly");

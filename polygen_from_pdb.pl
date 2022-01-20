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
my @atom_coords;
while (my ($k,$v)=each %$atomdata) {
   if ($k == 0){
      next;
   }
   #print $k, ":      ", $v, "\n";
   my $new_atom_coord = new Vector<Rational>([1,$v->[4],$v->[5],$v->[6]]);
   push @atom_coords,$new_atom_coord;

   if($k>$limit){
      last;
   }
}
print @atom_coords;
my $atom_coords_pm = new Array<Vector<Rational>>(\@atom_coords);
my $P = new Polytope<Rational>(POINTS=>$atom_coords_pm);
my $bary = new Vector<Rational>( @{barycenter($P->VERTICES)}[1..3]);
$P = translate($P,-$bary);
 $P->VISUAL;
print "\n",$P->VOLUME;
save($P,$dirname . "/" . $name . "-nonscaled-" . $limit . ".poly");

use application "polytope";
use strict;
use warnings;

use File::Basename;
my $dirname = dirname(__FILE__);

my $name = "7N7F";
my $fullname = $name . ".pdb";
my $gz_fullname = $fullname . ".gz";

my $atomdata = load($dirname . "/" . $name . "-map.txt");

my $limit = 0;
my @atom_coords;

my $total_num_aminos = 1;
my $num_current_aminos = 0;
my $amino_acid = "";
my @geodesic_vertices;
my $bounding_geodesic_hull = new Polytope<Rational>();
my $atom_num = 0;

while (my ($k,$v)=each %$atomdata) {
  if ($k == 0) {
    next;
  }
  
  my $new_atom_coord = new Vector<Rational>([1,$v->[4],$v->[5],$v->[6]]);
  

  if ($amino_acid ne $v->[1]) {
    $amino_acid = $v->[1];
    $num_current_aminos += 1;

  
    if ($num_current_aminos > $total_num_aminos) {

      push @geodesic_vertices, $new_atom_coord;
      print "END ", $new_atom_coord, " * \n";

      my $geodesic_coords_pm = new Array<Vector<Rational>>(\@geodesic_vertices);

      $bounding_geodesic_hull = new Polytope<Rational>(POINTS=>$geodesic_coords_pm);
      last;
      
    } else {
      print $amino_acid, " AMINO ACID NUM ", $num_current_aminos, " *** \n";
      print "START ", $new_atom_coord, " * \n";
      push @geodesic_vertices, $new_atom_coord;
      print @geodesic_vertices, "\n";

    }
  }


# take approximate mid point
  if ($atom_num == 3 || $atom_num == 8) {
    $new_atom_coord = new Vector<Rational>([1,$v->[4],$v->[5],$v->[6]]);
    print "MID ", $new_atom_coord, " * \n";
    push @geodesic_vertices, $new_atom_coord;
  }
  
  push @atom_coords,$new_atom_coord;

  $atom_num += 1;
  
  #print $v->[0], " ", $new_atom_coord, " ** \n ";

   
   #print $v,  "\n";
   #if($k>$limit){
   #   last;
   #}
}
#print @atom_coords;
my $atom_coords_pm = new Array<Vector<Rational>>(\@atom_coords);
my $P = new Polytope<Rational>(POINTS=>$atom_coords_pm);
my $bary = new Vector<Rational>( @{barycenter($P->VERTICES)}[1..3]);
$P = translate($P,-$bary);
$bounding_geodesic_hull = translate($bounding_geodesic_hull ,-$bary);
compose($P->VISUAL,$bounding_geodesic_hull->VISUAL(FacetColor=>'blue'));


#print "\n",$P->VOLUME;
#2save($P,$dirname . "/" . $name . "-nonscaled-" . $limit . ".poly");

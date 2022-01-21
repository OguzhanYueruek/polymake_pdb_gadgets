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

sub seperate_amino {
  my $atom_seq = $_[0];
  my $new_atom_coord = new Vector<Rational>();
  my $amino_acid = "";
  my @atom_vertices;
  my @aminos = [];
  
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
	#print $amino_hull;
	push @aminos, $amino_hull;
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

my $amino_hulls = seperate_amino($atomdata);
print $amino_hulls;


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

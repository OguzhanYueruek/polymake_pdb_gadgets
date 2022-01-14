use application "polytope";
use strict;
use warnings;

use File::Basename;
use IO::Zlib;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);


my $dirname = dirname(__FILE__);

my $name = "7N7F";
my $fullname = $name . ".pdb";
my $gz_fullname = $fullname . ".gz";

my $url = "https://files.rcsb.org/download/".$fullname;
print "Name is: $name\n";
print "Using url: $url\n";

my $path_to_gz = $dirname . "/" . $gz_fullname;
`curl $url > $path_to_gz`;
print "Stored as $path_to_gz\n";

my $path_to_pdb = $dirname . "/" . $fullname;
gunzip $path_to_gz => $path_to_pdb
    or die "gunzip failed: $GunzipError\n";
print "Extracted $path_to_pdb\n";






chomp $path_to_pdb;

unless (open(INPUTFILE, $path_to_pdb)) {
    print "Cannot read from '$path_to_pdb'.\nProgram closing.\n";
    <STDIN>;
    exit;
}


# load the file into an array
chomp(my @data = <INPUTFILE>);

my $total = scalar @data;
print $total;

# close the file
close(INPUTFILE);


my $atomdata = new Array<Array<String>>($total);
my $count = 0;

my @propnames = ("Atom No","AtomPos","Aminoacid","Chain","Sequence No","X-axis","Y-axis","Z-axis","Element");
$atomdata->[0] = new Array<String>(@propnames);
$count++;

for (my $line = 0; $line <scalar @data; $line++) {
    if ($data[$line]=~/^HEADER\s+(.*?)$/) {
       my $header = $1;
       print $header,"\n";
       print "$data[$line]\n";
    }
    if ($data[$line]=~/^TITLE\s+(.*?)$/) {
        my %parsing;
        $parsing{$line} = $1;
        print "$parsing{$line}\n";
    }

    if ($data[$line] =~ m/ATOM\s+(\d+)\s+(\S+)\s+(\w{3})\s+(.+)\s+(\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+).+(\w{1})/i) {

        my @linedata = ($1,$2,$3,$4,$5,$6,$7,$8,$9);
        $atomdata->[$count] = new Array<String>(@linedata);
        $count++;
        #$atom_data{"Atom No"} = $1;
    }


}
my @atomdata=@{$atomdata};
@atomdata = @atomdata[0..$count];
$atomdata = new Array<String>(@atomdata);

save($atomdata,$dirname . "/" . $name . "-polymake.txt");
for (my $i=0 ; $i<$count ; $i++) {
    print $atomdata->[$i],"\n";
}

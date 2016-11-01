#!/usr/bin/perl
# Lattice plotting utility
#
# the number of an atom to start with in each layer
my $afrom = -6;

# the number an atom to stop in each layer
my $ato = 35;

# the number of layers to expand
my $nlay = 12;

while(<>) {
  if(/\[(.*)\]/) {
     my @a = split(/,/, $1);
     print "# ", join(':', @a), "\n";    
     my $s = @a;
     print "# array of size ", $s, "\n";               
     my $nl = $s / 3;
     my $j;
     print "# plot ";
     for($j = 1; $j <= $nlay; $j ++) {       
       print "\"\" u ", 2 * $j - 1, ":", 2 * $j, " pt 7 ps 2 lc -1";
       if($j == $nlay) {
         print "\n\n\n";
       } else {
         print ", ";
       }
     }
     
     for($j = $afrom; $j <= $ato; $j ++) { 
       my $i;
       my $y;
       for($i = 0; $i < $nlay; $i ++) {
         my $patl = $a[0] + $a[3*($nl - 1)];
         my $basel = $i % $nl;
         my $npats = ($i - $basel) / $nl;
         if($i == 0) {
           $y = 0;
         } else {
           $y += $a[3 * $basel]
         }
         my $d = $a[3 * $basel + 1];
         my $s = $a[3 * $basel + 2];
         my $x = $d + $j * $s;
         print $x, " ", $y, " ";
       }
       print "\n";
     }
  }
}


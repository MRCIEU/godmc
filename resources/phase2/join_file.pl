#!/software/bin/perl-5.14.2 -w

use Getopt::Std;

getopts("i:o:a:",\%opts);
(my $input1, my $input2) = split(/\s+/, $opts{i});
my $out_file = $opts{o};
my $all = $opts{a}? $opts{a}: 0;
(my $file1, my $separator1, my $col1) = split(/,/, $input1);
(my $file2, my $separator2, my $col2) = split(/,/, $input2);
my $sep="";

my %id_HASH = ();

if ($file2 =~ /gz$/) {
	open FILE2, "zcat $file2 |" or die $!;
} else {
	open FILE2, "< $file2" or die $!;
}
while (<FILE2>) {
	next if /^$/;
        chomp;
        my @temp = ();
        if ($separator2 eq "COMMA") {
                @temp = split(/,/);
		$sep=",";
        } elsif ($separator2 eq "TAB") {
                @temp = split(/\t/);
		$sep="\t";
        } elsif ($separator2 eq "SPACE") {
		@temp = split(/\s+/);
		$sep=" ";
	}
        $id_HASH{$temp[$col2]} = $_;
}
close FILE2;

my $file2_cols =0;
if ($file2 =~ /gz$/) {
	$file2_cols = `zcat $file2 | awk 'NR==1 {printf NF}'`;
} else {
	$file2_cols = `cat $file2 | awk 'NR==1 {printf NF}'`;
}

if ($file1 =~ /gz$/) {
	open FILE1, "zcat $file1 |" or die $!;
} else {
	open FILE1, "< $file1" or die $!;
}
open OUT, "> $out_file" or die $!;
while (<FILE1>) {
	next if /^$/;
        chomp;
        my @temp = ();
        if ($separator1 eq "COMMA") {
                @temp = split(/,/);
	} elsif ($separator2 eq "TAB") {
                @temp = split(/\t/);
        } elsif ($separator2 eq "SPACE") {
                @temp = split(/\s+/);
        }
        if (exists $id_HASH{$temp[$col1]}) {
                print OUT $_ . $sep. $id_HASH{$temp[$col1]} . "\n";
        } elsif ($all==1) {
                print OUT "$_" . $sep; 
		for (my $i=1; $i<=$file2_cols; $i++) {
			print OUT $sep . "NA";
		}
		print OUT "\n";
        }
}
close FILE1;
close OUT;



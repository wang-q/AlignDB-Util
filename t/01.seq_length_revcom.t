#!/usr/bin/perl
use strict;
use warnings;

use Test::More qw(no_plan);    #tests => 33;

BEGIN {
    use_ok( 'AlignDB::Util', qw(seq_length revcom) );
}

{
    print "#seq_length\n";

    my @seqs = (
        [qw{ AAAA 4 }], [qw{ CCCC 4 }],
        [qw{ TAGGGATAACAGGGTAAT 18 }],
        [qw{ GCAN--NN--NNTGC 11 }],
    );

    foreach my $i ( 0 .. @seqs - 1 ) {
        my ( $ori, $expected ) = @{ $seqs[$i] };
        my $result1 = seq_length($ori);
        my $result2 = seq_length( \$ori );
        print "\n";
        print "original: $ori\n";
        is( $result1, $expected, "seq_length1_$i" );
        is( $result2, $expected, "seq_length2_$i" );
    }
}

{
    print "#revcom\n";

    my @seqs = (
        [qw{ AAAA TTTT }],
        [qw{ CCCC GGGG }],
        [qw{ TAGGGATAACAGGGTAAT ATTACCCTGTTATCCCTA }],    # I-Sce I endonuclease
        [qw{ GCANNNNNTGC GCANNNNNTGC }],                  # BstAP I
    );

    foreach my $i ( 0 .. @seqs - 1 ) {
        my ( $ori, $expected ) = @{ $seqs[$i] };
        my $result1 = revcom($ori);
        my $result2 = revcom( \$ori );
        print "\n";
        print "original: $ori\n";
        is( $result1, $expected, "revcom1_$i" );
        is( $result2, $expected, "revcom2_$i" );
    }
}

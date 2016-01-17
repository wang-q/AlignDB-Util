#!/usr/bin/perl
use strict;
use warnings;

use Test::More tests => 33;

BEGIN {
    use_ok( 'AlignDB::Util', qw(ref_indel_type) );
}

{
    print "#indel type\n";

    my @indel1  = ( [qw{ -- -- AT }],       "Q", "I" );
    my @indel2  = ( [qw{ -- AT --}],        "T", "I" );
    my @indel3  = ( [qw{ AT -- AT }],       "T", "D" );
    my @indel4  = ( [qw{ AT AT -- }],       "Q", "D" );
    my @indel5  = ( [qw{ AAA --- A-- }],    "C", "C" );
    my @indel6  = ( [qw{ AAA --- --- }],    "N", "N" );
    my @indel7  = ( [qw{ AAA A-- A-- }],    "N", "N" );
    my @indel8  = ( [qw{ AA- --- AA- }],    "T", "D" );
    my @indel9  = ( [qw{ --- T-- T-- }],    "N", "N" );
    my @indel10 = ( [qw{ T--- T--A ---A }], "C", "C" );
    my @indel11 = ( [qw{ ---- -AT- ----}],  "T", "I" );
    my @indel12 = ( [qw{ A--T AATT A--T}],  "T", "I" );
    my @indel13 = ( [qw{ -AT- -AT- ---- }], "Q", "D" );
    my @indel14 = ( [qw{ AATT AATT A--T }], "Q", "D" );
    my @indel15 = ( [qw{ - - - }],          "N", "N" );
    my @indel16 = ( [qw{ ---- ---- ---- }], "N", "N" );

    my @indels = (
        \@indel1, \@indel2,  \@indel3,  \@indel4,  \@indel5,  \@indel6,  \@indel7,  \@indel8,
        \@indel9, \@indel10, \@indel11, \@indel12, \@indel13, \@indel14, \@indel15, \@indel16,
    );

    foreach my $i ( 0 .. @indels - 1 ) {
        my ( $indel, $expect_occured, $expect_type ) = @{ $indels[$i] };
        my ( $indel_occured, $indel_type ) = &ref_indel_type(@$indel);
        print "\n";
        print "occured: $indel_occured\n";
        print "type:    $indel_type\n";
        is( $indel_occured, $expect_occured, "occured_$i" );
        is( $indel_type,    $expect_type,    "type_$i" );
    }
}

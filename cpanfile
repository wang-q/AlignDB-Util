requires 'Moose';
requires 'Carp';
requires 'YAML';
requires 'Path::Tiny';
requires 'List::Util';
requires 'List::MoreUtils';
requires 'Math::Combinatorics';
requires 'Statistics::Descriptive';
requires 'Bio::Seq';
requires 'Bio::Tools::Run::Alignment::Clustalw';
requires 'AlignDB::IntSpan';
requires 'perl', '5.008001';

on test => sub {
    requires 'Test::More', 0.88;
    requires 'Test::Number::Delta';
};

package STAR::STARClient;

use JSON::RPC::Client;
use POSIX;
use strict;
use Data::Dumper;
use URI;
use Bio::KBase::Exceptions;
my $get_time = sub { time, 0 };
eval {
    require Time::HiRes;
    $get_time = sub { Time::HiRes::gettimeofday() };
};

use Bio::KBase::AuthToken;

# Client version should match Impl version
# This is a Semantic Version number,
# http://semver.org
our $VERSION = "0.1.0";

=head1 NAME

STAR::STARClient

=head1 DESCRIPTION


Name of module: STAR

This KBase module wraps the free open source software STAR: ultrafast universal RNA-seq aligner.
STAR-2.6.1a

References:
https://github.com/alexdobin/STAR/
https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf


=cut

sub new
{
    my($class, $url, @args) = @_;
    

    my $self = {
	client => STAR::STARClient::RpcClient->new,
	url => $url,
	headers => [],
    };

    chomp($self->{hostname} = `hostname`);
    $self->{hostname} ||= 'unknown-host';

    #
    # Set up for propagating KBRPC_TAG and KBRPC_METADATA environment variables through
    # to invoked services. If these values are not set, we create a new tag
    # and a metadata field with basic information about the invoking script.
    #
    if ($ENV{KBRPC_TAG})
    {
	$self->{kbrpc_tag} = $ENV{KBRPC_TAG};
    }
    else
    {
	my ($t, $us) = &$get_time();
	$us = sprintf("%06d", $us);
	my $ts = strftime("%Y-%m-%dT%H:%M:%S.${us}Z", gmtime $t);
	$self->{kbrpc_tag} = "C:$0:$self->{hostname}:$$:$ts";
    }
    push(@{$self->{headers}}, 'Kbrpc-Tag', $self->{kbrpc_tag});

    if ($ENV{KBRPC_METADATA})
    {
	$self->{kbrpc_metadata} = $ENV{KBRPC_METADATA};
	push(@{$self->{headers}}, 'Kbrpc-Metadata', $self->{kbrpc_metadata});
    }

    if ($ENV{KBRPC_ERROR_DEST})
    {
	$self->{kbrpc_error_dest} = $ENV{KBRPC_ERROR_DEST};
	push(@{$self->{headers}}, 'Kbrpc-Errordest', $self->{kbrpc_error_dest});
    }

    #
    # This module requires authentication.
    #
    # We create an auth token, passing through the arguments that we were (hopefully) given.

    {
	my %arg_hash2 = @args;
	if (exists $arg_hash2{"token"}) {
	    $self->{token} = $arg_hash2{"token"};
	} elsif (exists $arg_hash2{"user_id"}) {
	    my $token = Bio::KBase::AuthToken->new(@args);
	    if (!$token->error_message) {
	        $self->{token} = $token->token;
	    }
	}
	
	if (exists $self->{token})
	{
	    $self->{client}->{token} = $self->{token};
	}
    }

    my $ua = $self->{client}->ua;	 
    my $timeout = $ENV{CDMI_TIMEOUT} || (30 * 60);	 
    $ua->timeout($timeout);
    bless $self, $class;
    #    $self->_validate_version();
    return $self;
}




=head2 run_star

  $returnVal = $obj->run_star($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a STAR.AlignReadsParams
$returnVal is a STAR.AlignReadsResult
AlignReadsParams is a reference to a hash where the following keys are defined:
	readsset_ref has a value which is a STAR.obj_ref
	genome_ref has a value which is a STAR.obj_ref
	output_workspace has a value which is a string
	output_name has a value which is a string
	alignment_suffix has a value which is a string
	condition has a value which is a string
	concurrent_njsw_tasks has a value which is an int
	concurrent_local_tasks has a value which is an int
	outSAMunmapped has a value which is a string
	create_report has a value which is a STAR.bool
	alignmentset_suffix has a value which is a string
	alignIntronMin has a value which is an int
	alignIntronMax has a value which is an int
	alignMatesGapMax has a value which is an int
	alignSJoverhangMin has a value which is an int
	alignSJDBoverhangMin has a value which is an int
	quantMode has a value which is a string
	outFilterType has a value which is a string
	outFilterMultimapNmax has a value which is an int
	outSAMtype has a value which is a string
	outSAMattrIHstart has a value which is an int
	outSAMstrandField has a value which is a string
	outFilterMismatchNmax has a value which is an int
	outFileNamePrefix has a value which is a string
	runThreadN has a value which is an int
obj_ref is a string
bool is an int
AlignReadsResult is a reference to a hash where the following keys are defined:
	output_directory has a value which is a string
	report_name has a value which is a string
	report_ref has a value which is a STAR.obj_ref
	alignmentset_ref has a value which is a STAR.obj_ref
	alignment_objs has a value which is a reference to a hash where the key is a STAR.obj_ref and the value is a STAR.AlignmentObj
AlignmentObj is a reference to a hash where the following keys are defined:
	ref has a value which is a STAR.obj_ref
	name has a value which is a string

</pre>

=end html

=begin text

$params is a STAR.AlignReadsParams
$returnVal is a STAR.AlignReadsResult
AlignReadsParams is a reference to a hash where the following keys are defined:
	readsset_ref has a value which is a STAR.obj_ref
	genome_ref has a value which is a STAR.obj_ref
	output_workspace has a value which is a string
	output_name has a value which is a string
	alignment_suffix has a value which is a string
	condition has a value which is a string
	concurrent_njsw_tasks has a value which is an int
	concurrent_local_tasks has a value which is an int
	outSAMunmapped has a value which is a string
	create_report has a value which is a STAR.bool
	alignmentset_suffix has a value which is a string
	alignIntronMin has a value which is an int
	alignIntronMax has a value which is an int
	alignMatesGapMax has a value which is an int
	alignSJoverhangMin has a value which is an int
	alignSJDBoverhangMin has a value which is an int
	quantMode has a value which is a string
	outFilterType has a value which is a string
	outFilterMultimapNmax has a value which is an int
	outSAMtype has a value which is a string
	outSAMattrIHstart has a value which is an int
	outSAMstrandField has a value which is a string
	outFilterMismatchNmax has a value which is an int
	outFileNamePrefix has a value which is a string
	runThreadN has a value which is an int
obj_ref is a string
bool is an int
AlignReadsResult is a reference to a hash where the following keys are defined:
	output_directory has a value which is a string
	report_name has a value which is a string
	report_ref has a value which is a STAR.obj_ref
	alignmentset_ref has a value which is a STAR.obj_ref
	alignment_objs has a value which is a reference to a hash where the key is a STAR.obj_ref and the value is a STAR.AlignmentObj
AlignmentObj is a reference to a hash where the following keys are defined:
	ref has a value which is a STAR.obj_ref
	name has a value which is a string


=end text

=item Description



=back

=cut

 sub run_star
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function run_star (received $n, expecting 1)");
    }
    {
	my($params) = @args;

	my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to run_star:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'run_star');
	}
    }

    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
	    method => "STAR.run_star",
	    params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'run_star',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method run_star",
					    status_line => $self->{client}->status_line,
					    method_name => 'run_star',
				       );
    }
}
 
  
sub status
{
    my($self, @args) = @_;
    if ((my $n = @args) != 0) {
        Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
                                   "Invalid argument count for function status (received $n, expecting 0)");
    }
    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
        method => "STAR.status",
        params => \@args,
    });
    if ($result) {
        if ($result->is_error) {
            Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
                           code => $result->content->{error}->{code},
                           method_name => 'status',
                           data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
                          );
        } else {
            return wantarray ? @{$result->result} : $result->result->[0];
        }
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method status",
                        status_line => $self->{client}->status_line,
                        method_name => 'status',
                       );
    }
}
   

sub version {
    my ($self) = @_;
    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
        method => "STAR.version",
        params => [],
    });
    if ($result) {
        if ($result->is_error) {
            Bio::KBase::Exceptions::JSONRPC->throw(
                error => $result->error_message,
                code => $result->content->{code},
                method_name => 'run_star',
            );
        } else {
            return wantarray ? @{$result->result} : $result->result->[0];
        }
    } else {
        Bio::KBase::Exceptions::HTTP->throw(
            error => "Error invoking method run_star",
            status_line => $self->{client}->status_line,
            method_name => 'run_star',
        );
    }
}

sub _validate_version {
    my ($self) = @_;
    my $svr_version = $self->version();
    my $client_version = $VERSION;
    my ($cMajor, $cMinor) = split(/\./, $client_version);
    my ($sMajor, $sMinor) = split(/\./, $svr_version);
    if ($sMajor != $cMajor) {
        Bio::KBase::Exceptions::ClientServerIncompatible->throw(
            error => "Major version numbers differ.",
            server_version => $svr_version,
            client_version => $client_version
        );
    }
    if ($sMinor < $cMinor) {
        Bio::KBase::Exceptions::ClientServerIncompatible->throw(
            error => "Client minor version greater than Server minor version.",
            server_version => $svr_version,
            client_version => $client_version
        );
    }
    if ($sMinor > $cMinor) {
        warn "New client version available for STAR::STARClient\n";
    }
    if ($sMajor == 0) {
        warn "STAR::STARClient version is $svr_version. API subject to change.\n";
    }
}

=head1 TYPES



=head2 assembly_ref

=over 4



=item Description

A 'typedef' allows you to provide a more specific name for
a type.  Built-in primitive types include 'string', 'int',
'float'.  Here we define a type named assembly_ref to indicate
a string that should be set to a KBase ID reference to an
Assembly data object.


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 bool

=over 4



=item Description

A boolean - 0 for false, 1 for true.
@range (0, 1)


=item Definition

=begin html

<pre>
an int
</pre>

=end html

=begin text

an int

=end text

=back



=head2 obj_ref

=over 4



=item Description

An X/Y/Z style reference


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 AlignReadsParams

=over 4



=item Description

Will align the input reads (or set of reads specified in a SampleSet) to the specified
assembly or assembly for the specified Genome (accepts Assembly, ContigSet, or Genome types)
and produces a ReadsAlignment object, or in the case of a SampleSet, a ReadsAlignmentSet object

obj_ref genome_ref: KBase workspace reference Genome
obj_ref readsset_ref: the workspace reference for the set of reads to align, referring to 
                    either a SingleEnd/PairedEnd reads, or a ReadsSet input
string output_workspace - name or id of the WS to save the results to, provided by the narrative for housing output in KBase
string output_name - name of the output ReadsAlignment or ReadsAlignmentSet object
int runThreadN - the number of threads for STAR to use (default to 2)
string outFileNamePrefix: you can change the file prefixes using --outFileNamePrefix /path/to/output/dir/prefix
                        By default, this parameter is ./, i.e. all output files are written in current directory without a prefix
string quantMode: types of quantification requested--none/TranscriptomeSAM/GeneCounts
int outFilterMultimapNmax: max number of multiple alignments allowed for a read: if exceeded,
                        the read is considered unmapped, default to 20
int alignSJoverhangMin: minimum overhang for unannotated junctions, default to 8
int alignSJDBoverhangMin: minimum overhang for annotated junctions, default to 1
int outFilterMismatchNmax: maximum number of mismatches per pair, large number switches off this filter, default to 999
int alignIntronMin: minimum intron length, default to 20
int alignIntronMax: maximum intron length, default to 1000000
int alignMatesGapMax: maximum genomic distance between mates, default to 1000000
int create_report: = 1 if we build a report, 0 otherwise. (default 1) (shouldn not be user set - mainly used for subtasks)

@optional alignmentset_suffix
@optional alignIntronMin
@optional alignIntronMax
@optional alignMatesGapMax
@optional alignSJoverhangMin
@optional alignSJDBoverhangMin
@optional quantMode
@optional outFilterType
@optional outFilterMultimapNmax
@optional outSAMtype
@optional outSAMattrIHstart
@optional outSAMstrandField
@optional outFilterMismatchNmax
@optional outFileNamePrefix
@optional runThreadN


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
readsset_ref has a value which is a STAR.obj_ref
genome_ref has a value which is a STAR.obj_ref
output_workspace has a value which is a string
output_name has a value which is a string
alignment_suffix has a value which is a string
condition has a value which is a string
concurrent_njsw_tasks has a value which is an int
concurrent_local_tasks has a value which is an int
outSAMunmapped has a value which is a string
create_report has a value which is a STAR.bool
alignmentset_suffix has a value which is a string
alignIntronMin has a value which is an int
alignIntronMax has a value which is an int
alignMatesGapMax has a value which is an int
alignSJoverhangMin has a value which is an int
alignSJDBoverhangMin has a value which is an int
quantMode has a value which is a string
outFilterType has a value which is a string
outFilterMultimapNmax has a value which is an int
outSAMtype has a value which is a string
outSAMattrIHstart has a value which is an int
outSAMstrandField has a value which is a string
outFilterMismatchNmax has a value which is an int
outFileNamePrefix has a value which is a string
runThreadN has a value which is an int

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
readsset_ref has a value which is a STAR.obj_ref
genome_ref has a value which is a STAR.obj_ref
output_workspace has a value which is a string
output_name has a value which is a string
alignment_suffix has a value which is a string
condition has a value which is a string
concurrent_njsw_tasks has a value which is an int
concurrent_local_tasks has a value which is an int
outSAMunmapped has a value which is a string
create_report has a value which is a STAR.bool
alignmentset_suffix has a value which is a string
alignIntronMin has a value which is an int
alignIntronMax has a value which is an int
alignMatesGapMax has a value which is an int
alignSJoverhangMin has a value which is an int
alignSJDBoverhangMin has a value which is an int
quantMode has a value which is a string
outFilterType has a value which is a string
outFilterMultimapNmax has a value which is an int
outSAMtype has a value which is a string
outSAMattrIHstart has a value which is an int
outSAMstrandField has a value which is a string
outFilterMismatchNmax has a value which is an int
outFileNamePrefix has a value which is a string
runThreadN has a value which is an int


=end text

=back



=head2 AlignmentObj

=over 4



=item Description

Created alignment object returned.
ref = the workspace reference of the new alignment object
name = the name of the new object, for convenience.


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
ref has a value which is a STAR.obj_ref
name has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
ref has a value which is a STAR.obj_ref
name has a value which is a string


=end text

=back



=head2 AlignReadsResult

=over 4



=item Description

Here is the definition of the output of the function.  The output
can be used by other SDK modules which call your code, or the output
visualizations in the Narrative.  'report_name' and 'report_ref' are
special output fields- if defined, the Narrative can automatically
render your Report.

output_directory: folder path that holds all output files generated by run_star
alignmentset_ref: if an alignment set is created
alignment_objs: for each individual alignment created. The keys are the references to the reads
                object being aligned.
report_name: report name generated by KBaseReport
report_ref: report reference generated by KBaseReport


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
output_directory has a value which is a string
report_name has a value which is a string
report_ref has a value which is a STAR.obj_ref
alignmentset_ref has a value which is a STAR.obj_ref
alignment_objs has a value which is a reference to a hash where the key is a STAR.obj_ref and the value is a STAR.AlignmentObj

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
output_directory has a value which is a string
report_name has a value which is a string
report_ref has a value which is a STAR.obj_ref
alignmentset_ref has a value which is a STAR.obj_ref
alignment_objs has a value which is a reference to a hash where the key is a STAR.obj_ref and the value is a STAR.AlignmentObj


=end text

=back



=cut

package STAR::STARClient::RpcClient;
use base 'JSON::RPC::Client';
use POSIX;
use strict;

#
# Override JSON::RPC::Client::call because it doesn't handle error returns properly.
#

sub call {
    my ($self, $uri, $headers, $obj) = @_;
    my $result;


    {
	if ($uri =~ /\?/) {
	    $result = $self->_get($uri);
	}
	else {
	    Carp::croak "not hashref." unless (ref $obj eq 'HASH');
	    $result = $self->_post($uri, $headers, $obj);
	}

    }

    my $service = $obj->{method} =~ /^system\./ if ( $obj );

    $self->status_line($result->status_line);

    if ($result->is_success) {

        return unless($result->content); # notification?

        if ($service) {
            return JSON::RPC::ServiceObject->new($result, $self->json);
        }

        return JSON::RPC::ReturnObject->new($result, $self->json);
    }
    elsif ($result->content_type eq 'application/json')
    {
        return JSON::RPC::ReturnObject->new($result, $self->json);
    }
    else {
        return;
    }
}


sub _post {
    my ($self, $uri, $headers, $obj) = @_;
    my $json = $self->json;

    $obj->{version} ||= $self->{version} || '1.1';

    if ($obj->{version} eq '1.0') {
        delete $obj->{version};
        if (exists $obj->{id}) {
            $self->id($obj->{id}) if ($obj->{id}); # if undef, it is notification.
        }
        else {
            $obj->{id} = $self->id || ($self->id('JSON::RPC::Client'));
        }
    }
    else {
        # $obj->{id} = $self->id if (defined $self->id);
	# Assign a random number to the id if one hasn't been set
	$obj->{id} = (defined $self->id) ? $self->id : substr(rand(),2);
    }

    my $content = $json->encode($obj);

    $self->ua->post(
        $uri,
        Content_Type   => $self->{content_type},
        Content        => $content,
        Accept         => 'application/json',
	@$headers,
	($self->{token} ? (Authorization => $self->{token}) : ()),
    );
}



1;

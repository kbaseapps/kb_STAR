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
STAR-2.5.3a

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

  $output = $obj->run_star($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a STAR.STARParams
$output is a STAR.STARResults
STARParams is a reference to a hash where the following keys are defined:
	reads_ref has a value which is a string
	assembly_ref has a value which is a string
	genome_ref has a value which is a string
	workspace_name has a value which is a string
	runMode has a value which is a string
	runThreadN has a value which is an int
	genomeFastaFiles has a value which is a reference to a list where each element is a string
	sjdbGTFfile has a value which is a string
	sjdbOverhang has a value which is an int
	readFilesIn has a value which is a reference to a list where each element is a string
	outFileNamePrefix has a value which is a string
STARResults is a reference to a hash where the following keys are defined:
	reads_alignment_ref has a value which is a string
	report_name has a value which is a string
	report_ref has a value which is a string

</pre>

=end html

=begin text

$params is a STAR.STARParams
$output is a STAR.STARResults
STARParams is a reference to a hash where the following keys are defined:
	reads_ref has a value which is a string
	assembly_ref has a value which is a string
	genome_ref has a value which is a string
	workspace_name has a value which is a string
	runMode has a value which is a string
	runThreadN has a value which is an int
	genomeFastaFiles has a value which is a reference to a list where each element is a string
	sjdbGTFfile has a value which is a string
	sjdbOverhang has a value which is an int
	readFilesIn has a value which is a reference to a list where each element is a string
	outFileNamePrefix has a value which is a string
STARResults is a reference to a hash where the following keys are defined:
	reads_alignment_ref has a value which is a string
	report_name has a value which is a string
	report_ref has a value which is a string


=end text

=item Description

The actual function is declared using 'funcdef' to specify the name
and input/return arguments to the function.  For all typical KBase
Apps that run in the Narrative, your function should have the 
'authentication required' modifier.

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



=head2 boolean

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



=head2 STARParams

=over 4



=item Description

Arguments for star_generate_indexes

string reads_ref, assembly_ref and genome_ref: KBase style variable references
string runMode: default: alignReads
        type of the run:
        alignReads => map reads
        genomeGenerate => generate genome files
        inputAlignmentsFromBAM => input alignments from BAM. Presently only works with -outWigType
                and -bamRemoveDuplicates.
        liftOver => lift-over of GTF files (-sjdbGTFfile) between genome assemblies using
                chain file(s) from -genomeChainFiles.
int runThreadN: default: 1
        number of threads to run STAR
list<string> genomeFastaFiles: path(s) to the fasta files with genomic sequences for genome generation. 
        Only used if runMode==genomeGenerate.These files should be plain text FASTA files, they *cannot* be zipped.
list<string> readFilesIn: default: Read1 Read2
        paths to files that contain input read1 (and, if needed, read2)

string sjdbGTFfile: default: -; path to the file with annotated transcripts in the standard GTF format
int sjdbOverhang: default: 100; int>0: length of the donor/acceptor sequence on each side of the junctions,
        ideally = (ReadLength - 1)
string outFileNamePrefix: you can change the file prefixes using --outFileNamePrefix /path/to/output/dir/prefix.
        By default, this parameter is ./, i.e. all output files are written in the current directory

@optional sjdbGTFfile
@optional sjdbOverhang


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
reads_ref has a value which is a string
assembly_ref has a value which is a string
genome_ref has a value which is a string
workspace_name has a value which is a string
runMode has a value which is a string
runThreadN has a value which is an int
genomeFastaFiles has a value which is a reference to a list where each element is a string
sjdbGTFfile has a value which is a string
sjdbOverhang has a value which is an int
readFilesIn has a value which is a reference to a list where each element is a string
outFileNamePrefix has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
reads_ref has a value which is a string
assembly_ref has a value which is a string
genome_ref has a value which is a string
workspace_name has a value which is a string
runMode has a value which is a string
runThreadN has a value which is an int
genomeFastaFiles has a value which is a reference to a list where each element is a string
sjdbGTFfile has a value which is a string
sjdbOverhang has a value which is an int
readFilesIn has a value which is a reference to a list where each element is a string
outFileNamePrefix has a value which is a string


=end text

=back



=head2 STARResults

=over 4



=item Description

Here is the definition of the output of the function.  The output
can be used by other SDK modules which call your code, or the output
visualizations in the Narrative.  'report_name' and 'report_ref' are
special output fields- if defined, the Narrative can automatically
render your Report.


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
reads_alignment_ref has a value which is a string
report_name has a value which is a string
report_ref has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
reads_alignment_ref has a value which is a string
report_name has a value which is a string
report_ref has a value which is a string


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

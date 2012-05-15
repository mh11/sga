#! /usr/bin/perl
use strict;
use Carp;
use List::Util qw[min];
use Getopt::Long;
use Pod::Usage;
use threads;
use Thread::Queue qw( );
use Thread::Semaphore;
use File::Basename;
use threads::shared;

my $help = 0;
my $num_threads = 8;
my $output_prefix = "final";
my $no_reverse = 0;
my $no_index = 0;
my $no_merge = 0;
my $dry_run = 0;

my $error :shared;
$error = 0;


GetOptions('threads=i'    => \$num_threads,
           'no-reverse'   => \$no_reverse,
           'no-index'   => \$no_index,
           'no-merge'   => \$no_merge,
           'dry-run'   => \$dry_run,
	   'help|?' => \$help
) or pod2usage(2);

pod2usage(1) if $help;

=head1 NAME

  sga-pims.pl - Parallel Indexing and Merging Script

=head1 SYNOPSIS

sga-pims.pl [options] [files ...]

=head1 OPTIONS

=over 8

=item B<--help>
Print a brief help message and exits.

=item B<--threads>
threads to use [default: 8] 

=item B<--dry-run>
don't execute commands - just print them

=item B<--no-reverse>
don't build reverse index_jobs

=item B<--no-index>
don't run indexing

=item B<--no-merge>
don't run index merge

=back

=cut

my @files = @ARGV;

my @fList = @files;

if(!$no_index){
  @fList = indexFiles(@fList);
}
if(!$no_merge){
  my $lastFile = mergeFiles(@fList);
  renameFiles($lastFile,$output_prefix);
}

print "Finished processing files !!! \n" ;

exit 0;

sub indexFiles
{
  my @files = @_;
  if(scalar(@files) == 0){
    print STDERR "No files to merge !!!\n";
    exit 1;
  }
  my $q = Thread::Queue->new();
  my @threads;
  my $index_jobs = min($num_threads, scalar(@files));

  print "Index files ",scalar(@files)," files using ",$index_jobs," threads ... \n";

  $q->enqueue(@files);

  for (1..$index_jobs) {
    push @threads, async {
      sleep 1; # to wait for other threads to start up
      my $thr = threads->self();
      while (!$error && (my $item = $q->dequeue_nb())) {
	my $threads_running = scalar(threads->list(threads::running))+1;
	# Calculate the amount of threads free to use
	my $threads_to_use = 0;
	{ use integer;
	  $threads_to_use = $num_threads/$threads_running;
	}

	printf("Thread[%d] indexing %s using %d threads ... \n", $thr->tid(), $item, $threads_to_use);
	my $no_reverse_opt = $no_reverse ? "--no-reverse" : "";
	run("sga index -t $threads_to_use -a BCR -d 5000000 $no_reverse_opt $item");
      }
      printf("Thread[%d] finished!\n", $thr->tid());
      return 0;
    };
  }
  do_shutdown(@threads);
  return @files;
}

sub mergeFiles
{
  my @files_to_merge = @_;
  if(scalar(@files_to_merge) == 0){
    print STDERR "No files to merge !!!\n";
    exit 1;
  } elsif(scalar(@files_to_merge) ==1){
    print STDERR "Nothing to merge with one file ",$files_to_merge[0],"!!!\n";
    return $files_to_merge[0];
  }

  my $q = Thread::Queue->new();
  my @threads;
  my $merge_jobs = min($num_threads, scalar(@files_to_merge));

  print "Merge ",scalar(@files_to_merge)," files using ",$merge_jobs," threads ...\n";
 
  $q->enqueue(@files_to_merge);

  for (1..$merge_jobs) {
    push @threads, async {
	sleep 1; # to wait for other threads to start up
	my $thr = threads->self();
	my $i_merged = 0;
	# if one or no file is pending, stop Thread -> requires two files to merge
	while (!$error && ($q->pending() > 1) && (my @arr = $q->dequeue_nb(2))) {
	  if(scalar(@arr)== 2){
	    my $out = sprintf("%s.merged.%d_%d", $output_prefix, $thr->tid(),$i_merged++);
	    my $outf = "$out.fa";
	    my $threads_running = scalar(threads->list(threads::running))+1;

	    # Calculate the amount of threads free to use
	    my $threads_to_use = 0;
	    { use integer;
	      $threads_to_use = $num_threads/$threads_running;
	    }
	    printf("Thread[%d] merging %s into %s using %d threads with %d running ... \n", $thr->tid(), join(" ",@arr),$outf,$threads_to_use,$threads_running);
	    
	    # run SGA merge
	    run("sga merge -r -t $threads_to_use -p $out $arr[0] $arr[1] ");
	    printf("Thread[%d] merging into %s done.\n", $thr->tid(), $outf);

	    # push the result file back onto the queue
	    $q->enqueue($outf);

	  } else {
	    # put it back on the queue, if there is only one file left - try luck again next time round.
	    $q->enqueue(@arr);
	    printf("Thread[%d] requeue %d items!\n", $thr->tid(),scalar(@arr));
	  }
	}
	printf("Thread[%d] finished!\n", $thr->tid());
      return 0;
    };
  }
  print "[Main] Waiting for ",scalar(@threads)," threads to finish ... \n";

  do_shutdown(@threads);

  if($q->pending() != 1){
    print STDERR "Unexpected ",$q->pending()," files left: ",join(",",$q->dequeue_nb($q->pending())),"\n";
    exit 1;
  }
  my $return_file = $q->dequeue_nb();
  return $return_file;
}

#rename all fa and index files files 
sub renameFiles
{
  my ($curr,$target) = @_;
  print "Rename ",$curr," files to prefix $target ... \n";
  my ($file, $path, $suffix) = fileparse($curr, (".fa", ".fastq",".fq"));
  my $inbase = $path . $file;
  my $outbase = $path . $target;

  foreach my $ext (".sai", ".bwt", $suffix)
  {
      run("mv $inbase$ext $outbase$ext");
  }
  if(!$no_reverse){
    foreach my $ext (".rsai", ".rbwt")
    {
	run("mv $inbase$ext $outbase$ext");
    }
  }
}

sub do_shutdown
{
  my @threads = @_;
  my $hasError = 0;
  foreach my $t (@threads){
    my @t_ret = $t->join();
    # expect 0 as return value
    if($t_ret[0] ne 0){
      print STDERR "Thread ",$t->tid()," returned an unexpected value: ",join(",",@t_ret),"\n";
      $hasError = 1;
    }
  }
  # if any error reported or found
  if($hasError || $error){	
    print STDERR "Problem detected during execution -> Exit!!!\n";
    exit 1;
  }
}

# Run a command                                                                                                                                              
sub run
{
    my($cmd) = @_;
    print $cmd . "\n";
	  sleep 1;
    if($dry_run){
	  sleep 1;
    } else {
      my $returnValue = system($cmd);
      if($returnValue != 0){
	$error = 1;
        croak("Failed to execute >$cmd<: $!\n");	
      }
    }
}



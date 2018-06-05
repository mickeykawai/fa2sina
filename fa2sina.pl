#!/usr/bin/env perl

###############################################################################
#
#    fa2sina.pl
#
#	 See usage by calling fa2sina.pl without arguments. 
#    
#    Copyright (C) 2018 Mikihiko Kawai
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

use strict;
use warnings 'all';

use CGI;
#use Getopt::Long;

use File::Basename;
use Data::Dumper;
use File::Path;
use Carp;
use FileHandle;

use Cwd;

sub ___________________Object_oriented___________________{
}

sub _____________________Base___________________{
}
sub _write_conf_make{
		#$e = $obj->_write_conf_make($f_conf_make, $e);
	my $obj = shift;
	my $f_conf_make = shift;
	my $e = shift;
	my $key_list_cmd = shift || [];	#cmd_init, cmd_link etc (to permit multiple cmds)
	
	my $date_time = date_time(5);
	
	my $key_list = [qw(macro_param macro_cmdname macro_cmd dependency_target)];
	
	my $seen_ks = {};
	if (scalar @$key_list_cmd){
		foreach my $k (@$key_list_cmd){
			if ($k !~ /^cmd_/){
				confess "given key list for cmds must be 'cmd_<cmd_type>'. given '$k' cannot be accepted. \n";
			}
			push @$key_list, $k if !$seen_ks->{$k}++;
		}
		#push @$key_list, @$key_list_cmd;
	}
	
	my $con_dependency_target = $e->{dependency_target};
	my $es_dependency_target = get_hash($con_dependency_target, undef, {skip_empty_line_flag => 1});
	my $ks_cmd_list1 = [map {'cmd_' . $_->{cmd_type}} @$es_dependency_target];
	foreach my $k (@$ks_cmd_list1){
		if ($e->{$k}){
				#e.g. cmd_init
			push @$key_list, $k if !$seen_ks->{$k}++;
		}
	}
	
	my $con_conf_make = << "CON_CONF_MAKE";
##
 This is a configuration file to make a Makefile by the following command:
 conf_make2makefile.pl in=<this_conf_make_file> (out=<Makefile>)
 Then use the Makefile by the following way as standard make:
 cd <dir_of_the_Makefile>;make (or equivalently make -f <makefile>).
 The format of this configuration file itself is specifically defined (not general format), 
 but Makefile follows the format of gnu makefile. 
 
 This file is created on $date_time by script '$0'. 
##
//
CON_CONF_MAKE
	$con_conf_make =~ s/^/#/mg;
	
	foreach my $key (@$key_list){
	#foreach my $key (qw(macro_param macro_cmdname macro_cmd dependency_target)){
			#$e->{macro_cmd} = \$con_macro_cmd;
			#can be with empty line
		if ($e->{$key}){
			$con_conf_make .= << "CON_CONF_MAKE";
##
# $key
##
${$e->{$key}}
#//
CON_CONF_MAKE
		}
	}
	
	mkpath(dirname($f_conf_make));
	open (OUT, ">$f_conf_make") or confess "cannot open file '$f_conf_make': $!";
	print OUT $con_conf_make;
	close (OUT) or confess "cannot close file '$f_conf_make': $!";
	#system "chmod 0775 $f_sh";
	
	#$e->{content_conf_make} = $con_conf_make;
	
	return $e;
}

sub parse_fa_sina{
		#$obj->parse_fa_sina();
	my $obj = shift;
	
#	my $fa = $obj->{fa};
#	my $arb = $obj->{arb};
	my $lca_field = $obj->{lca_field};
	my $out = $obj->{out};
	my $is_name_till_semicolon = $obj->{is_name_till_semicolon};
	
	my $f_tax = $out . '.tax.txt';
	my $f_log_sina = $out . '.log_sina.txt';	#this is parsed log of sina output fasta, not log file of sina itself. 
#	if (
#		-e $f_tax and -s $f_tax and -e $f_log_sina and -s $f_log_sina
#			and
#		(-f $out and -s $out and (stat $out)[9] <= (stat $f_tax)[9])
#	){
#		warn "tax list file '$f_tax' and sina parsed log file '$f_log_sina' exists, skip. \n";
#	}else{
		open (IN, $out) || die "cannot open $out: $!";
		open (OUT_TAX, ">$f_tax") or confess "cannot open $f_tax: $!";
		open (OUT_LOG, ">$f_log_sina") or confess "cannot open $f_log_sina: $!";
		#my $i = -1;
		
		my $ks_tax = [qw(name lca_field tax align_ident_slv low_confident)];
		my $str_k_tax = join "\t", @$ks_tax;
		print OUT_TAX "$str_k_tax\n";
		
		my $ks_log = [qw(name)];
		my $seen_ks = {};
		my $is_to_out_ks_log = 1;
		while (<IN>){
			if (/^>(.+)$/){
				chomp;
				#$i++;
				my $fasta_name = $1;
				
				my $e = {
					name => $fasta_name,
				};
				$e->{name_full} = $e->{name};
				
				#my $pairs = [
				my $pairs = [];
				push @$pairs, $1, $2 while (
					$e->{name_full} =~ /
						\[
							([^=]+?)[=](
								(?&WORD) | (?&BRACKETED)
								(?(DEFINE)
									(?<WORD>      \s* [^\[\]]+ )
									(?<BRACKETED> \s* \[ (?&TEXT)? \s* \] )
									(?<TEXT>      (?: (?&WORD) | (?&BRACKETED) )+ )
								)
							)
						\]
					/xg
				);
				#];
					#[nearest_slv=ARB_25BABFAF...~0.951049 X89043.[EBI] X89043.1..~0.941725 ARB_C313EA80...~0.939394 ARB_14293161...~0.937063 ARB_8FADD83B...~0.934732 ARB_1DEA9552...~0.934732 ARB_C41E8E79...~0.927739 X79548.[EBI] X79548.1..~0.925408 ARB_72F20AB7...~0.925408 ARB_B2DE4B46...~0.923077 ] [nuc=429] [nuc_gene_slv=429] [turn_slv=none]
				#my $pairs = [$e->{name_full} =~ /\[([^=]+?)[=](.*?)\]/g];
				#print STDERR Dumper $pairs;die;#;stop;

#push @$pairs, $1
#	while /
#		\G \s*+ ( (?&WORD) | (?&BRACKETED) )
#		
#		(?(DEFINE)
#			(?<WORD>      \s* \w+ )
#			(?<BRACKETED> \s* \[ (?&TEXT)? \s* \] )
#			(?<TEXT>      (?: (?&WORD) | (?&BRACKETED) )+ )
#		)
#	/xg;

				my $h = {};
				for (my $i = 0; $i <= $#$pairs; $i += 2){
					my ($k, $v) = ($pairs->[$i], $pairs->[$i+1]);
					$h->{$k} = $v;
					push @$ks_log, $k if !$seen_ks->{$k}++;
				}
				#$h = {$e->{name_full} =~ /\[([^=]+?)[=](.+?)\]/g};
				
				$e->{name} =~ s/ .+//;
				$e->{name} =~ s/;.+// if $is_name_till_semicolon;	#Hak16S1_1;size=9279; -> Hak16S1_1
				$e->{lca_field} = $lca_field;
				$e->{tax} = $h->{'lca_' . $lca_field};	#lca_tax_slv, lca_tax_yh, ...
				$e->{align_ident_slv} = $h->{align_ident_slv};
				if ($e->{align_ident_slv} < 95 and $e->{align_ident_slv} >= 90){
					$e->{low_confident} = '-';
				}
				elsif ($e->{align_ident_slv} < 90 and $e->{align_ident_slv} >= 80){
					$e->{low_confident} = '--';
				}
				elsif ($e->{align_ident_slv} < 80){
					$e->{low_confident} = '---';
				}
				#$e->{align_startpos_slv} = $h->{align_startpos_slv};
				#$e->{align_stoppos_slv} = $h->{align_stoppos_slv};
				#[align_startpos_slv=988] [align_stoppos_slv=43283]
					#already recorded in log file as tab format.
				
				$h->{name} = $e->{name};
				#push @$es_log, $h;
				#print STDERR Dumper $e,$h,$ks_log,$ks_tax;die;#;stop;
				#my $f_tax = $out . '.tax.txt';
				foreach my $k (@$ks_tax){
					defined $e->{$k} or $e->{$k} = '';
				}
				my $str_tax = join "\t", @{$e}{@$ks_tax};
				print OUT_TAX "$str_tax\n";
				
				#my $f_log_sina = $out . '.log_sina.txt';
				foreach my $k (@$ks_log){
					defined $h->{$k} or $h->{$k} = '';
				}
				if ($is_to_out_ks_log){
					my $str_k_log = join "\t", @$ks_log;
					print OUT_LOG "$str_k_log\n";
					$is_to_out_ks_log = 0;
				}
				my $str_log = join "\t", @{$h}{@$ks_log};
				print OUT_LOG "$str_log\n";
			}
			else{
			}
		}
		close (OUT_TAX) or confess "cannot close $f_tax: $!";
		close (OUT_LOG) or confess "cannot close $f_log_sina: $!";
		close (IN) || die "cannot close $out: $!";
	
	return;
}

sub summary_tax{
		#$obj->summary_tax();
	my $obj = shift;
	
	my $fa = $obj->{fa};
	my $arb = $obj->{arb};
#	my $lca_field = $obj->{lca_field};
	my $out = $obj->{out};
	#my $is_sum_size = $obj->{is_sum_size};
	my $is_name_till_semicolon = $obj->{is_name_till_semicolon};
	
	mkpath(dirname($out));
	
	my $search_min_sim = 0.80;	#for sina
	my $idents_check = [80, 90, 95, 75];	#for summary
		#must be ascending order. 
	my $main_threshold_ident = $idents_check->[0];	#80; for summary
	#my $main_threshold_ident = 80;	#for summary
	
	
	my $f_tax = $out . '.tax.txt';
#	my $f_log_sina = $out . '.log_sina.txt';	#this is parsed log of sina output fasta, not log file of sina itself. 
	
	my $f_summary = $out . '.tax.summary.txt';
#	if (
#		-e $f_summary and -s $f_summary
#			and
#		(-f $f_tax and -s $f_tax and (stat $f_tax)[9] <= (stat $f_summary)[9])
#	){
#	#if (0 and -e $f_summary and -s $f_summary){
#		warn "tax summary file '$f_summary' exists, skip. \n";
#	}else{
		#my $idents_check = [qw(80, 90, 95)];
		my $idents = [];
		my $is_same_found = 0;
		my $ident0 = int($search_min_sim * 100);
		foreach my $ident (@$idents_check){
			if ($ident == $ident0){
				$is_same_found = 1;
				#last;
			}
			#if ($ident >= $ident0){	#in a few cases, there are hits with ident lower than search_min_sim are detected. 
				push @$idents, $ident;
			#}
		}
		if (!$is_same_found){
			push @$idents, int($search_min_sim * 100);
		}
		
		open (IN, $f_tax) || die "cannot open $f_tax: $!";
		my $ks = [];
		my $is_first = 1;
		my $tax2n = {};
		while (<IN>){
			next if !$is_first and /^#/;
			
			chomp;
			#warn "$_\n";
			if ($is_first){
				$ks = [split /\t/, $_];
				$is_first = 0;
			}else{
				my $els = [split /\t/, $_];
				my $e = {};
				@{$e}{@$ks} = @$els;
				
				foreach my $threshold_ident (@$idents){
					#foreach my $e (@$es){
					if ($e->{align_ident_slv} >= $threshold_ident){
						#if ($is_sum_size){
						
						#}else{
							$tax2n->{$threshold_ident}{$e->{tax}}++;
						#}
					}
						#$tax2n->{$e->{tax}}++ if $e->{align_ident_slv} >= $threshold_ident;
					#}
				}
				#print STDERR Dumper $e,$tax2n;die;#;stop;
			}
		}
		close (IN) || die "cannot close $f_tax: $!";
		#print STDERR Dumper $idents;die;#;stop;
		my $ks_count = [qw(tax count)];
		foreach my $threshold_ident (@$idents){
			warn "$threshold_ident..\n";
			my $es_s = [];
			foreach my $k (sort keys %{$tax2n->{$threshold_ident}}){
				push @$es_s, {tax => $k, count => $tax2n->{$threshold_ident}{$k}};
			}
			my $f_summary_w_ident = $out . '.tax.summary' . ".ident.$threshold_ident.txt";
			#my $es_s1 = [grep {$_->{align_ident_slv} >= $threshold_ident} @$es_s];
			output_hash_my($f_summary_w_ident, $es_s, K => $ks_count);
			#output_hash_my($f_summary_w_ident, $es_s1, K => $ks);
		}
		
		#my $main_threshold_ident = 80;
		##my $main_threshold_ident = 95;
		my $f_summary_w_ident = $out . '.tax.summary' . ".ident.$main_threshold_ident.txt";
		#my $f_summary = $out . '.tax.summary.txt';
		
		my $wd = dirname($f_summary);
		my $base_summary_w_ident = basename($f_summary_w_ident);
		my $base_summary = basename($f_summary);
		system "cd $wd;ln -s -f $base_summary_w_ident $base_summary";
		#system "ln -s -f $f_summary_w_ident $f_summary";
		#output_hash_my(\*STDERR, $es_s, K => $ks);
#	}
	
	my ($ks1, $es_s) = get_hash_my($f_summary);
		#my $ks = [qw(tax count)];
	output_hash_my(\*STDERR, $es_s, K => $ks1);
	
#	my $fa_arb = $out . '.arb.fa';
#	if (-e $fa_arb and -s $fa_arb){
#		warn "arb-compatible output file '$fa_arb' exists, skip. \n";
#	}else{
#		my $cmd_sina2arb = "sina_fas2arb_fas.pl fa=$out out=$fa_arb";
#		warn "\n", "when you load the sequences on ARB as pre-aligned, use following file.\n",
#			"(gap characters are converted to '.' from '-', to be compatible with aligned ARB, at 5' and 3' ends.)",
#			"$cmd_sina2arb\n";
#		system "$cmd_sina2arb";
#	}
	
	return;
}

sub sina_fas2arb_fas{
		#$obj->sina_fas2arb_fas();
	my $obj = shift;
	
	my $fa = $obj->{fa};
	my $arb = $obj->{arb};
	my $lca_field = $obj->{lca_field};
	my $out = $obj->{out};
	#my $is_sum_size = $obj->{is_sum_size};
	my $is_name_till_semicolon = $obj->{is_name_till_semicolon};
	
	mkpath(dirname($out));
	
	my $fa_arb = $out . '.arb.fa';
#	if (-e $fa_arb and -s $fa_arb){
#		warn "arb-compatible output file '$fa_arb' exists, skip. \n";
#	}else{
		my $cmd_sina2arb = "sina_fas2arb_fas.pl fa=$out out=$fa_arb";
		warn "\n", "when you load the sequences on ARB as pre-aligned, use following file.\n",
			"(gap characters are converted to '.' from '-', to be compatible with aligned ARB, at 5' and 3' ends.)",
			"$cmd_sina2arb\n";
		system "$cmd_sina2arb";
#	}
	
	return;
}
sub fa2sina{
		#$obj->fa2sina();
	my $obj = shift;
	
	my $class = ref $obj;
	(my $class_str = lc $class) =~ s/::/__/g;
	$obj->init_parameters_default(
		"set_obj_key_list__$class_str" => [qw(
		)], 
		"set_obj_key_list_opt__$class_str" => [qw(
			search_min_sim
			idents_check
			main_threshold_ident
			
			is_renew
			is_only_set_params
			make_targets
		)], 
		search_min_sim => 0.80,
		idents_check => '80,90,95,75',	#for summary
			#must be ascending order. 
		main_threshold_ident => '80',
		
		is_renew => 0, 
		is_only_set_params => 0,
		
		make_targets => 'sina parse_fa_sina summary_tax sina_fas2arb_fas',
			#space-delimited (not comma-delimited)
	);	#in 
	
	my $fa = $obj->{fa};
	my $arb = $obj->{arb};
	my $lca_field = $obj->{lca_field};
	my $out = $obj->{out};
	#my $is_sum_size = $obj->{is_sum_size};
	my $is_name_till_semicolon = $obj->{is_name_till_semicolon};
	my $is_renew = $obj->{is_renew};
	
	if (!-e $fa or -z $fa){
		warn "query not exists, skip.\n";
		return;
	}
	
	my $cmdname_fa2sina = "fa2sina.pl";
	
	my $outdir = dirname($out);
	mkpath($outdir);
	#mkpath(dirname($out));
	
	my $search_min_sim = $obj->{search_min_sim};	#for sina
	my $ident_check_str = $obj->{idents_check};	#for summary
		#must be ascending order. 
	my $idents_check = [split /,/, $ident_check_str];
	my $main_threshold_ident = $obj->{main_threshold_ident};	#80; for summary
	#my $main_threshold_ident = 80;	#for summary
#	my $search_min_sim = 0.80;	#for sina
#	my $idents_check = [80, 90, 95, 75];	#for summary
#		#must be ascending order. 
#	my $main_threshold_ident = $idents_check->[0];	#80; for summary
#	#my $main_threshold_ident = 80;	#for summary
	
	my $logfile_sina = $out . '.sina.log.txt';
		#  --log-file arg (=/dev/stderr) log to file (default: STDERR)
		#to avoid error:
		#Fatal error: Unable to open file "/dev/stderr" for writing.
	
	my $f_tax = $out . '.tax.txt';
	my $f_log_sina = $out . '.log_sina.txt';	#this is parsed log of sina output fasta, not log file of sina itself. 
	
	my $f_summary = $out . '.tax.summary.txt';
	
	my $fa_arb = $out . '.arb.fa';
	
	if (-e $out and -s $out){
		warn "sina output file '$out' exists, skip. \n";
	}else{
		warn "sina output file '$out' not exists, do sina.. \n";
	}
		#my $cmd_sina = "sina --search-min-sim $search_min_sim -i $fa -o $out --ptdb $arb --search --search-db $arb --lca-fields $lca_field --meta-fmt header --intype fasta --log-file $logfile_sina";
		##my $cmd_sina = "sina --search-min-sim $search_min_sim -i $fa -o $out --ptdb $arb --search --search-db $arb --lca-fields $lca_field --meta-fmt header --intype fasta";
		##my $cmd_sina = "sina --search-min-sim 0.80 -i $fa -o $out --ptdb $arb --search --search-db $arb --lca-fields $lca_field --meta-fmt header --intype fasta";
		##my $cmd_sina = "sina --search-min-sim 0.90 -i $fa -o $out --ptdb $arb --search --search-db $arb --lca-fields $lca_field --meta-fmt header --intype fasta";
		#	#--search-min-sim arg
		##my $cmd_sina = "sina --search-all --search-no-fast -i $fa -o $out --ptdb $arb --search --search-db $arb --lca-fields $lca_field --meta-fmt header --intype fasta";
		#	#--search-all
		#	#--search-no-fast
		##my $cmd_sina = "sina --lca-quorum 0.5 -i $fa -o $out --ptdb $arb --search --search-db $arb --lca-fields $lca_field --meta-fmt header --intype fasta";
		#	#--lca-quorum 0.5
		##my $cmd_sina = "sina -i $fa -o $out --ptdb $arb --search --search-db $arb --lca-fields $lca_field --meta-fmt header --intype fasta";
		##my $cmd_sina = "sina -i $fa -o $out --ptdb $arb --search --search-db $arb --lca-fields $lca_field --meta-fmt header --intype fasta > $fa.log.sina.txt 2>&1";
		#	#Fatal error: Unable to open file "/dev/stderr" for writing.
		##my $cmd_sina = "sina -i $fa -o $out --ptdb $arb --search --search-db $arb --lca-fields $lca_field --meta-fmt header --intype fasta 2>&1 | tee $fa.log.sina.txt";
		##my $cmd_sina = "sina -i $fa -o $out --ptdb $arb --search --search-db $arb --lca-fields $lca_field --meta-fmt csv --intype fasta | tee $fa.log.sina.txt";
		##sina -i $f -o all_prinseq_good_igMu.fa.160607_0206ssuARB_HNG2_MK3.fasta --ptdb $arb --search --search-db $arb --lca-fields tax_yh --meta-fmt csv --intype fasta | tee log.arb.morikawa.txt
		##--ptport :/tmp/sina_pt_68329 --search-port :/tmp/sina_pt_68329
	
#	mkpath($outdir);
	my $fa_base = basename($fa);
	
	my $wd = $outdir;
	my $f_sh = join '/', $wd, 'fa2sina' . '.' . $fa_base . '.sh';
	
#	my $wd = join '/', $e->{wd};#, $opt_base_prokka;
#	mkpath($wd);
#	$f_sh = join '/', $wd, basename($f_sh);
	
	my $cmd;
	CONF_MAKE2MAKECONF:
	{
		my $f_checkpoint_start = $f_sh . '.start';
		my $f_checkpoint = $f_sh . '.end';
		if ($is_renew and -f $f_checkpoint_start){
			system "touch $f_checkpoint_start";
		}
		
		( my $base_program = basename($0) ) =~ s/\..+$//;;
			#fa2sina
		#my $base_program = $type;
		my $f_conf_make = join '/', $wd, $base_program . '.conf.make.txt';
		my $f_makefile = join '/', $wd, 'Makefile.' . $base_program;
		
		my $e0 = {};
		my $con_key2param = << "CON_KEY2PARAM_ENTRIES";
key						param							is_rm_empty
wd						${wd}

f_checkpoint_start		${f_checkpoint_start}			0
f_checkpoint			${f_checkpoint}					0

f_fa					$fa					
arb						$arb
lca_field				$lca_field

is_name_till_semicolon	$is_name_till_semicolon
search_min_sim			$search_min_sim
idents_check			$ident_check_str
main_threshold_ident	$main_threshold_ident

f_sina					$out							1
f_log_sina				$logfile_sina					1

f_tax					$f_tax							1
f_log_sina				$f_log_sina						1

f_summary				$f_summary						1
f_arb					$fa_arb							1

cmdname_fa2sina			$cmdname_fa2sina
CON_KEY2PARAM_ENTRIES
#name					$e->{name}
		
		$con_key2param =~ s/\t+/\t/g;	#mutliple tabs in the content above between cells for readablity. trim them to be one tab. 
		if ($obj->{is_only_set_params}){
			print $con_key2param;
			exit(0);
		}
		$e0->{macro_param} = \$con_key2param;
		
		my $con_dependency_target;
		$con_dependency_target = << 'CON_CMD_ENTRIES';
cmd_type			dependency_str						target_str						cmd
init				f_fa								f_checkpoint_start				touch $(f_checkpoint_start)
#link				f_fa_given,f_checkpoint_start		f_fa							ln -s -f $(f_fa_given) $(f_fa)

sina				f_fa,f_checkpoint_start				f_sina							sina --search-min-sim $(search_min_sim) -i $(f_fa) -o $(f_sina) --ptdb $(arb) --search --search-db $(arb) --lca-fields $(lca_field) --meta-fmt header --intype fasta --log-file $(f_log_sina)
parse_fa_sina		f_sina								f_tax,f_log_sina				$(cmdname_fa2sina) methods=parse_fa_sina out=$(f_sina) fa=$(f_fa) arb=$(arb) is_name_till_semicolon=$(is_name_till_semicolon)
summary_tax			f_tax								f_summary						$(cmdname_fa2sina) methods=summary_tax out=$(f_sina) fa=$(f_fa) arb=$(arb) search_min_sim=$(search_min_sim) main_threshold_ident=$(main_threshold_ident) idents_check=$(idents_check)
sina_fas2arb_fas	f_sina								f_arb							$(cmdname_fa2sina) methods=sina_fas2arb_fas out=$(f_sina) fa=$(f_fa) arb=$(arb) 

end		f_sina,f_tax,f_log_sina,f_summary,f_arb			f_checkpoint					touch $(f_checkpoint)
CON_CMD_ENTRIES
		
		$con_dependency_target =~ s/\t+/\t/g;	#mutliple tabs in the content above between cells for readablity. trim them to be one tab.
		$con_dependency_target =~ s/#dummy#//g;
		$e0->{dependency_target} = \$con_dependency_target;
		
#				$e0->{cmd_init} = \ << 'CON_CMD';
#. ~/.bash_profile
#mkdir -p $(wd)
#cd $(wd)
#touch $(f_checkpoint_start)
#CON_CMD
		
		$e0 = $obj->_write_conf_make($f_conf_make, $e0);
		
		my $cmd_conf_make2makefile = "conf_make2makefile.pl in=$f_conf_make out=$f_makefile";
		warn "$cmd_conf_make2makefile\n\n";
		system "$cmd_conf_make2makefile";
		
		$cmd = "cd $wd;make -f $f_makefile";
		if ($obj->{make_targets}){
			$cmd .= " $obj->{make_targets}";
		}
	};
	
	open (OUT, ">$f_sh") or confess "cannot open file '$f_sh': $!";
	print OUT $cmd;
	close (OUT) or confess "cannot close file '$f_sh': $!";
	system "chmod 0775 $f_sh";
	
	warn "f_sh1 '$f_sh'\n";
	my $cmd11 = "sh $f_sh";
	
	system "$cmd11";
	
	return;
}

sub fa2sina_ver_sh{
		#$obj->fa2sina_ver_sh();
	my $obj = shift;
	
	my $fa = $obj->{fa};
	my $arb = $obj->{arb};
	my $lca_field = $obj->{lca_field};
	my $out = $obj->{out};
	#my $is_sum_size = $obj->{is_sum_size};
	my $is_name_till_semicolon = $obj->{is_name_till_semicolon};
	
	mkpath(dirname($out));
	
	my $search_min_sim = 0.80;	#for sina
	my $idents_check = [80, 90, 95, 75];	#for summary
		#must be ascending order. 
	my $main_threshold_ident = $idents_check->[0];	#80; for summary
	#my $main_threshold_ident = 80;	#for summary
	
	if (!-e $fa or -z $fa){
		warn "query not exists, skip.\n";
		return;
	}
	
	if (-e $out and -s $out){
		warn "sina output file '$out' exists, skip. \n";
	}else{
		warn "sina output file '$out' not exists, do sina.. \n";
		my $logfile_sina = $out . '.sina.log.txt';
			#  --log-file arg (=/dev/stderr) log to file (default: STDERR)
			#to avoid error:
			#Fatal error: Unable to open file "/dev/stderr" for writing.
		my $cmd_sina = "sina --search-min-sim $search_min_sim -i $fa -o $out --ptdb $arb --search --search-db $arb --lca-fields $lca_field --meta-fmt header --intype fasta --log-file $logfile_sina";
		#my $cmd_sina = "sina --search-min-sim $search_min_sim -i $fa -o $out --ptdb $arb --search --search-db $arb --lca-fields $lca_field --meta-fmt header --intype fasta";
		#my $cmd_sina = "sina --search-min-sim 0.80 -i $fa -o $out --ptdb $arb --search --search-db $arb --lca-fields $lca_field --meta-fmt header --intype fasta";
		#my $cmd_sina = "sina --search-min-sim 0.90 -i $fa -o $out --ptdb $arb --search --search-db $arb --lca-fields $lca_field --meta-fmt header --intype fasta";
			#--search-min-sim arg
		#my $cmd_sina = "sina --search-all --search-no-fast -i $fa -o $out --ptdb $arb --search --search-db $arb --lca-fields $lca_field --meta-fmt header --intype fasta";
			#--search-all
			#--search-no-fast
		#my $cmd_sina = "sina --lca-quorum 0.5 -i $fa -o $out --ptdb $arb --search --search-db $arb --lca-fields $lca_field --meta-fmt header --intype fasta";
			#--lca-quorum 0.5
		#my $cmd_sina = "sina -i $fa -o $out --ptdb $arb --search --search-db $arb --lca-fields $lca_field --meta-fmt header --intype fasta";
		#my $cmd_sina = "sina -i $fa -o $out --ptdb $arb --search --search-db $arb --lca-fields $lca_field --meta-fmt header --intype fasta > $fa.log.sina.txt 2>&1";
			#Fatal error: Unable to open file "/dev/stderr" for writing.
		#my $cmd_sina = "sina -i $fa -o $out --ptdb $arb --search --search-db $arb --lca-fields $lca_field --meta-fmt header --intype fasta 2>&1 | tee $fa.log.sina.txt";
		#my $cmd_sina = "sina -i $fa -o $out --ptdb $arb --search --search-db $arb --lca-fields $lca_field --meta-fmt csv --intype fasta | tee $fa.log.sina.txt";
		#sina -i $f -o all_prinseq_good_igMu.fa.160607_0206ssuARB_HNG2_MK3.fasta --ptdb $arb --search --search-db $arb --lca-fields tax_yh --meta-fmt csv --intype fasta | tee log.arb.morikawa.txt
		#--ptport :/tmp/sina_pt_68329 --search-port :/tmp/sina_pt_68329
		
		warn "$cmd_sina\n";
		system "$cmd_sina";
	}
	
	my $f_tax = $out . '.tax.txt';
	my $f_log_sina = $out . '.log_sina.txt';	#this is parsed log of sina output fasta, not log file of sina itself. 
	if (
		-e $f_tax and -s $f_tax and -e $f_log_sina and -s $f_log_sina
			and
		(-f $out and -s $out and (stat $out)[9] <= (stat $f_tax)[9])
	){
		warn "tax list file '$f_tax' and sina parsed log file '$f_log_sina' exists, skip. \n";
	}else{
		$obj->parse_fa_sina();
	}
	
	my $f_summary = $out . '.tax.summary.txt';
	if (
		-e $f_summary and -s $f_summary
			and
		(-f $f_tax and -s $f_tax and (stat $f_tax)[9] <= (stat $f_summary)[9])
	){
		warn "tax summary file '$f_summary' exists, skip. \n";
	}else{
		$obj->summary_tax();
	}
	
	my $fa_arb = $out . '.arb.fa';
	if (-e $fa_arb and -s $fa_arb){
		warn "arb-compatible output file '$fa_arb' exists, skip. \n";
	}else{
		$obj->sina_fas2arb_fas();
	}
	
	return;
}
sub fa2sina_ver_sh_orig{
		#$obj->fa2sina_ver_sh();
	my $obj = shift;
	
	my $fa = $obj->{fa};
	my $arb = $obj->{arb};
	my $lca_field = $obj->{lca_field};
	my $out = $obj->{out};
	#my $is_sum_size = $obj->{is_sum_size};
	my $is_name_till_semicolon = $obj->{is_name_till_semicolon};
	
	mkpath(dirname($out));
	
	my $search_min_sim = 0.80;	#for sina
	my $idents_check = [80, 90, 95, 75];	#for summary
		#must be ascending order. 
	my $main_threshold_ident = $idents_check->[0];	#80; for summary
	#my $main_threshold_ident = 80;	#for summary
	
	if (-e $out and -s $out){
		warn "sina output file '$out' exists, skip. \n";
	}else{
		warn "sina output file '$out' not exists, do sina.. \n";
		my $logfile_sina = $out . '.sina.log.txt';
			#  --log-file arg (=/dev/stderr) log to file (default: STDERR)
			#to avoid error:
			#Fatal error: Unable to open file "/dev/stderr" for writing.
		my $cmd_sina = "sina --search-min-sim $search_min_sim -i $fa -o $out --ptdb $arb --search --search-db $arb --lca-fields $lca_field --meta-fmt header --intype fasta --log-file $logfile_sina";
		#my $cmd_sina = "sina --search-min-sim $search_min_sim -i $fa -o $out --ptdb $arb --search --search-db $arb --lca-fields $lca_field --meta-fmt header --intype fasta";
		#my $cmd_sina = "sina --search-min-sim 0.80 -i $fa -o $out --ptdb $arb --search --search-db $arb --lca-fields $lca_field --meta-fmt header --intype fasta";
		#my $cmd_sina = "sina --search-min-sim 0.90 -i $fa -o $out --ptdb $arb --search --search-db $arb --lca-fields $lca_field --meta-fmt header --intype fasta";
			#--search-min-sim arg
		#my $cmd_sina = "sina --search-all --search-no-fast -i $fa -o $out --ptdb $arb --search --search-db $arb --lca-fields $lca_field --meta-fmt header --intype fasta";
			#--search-all
			#--search-no-fast
		#my $cmd_sina = "sina --lca-quorum 0.5 -i $fa -o $out --ptdb $arb --search --search-db $arb --lca-fields $lca_field --meta-fmt header --intype fasta";
			#--lca-quorum 0.5
		#my $cmd_sina = "sina -i $fa -o $out --ptdb $arb --search --search-db $arb --lca-fields $lca_field --meta-fmt header --intype fasta";
		#my $cmd_sina = "sina -i $fa -o $out --ptdb $arb --search --search-db $arb --lca-fields $lca_field --meta-fmt header --intype fasta > $fa.log.sina.txt 2>&1";
			#Fatal error: Unable to open file "/dev/stderr" for writing.
		#my $cmd_sina = "sina -i $fa -o $out --ptdb $arb --search --search-db $arb --lca-fields $lca_field --meta-fmt header --intype fasta 2>&1 | tee $fa.log.sina.txt";
		#my $cmd_sina = "sina -i $fa -o $out --ptdb $arb --search --search-db $arb --lca-fields $lca_field --meta-fmt csv --intype fasta | tee $fa.log.sina.txt";
		#sina -i $f -o all_prinseq_good_igMu.fa.160607_0206ssuARB_HNG2_MK3.fasta --ptdb $arb --search --search-db $arb --lca-fields tax_yh --meta-fmt csv --intype fasta | tee log.arb.morikawa.txt
		#--ptport :/tmp/sina_pt_68329 --search-port :/tmp/sina_pt_68329
		
		warn "$cmd_sina\n";
		system "$cmd_sina";
	}
	
	my $f_tax = $out . '.tax.txt';
	my $f_log_sina = $out . '.log_sina.txt';	#this is parsed log of sina output fasta, not log file of sina itself. 
	if (
		-e $f_tax and -s $f_tax and -e $f_log_sina and -s $f_log_sina
			and
		(-f $out and -s $out and (stat $out)[9] <= (stat $f_tax)[9])
	){
		warn "tax list file '$f_tax' and sina parsed log file '$f_log_sina' exists, skip. \n";
	}else{
		open (IN, $out) || die "cannot open $out: $!";
		open (OUT_TAX, ">$f_tax") or confess "cannot open $f_tax: $!";
		open (OUT_LOG, ">$f_log_sina") or confess "cannot open $f_log_sina: $!";
		#my $i = -1;
		
		my $ks_tax = [qw(name lca_field tax align_ident_slv low_confident)];
		my $str_k_tax = join "\t", @$ks_tax;
		print OUT_TAX "$str_k_tax\n";
		
		my $ks_log = [qw(name)];
		my $seen_ks = {};
		my $is_to_out_ks_log = 1;
		while (<IN>){
			if (/^>(.+)$/){
				chomp;
				#$i++;
				my $fasta_name = $1;
				
				my $e = {
					name => $fasta_name,
				};
				$e->{name_full} = $e->{name};
				
				#my $pairs = [
				my $pairs = [];
				push @$pairs, $1, $2 while (
					$e->{name_full} =~ /
						\[
							([^=]+?)[=](
								(?&WORD) | (?&BRACKETED)
								(?(DEFINE)
									(?<WORD>      \s* [^\[\]]+ )
									(?<BRACKETED> \s* \[ (?&TEXT)? \s* \] )
									(?<TEXT>      (?: (?&WORD) | (?&BRACKETED) )+ )
								)
							)
						\]
					/xg
				);
				#];
					#[nearest_slv=ARB_25BABFAF...~0.951049 X89043.[EBI] X89043.1..~0.941725 ARB_C313EA80...~0.939394 ARB_14293161...~0.937063 ARB_8FADD83B...~0.934732 ARB_1DEA9552...~0.934732 ARB_C41E8E79...~0.927739 X79548.[EBI] X79548.1..~0.925408 ARB_72F20AB7...~0.925408 ARB_B2DE4B46...~0.923077 ] [nuc=429] [nuc_gene_slv=429] [turn_slv=none]
				#my $pairs = [$e->{name_full} =~ /\[([^=]+?)[=](.*?)\]/g];
				#print STDERR Dumper $pairs;die;#;stop;

#push @$pairs, $1
#	while /
#		\G \s*+ ( (?&WORD) | (?&BRACKETED) )
#		
#		(?(DEFINE)
#			(?<WORD>      \s* \w+ )
#			(?<BRACKETED> \s* \[ (?&TEXT)? \s* \] )
#			(?<TEXT>      (?: (?&WORD) | (?&BRACKETED) )+ )
#		)
#	/xg;

				my $h = {};
				for (my $i = 0; $i <= $#$pairs; $i += 2){
					my ($k, $v) = ($pairs->[$i], $pairs->[$i+1]);
					$h->{$k} = $v;
					push @$ks_log, $k if !$seen_ks->{$k}++;
				}
				#$h = {$e->{name_full} =~ /\[([^=]+?)[=](.+?)\]/g};
				
				$e->{name} =~ s/ .+//;
				$e->{name} =~ s/;.+// if $is_name_till_semicolon;	#Hak16S1_1;size=9279; -> Hak16S1_1
				$e->{lca_field} = $lca_field;
				$e->{tax} = $h->{'lca_' . $lca_field};	#lca_tax_slv, lca_tax_yh, ...
				$e->{align_ident_slv} = $h->{align_ident_slv};
				if ($e->{align_ident_slv} < 95 and $e->{align_ident_slv} >= 90){
					$e->{low_confident} = '-';
				}
				elsif ($e->{align_ident_slv} < 90 and $e->{align_ident_slv} >= 80){
					$e->{low_confident} = '--';
				}
				elsif ($e->{align_ident_slv} < 80){
					$e->{low_confident} = '---';
				}
				#$e->{align_startpos_slv} = $h->{align_startpos_slv};
				#$e->{align_stoppos_slv} = $h->{align_stoppos_slv};
				#[align_startpos_slv=988] [align_stoppos_slv=43283]
					#already recorded in log file as tab format.
				
				$h->{name} = $e->{name};
				#push @$es_log, $h;
				#print STDERR Dumper $e,$h,$ks_log,$ks_tax;die;#;stop;
				#my $f_tax = $out . '.tax.txt';
				foreach my $k (@$ks_tax){
					defined $e->{$k} or $e->{$k} = '';
				}
				my $str_tax = join "\t", @{$e}{@$ks_tax};
				print OUT_TAX "$str_tax\n";
				
				#my $f_log_sina = $out . '.log_sina.txt';
				foreach my $k (@$ks_log){
					defined $h->{$k} or $h->{$k} = '';
				}
				if ($is_to_out_ks_log){
					my $str_k_log = join "\t", @$ks_log;
					print OUT_LOG "$str_k_log\n";
					$is_to_out_ks_log = 0;
				}
				my $str_log = join "\t", @{$h}{@$ks_log};
				print OUT_LOG "$str_log\n";
			}
			else{
			}
		}
		close (OUT_TAX) or confess "cannot close $f_tax: $!";
		close (OUT_LOG) or confess "cannot close $f_log_sina: $!";
		close (IN) || die "cannot close $out: $!";
	}
	
	my $f_summary = $out . '.tax.summary.txt';
	if (
		-e $f_summary and -s $f_summary
			and
		(-f $f_tax and -s $f_tax and (stat $f_tax)[9] <= (stat $f_summary)[9])
	){
	#if (0 and -e $f_summary and -s $f_summary){
		warn "tax summary file '$f_summary' exists, skip. \n";
	}else{
		#my $idents_check = [qw(80, 90, 95)];
		my $idents = [];
		my $is_same_found = 0;
		my $ident0 = int($search_min_sim * 100);
		foreach my $ident (@$idents_check){
			if ($ident == $ident0){
				$is_same_found = 1;
				#last;
			}
			#if ($ident >= $ident0){	#in a few cases, there are hits with ident lower than search_min_sim are detected. 
				push @$idents, $ident;
			#}
		}
		if (!$is_same_found){
			push @$idents, int($search_min_sim * 100);
		}
		
		open (IN, $f_tax) || die "cannot open $f_tax: $!";
		my $ks = [];
		my $is_first = 1;
		my $tax2n = {};
		while (<IN>){
			next if !$is_first and /^#/;
			
			chomp;
			#warn "$_\n";
			if ($is_first){
				$ks = [split /\t/, $_];
				$is_first = 0;
			}else{
				my $els = [split /\t/, $_];
				my $e = {};
				@{$e}{@$ks} = @$els;
				
				foreach my $threshold_ident (@$idents){
					#foreach my $e (@$es){
					if ($e->{align_ident_slv} >= $threshold_ident){
						#if ($is_sum_size){
						
						#}else{
							$tax2n->{$threshold_ident}{$e->{tax}}++;
						#}
					}
						#$tax2n->{$e->{tax}}++ if $e->{align_ident_slv} >= $threshold_ident;
					#}
				}
				#print STDERR Dumper $e,$tax2n;die;#;stop;
			}
		}
		close (IN) || die "cannot close $f_tax: $!";
		#print STDERR Dumper $idents;die;#;stop;
		my $ks_count = [qw(tax count)];
		foreach my $threshold_ident (@$idents){
			warn "$threshold_ident..\n";
			my $es_s = [];
			foreach my $k (sort keys %{$tax2n->{$threshold_ident}}){
				push @$es_s, {tax => $k, count => $tax2n->{$threshold_ident}{$k}};
			}
			my $f_summary_w_ident = $out . '.tax.summary' . ".ident.$threshold_ident.txt";
			#my $es_s1 = [grep {$_->{align_ident_slv} >= $threshold_ident} @$es_s];
			output_hash_my($f_summary_w_ident, $es_s, K => $ks_count);
			#output_hash_my($f_summary_w_ident, $es_s1, K => $ks);
		}
		
		#my $main_threshold_ident = 80;
		##my $main_threshold_ident = 95;
		my $f_summary_w_ident = $out . '.tax.summary' . ".ident.$main_threshold_ident.txt";
		#my $f_summary = $out . '.tax.summary.txt';
		system "ln -s -f $f_summary_w_ident $f_summary";
		#output_hash_my(\*STDERR, $es_s, K => $ks);
	}
	
	my ($ks, $es_s) = get_hash_my($f_summary);
		#my $ks = [qw(tax count)];
	output_hash_my(\*STDERR, $es_s, K => $ks);
	
	my $fa_arb = $out . '.arb.fa';
	if (-e $fa_arb and -s $fa_arb){
		warn "arb-compatible output file '$fa_arb' exists, skip. \n";
	}else{
		my $cmd_sina2arb = "sina_fas2arb_fas.pl fa=$out out=$fa_arb";
		warn "\n", "when you load the sequences on ARB as pre-aligned, use following file.\n",
			"(gap characters are converted to '.' from '-', to be compatible with aligned ARB, at 5' and 3' ends.)",
			"$cmd_sina2arb\n";
		system "$cmd_sina2arb";
	}
	
	return;
}
sub _____________________Init_obj___________________{
}
sub init_parameters_default{
		#($args) = $obj->init_parameters_default(%$args);	#in 
		#($args) = $this->init_parameters_default(%$args);	#in 
			#delete $args->{$_} foreach (qw());
	
	my $usage = << 'USAGE';
	my $class = ref $obj;
	(my $class_str = lc $class) =~ s/::/__/g;
	#warn "log_entries_subset for clustid_str '$clustid_str'..\n";
	$obj->init_parameters_default(
		"set_obj_key_list__$class_str" => [qw()], 
		"set_obj_key_list_opt__$class_str" => [qw(dirlist list yymm)], 
	);	#in 
USAGE
	
	my $this = shift;
	my $args = {
		@_,
	};
	
	my $class = ref $this;
	(my $class_str = lc $class) =~ s/::/__/g;
	$this->{class_str} or $this->{class_str} = $class_str;
	
	$class_str = $this->{class_str};
	#my $class_str = $this->{class_str};
	
	my ($set_args_key_list_key, $set_args_key_list_opt_key, $set_obj_key_list_key, $set_obj_key_list_opt_key);
#	if ($args->{set_args_obj_key_list_type}){
#		$set_args_key_list_key = join '_', 'set', $args->{set_args_obj_key_list_type}, 'args_key_list';
#		$set_obj_key_list_key = join '_', 'set', $args->{set_args_obj_key_list_type}, 'obj_key_list';
#		$set_obj_key_list_opt_key = join '_', 'set', $args->{set_args_obj_key_list_type}, 'obj_key_list_opt';
#		($inits_key, $inits_str_key) = ('inits' . '_' . , 'inits_str');
#	}else{
		$set_args_key_list_key = 'set_args_key_list' . '__' . $class_str;
		$set_args_key_list_opt_key = 'set_args_key_list_opt' . '__' . $class_str;
		$set_obj_key_list_key = 'set_obj_key_list' . '__' . $class_str;
		$set_obj_key_list_opt_key = 'set_obj_key_list_opt' . '__' . $class_str;
#	}
#print STDERR Dumper $args;stop;
	if ($args->{"${set_args_key_list_key}__str"}){
		$args->{$set_args_key_list_key} = [split /,/, $args->{"${set_args_key_list_key}__str"}];
	}
	if ($args->{"${set_args_key_list_opt_key}__str"}){
		$args->{$set_args_key_list_opt_key} = [split /,/, $args->{"${set_args_key_list_opt_key}__str"}];
	}
	if ($args->{"${set_obj_key_list_key}__str"}){
		$args->{$set_obj_key_list_key} = [split /,/, $args->{"${set_obj_key_list_key}__str"}];
	}
	if ($args->{"${set_obj_key_list_opt_key}__str"}){
		$args->{$set_obj_key_list_opt_key} = [split /,/, $args->{"${set_obj_key_list_opt_key}__str"}];
	}
	
	if ($args->{$set_args_key_list_key}){
		foreach my $k (@{$args->{$set_args_key_list_key}}){
			$args->{$k} = defined $args->{$k} ? $args->{$k} : set($0, uc $k);
		}
	}
	if ($args->{$set_args_key_list_opt_key}){
		foreach my $k (@{$args->{$set_args_key_list_opt_key}}){
			#defined $args->{$k} or $this->{$k} = $args->{$k};
			$args->{$k} = defined $args->{$k} ? $args->{$k} : undef;#set($0, uc $k);
		}
	}
	if ($args->{$set_obj_key_list_key}){
		#use Data::Dumper;print STDERR Dumper $this;die;
		foreach my $k (@{$args->{$set_obj_key_list_key}}){
			if (!defined $this->{$k}){
				$this->{$k} = (defined $this->{cgi} and defined $this->{cgi}->param($k))
					? $this->{cgi}->param($k)
					: (defined $args->{$k})
					? $args->{$k}
					: (warn "$k=<$k> required.\n" && exit(1));	#set($0, uc $k);
					#: set($0, uc $k);
			}
			#!defined $this->{$k} and $this->{$k} = defined $args->{$k} ? $args->{$k} : set($0, uc $k);
		}
	}
	if ($args->{$set_obj_key_list_opt_key}){
		foreach my $k (@{$args->{$set_obj_key_list_opt_key}}){
			if (!defined $this->{$k}){
				$this->{$k} = (defined $this->{cgi} and defined $this->{cgi}->param($k))
					? $this->{cgi}->param($k)
					: (defined $args->{$k}) 
					? $args->{$k}
					: undef;
					#: set($0, uc $k);
			}
			#defined $this->{$k} or $this->{$k} = $args->{$k};
		}
	}
	
	return ($args);
}

sub ___________________Non_object_oriented___________________{
}

sub _____________________Read_write___________________{
}
sub get_hash_my{
		#my $es = get_hash_my($f);
		#my ($ks, $es) = get_hash_my($f);
	my $f_list = shift or confess "file not given.\n";
	my $ks = shift || undef;
	
	open (IN, $f_list) or confess "cannot open file '$f_list': $!";
	my $es_set = [];
	my $is_first;
	if ($ks and scalar @$ks){
		$is_first = 0;
	}else{
		$ks = [];
		$is_first = 1;
	}
	#my $ks = [];
	#my $is_first = 1;
	#my $is_first = 0;
	while (<IN>){
		next if !$is_first and /^#/;
		
		chomp;
		#warn "$_\n";
		if ($is_first){
			$ks = [split /\t/, $_];
			$is_first = 0;
		}else{
			my $els = [split /\t/, $_];
			my $e = {};
			@{$e}{@$ks} = @$els;
			push @$es_set, $e;
			#print @{$e}{@$ks}, "\n";
		}
	}
	close (IN) or confess "cannot close file '$f_list': $!";
	#print STDERR Dumper $es_set;#;stop;
	
	return ($ks, $es_set);
}
sub _output_hash_my{
		#output_hash_my($f_out, $es, K => $ks);
	my $OUT = shift;
	#my $f_out = shift or confess "file not given.\n";
	my $es = shift;
	my $args = {
		#K => [], 
		@_, 	#K
	};
	
	my $ks = $args->{K} || confess "key list not given.\n";
	#output_hash($o_f, $es_all, K => [@$ks_all, @$bases, @$bases_len], R => 1);
	#my $ks = [@$ks_all, @$bases, @$bases_len];
	
	#open (OUT, ">$f_out") or confess "cannot open file '$f_out': $!";
	my $str_k = join "\t", @$ks;
	print $OUT "$str_k\n";
	foreach my $e (@$es){
		foreach my $k (@$ks){
			defined $e->{$k} or $e->{$k} = '';
		}
		#print STDERR Dumper $e;
		my $str = join "\t", @{$e}{@$ks};
		print $OUT "$str\n";
		#}
	}
	#close (OUT) or confess "cannot close file '$f_out': $!";
	
	return;
}
sub output_hash_my{
		#output_hash_my($f_out, $es, K => $ks);
	my $outputfile_or_OUT = shift;
	#my $f_out = shift or confess "file not given.\n";
	my $es = shift;
	my $args = {
		#K => [], 
		@_, 	#K
	};
	
	my $ks = $args->{K} || confess "key list not given.\n";
	
	if (ref $outputfile_or_OUT eq ''){
		open (OUT, ">$outputfile_or_OUT") or confess "cannot open $outputfile_or_OUT: $!";
		my $OUT = \*OUT;
		_output_hash_my($OUT, $es, %$args);
		close (OUT) or confess "cannot close $outputfile_or_OUT: $!";
	}
	elsif (ref $outputfile_or_OUT eq 'GLOB'){
		my $OUT = $outputfile_or_OUT;
		_output_hash_my($OUT, $es, %$args);
	}
	
	return;
}
sub _____________________Date_time___________________{
}
sub _date_time{
		#my $date_time_table = _date_time();	#in package Operate_dir::Operate_dir_ver3;
		#my $date_time_table = _date_time($localtime);	#in package Operate_dir::Operate_dir_ver3;
	#my $hm_s_flag = shift || 3;
	my $localtime = shift;
		#(0 or '') -> '', 1 -> hm, 2 -> hms,
		#3 -> ymd, 4 -> ymd_hm, 5 -> ymd_hms, 
		#6 -> md
			#before ver18: default=0 for ymd., 0: ymd, 1: ymd_hm, 2: ymd_hms.
			#before ver15: only 0 -> ymdhm, others ('' or 1) -> ymd
	
	my ($sec, $min, $hour, $day, $mon_minus_1, $year_from_1900, $wday);
	if ($localtime){
		($sec, $min, $hour, $day, $mon_minus_1, $year_from_1900, $wday) = localtime($localtime);
	}else{
		($sec, $min, $hour, $day, $mon_minus_1, $year_from_1900, $wday) = localtime;
	}
	
	my $two_year = $year_from_1900 - 100;
	my $four_year = $year_from_1900 + 1900;
	my $mon = $mon_minus_1 + 1;
	
	foreach ($two_year, $mon, $day, $hour, $min, $sec){
		s/(\d+)/sprintf ("%02d", $1)/e ;
			#show "03" not "3".
			#2014 will be 14, not 114
	}
	
	my $ymd = $two_year . $mon . $day;
	my $hm = $hour . $min;
	
	my $table = {
		0 => '', 
		1 => $hm, 
		2 => $hm . $sec, 
		3 => $ymd, 
		4 => $ymd . '_' . $hm, 
		5 => $ymd . '_' . $hm . $sec, 
		6 => $mon . $day, 
		
		ymd => $ymd,
		ymd_hm => $ymd . '_' . $hm,
		ymd_hms => $ymd . '_' . $hm . $sec, 
		ym => $two_year . $mon, 
		y4m => $four_year . $mon, 
		md => $mon . $day, 
		
		hm => $hm, 
		hms => $hm . $sec, 
		
		year => $four_year, 
		four_year => $four_year, 
		two_year => $two_year, 
		mon => $mon, 
		day => $day, 
		hour => $hour, 
		min => $min, 
		sec => $sec,
	};
	
	return $table;
}
sub date_time{
		#my $date_time = date_time(3);	#in package Operate_dir::Operate_dir_ver3;
	my $hm_s_flag = shift || 3;
	my $localtime = shift;
	
	my $date_time_table = _date_time($localtime);
	
	my $date_time = $date_time_table->{$hm_s_flag};
	
	return $date_time;
}
sub _____________________Multifasta___________________{
}
sub get_mfa{
		#get contents of file. accept multifasta format.
		#can take file or fasta_ref.
		#usage:	$multifasta_set = get_contents_of_fasta($file_or_fasta_ref).
		#return:	ref to array of each fasta set, which is ref to hash which keys are {fasta_name} and {seq}.
	my $file_or_fasta_ref = shift;
	my $args = {
		as_is_flag => 0, 
		@_, 
	};
	
	my $multifasta_set;
	
	if (ref $file_or_fasta_ref eq 'SCALAR') {	#$fasta_ref.
		while ($$file_or_fasta_ref =~ /^>(.+)\n((?:^[^>\n].*\n)+)/mg){
			my $fasta_set = {
				name => $1, 
				fasta_name => $1,
				seq => $2,
			};
			$fasta_set->{seq} =~ s/\n//g if !$args->{as_is_flag};
			push @$multifasta_set, $fasta_set;
		}
	}
	elsif (ref $file_or_fasta_ref eq '') {	#$seqfile.
		my $separator = $/;	#if changed outside, $/ is global so now may not be "\n".
		$/ = "\n";
		
		open (IN, $file_or_fasta_ref) || die "cannot open $file_or_fasta_ref: $!";
		my $i = -1;
		while (<IN>){
			if (/^>(.+)$/){
				chomp;
				$i++;
				my $fasta_name = $1;
				$multifasta_set->[$i]->{fasta_name} = $fasta_name;
				$multifasta_set->[$i]->{name} = $fasta_name;
			}
			else{
				chomp if !$args->{as_is_flag};
	#			tr/[A-Z]/[a-z]/;
	#			my $illegal_num = s/[^acgtryswmkbdhvnu]//g;
	#			$illegal_num and print "\n$illegal_num illegal letters (except return keys) were removed from seq'$multifasta_set->[$i]->{fasta_name}'.\n\n";
	#			/[^acgtu]/ and print "Warn: Seq'$multifasta_set->[$i]->{fasta_name}' contains not 'a', 'c', 'g', 't', 'u' nucleotide(s).\n",
	#									"Is it surely nucleotide seq?\n\n";
				next unless $multifasta_set;
				$multifasta_set->[$i]->{seq} .= $_;
				#$multifasta_set->[$i]->{seq} .= $/ if $args->{as_is_flag} and $/ ne "\n";
			}
		}
		close (IN) || die "cannot close $file_or_fasta_ref: $!";
		$/ = $separator;
	}
	$multifasta_set;
}
sub output_mfa{
	#use Carp;
	
	my $outputfile_or_OUT = shift;
	my $fasta_set = shift;
	my $seq_type_or_options = shift || 'seq';
		#'nuc_seq' or 'aa_seq', or '0' or '1' for renew_flag, or hash_ref {} of options.
	my $renew_flag_given = shift;
	
	my $renew_flag;
	my $seq_type;
	my $name_key = 'name';
	my $ret_as_is_flag = 0;
	my $strict_format_flag = 0;	#if 1, add '//' terminator line. 
	my $drop_entries_w_no_seq_flag = 0;
	if (ref $seq_type_or_options eq ''){
		if ($seq_type_or_options =~ /^0$|^1$/){
			$seq_type = 'seq';
			$renew_flag = $seq_type_or_options;
		}else{
			$seq_type = $seq_type_or_options;
			$renew_flag = $renew_flag_given || 0;	#default is append, not renew.
		}
	}
	elsif (ref $seq_type_or_options eq 'HASH'){	#hash ref of options.
		$renew_flag = $seq_type_or_options->{renew_flag} || 0;	#default is append, not renew.
		$seq_type = $seq_type_or_options->{seq_type} || 'seq';
		$name_key = $seq_type_or_options->{name_key} || 'name';
		$ret_as_is_flag = $seq_type_or_options->{ret_as_is_flag} || 0;	#0, double return; 1, as is (as given). 
		$strict_format_flag = $seq_type_or_options->{strict_format_flag} || 0;
		$drop_entries_w_no_seq_flag = $seq_type_or_options->{drop_entries_w_no_seq_flag} || 0;
	}
	
	my $ret;
	if ($strict_format_flag){
		$ret_as_is_flag = 0;
		$ret = "//\n";
	}else{
		$ret = "\n";	#only used when $ret_as_is_flag as 0
	}
	#warn "seq_type: $seq_type\n";
	
	if (ref $outputfile_or_OUT eq ''){
		if ($renew_flag){
			open (OUT, ">$outputfile_or_OUT") or confess "cannot open $outputfile_or_OUT: $!";
		}
		else{
			open (OUT, ">>$outputfile_or_OUT") or confess "cannot open $outputfile_or_OUT: $!";
		}
		if ($ret_as_is_flag){
			foreach my $set (@$fasta_set){
				my $name = $set->{$name_key} || $set->{fasta_name};
				
				if (! ($set->{$seq_type} || $set->{nuc_seq})){
					warn "no seq: '$name'\n";
					if (!$drop_entries_w_no_seq_flag){
						print OUT ">", $name, "\n";
						print OUT $set->{$seq_type} || $set->{nuc_seq};
					}
				}else{
					print OUT ">", $name, "\n";
					print OUT $set->{$seq_type} || $set->{nuc_seq};
				}
			}
		}else{
			foreach my $set (@$fasta_set){
				my $name = $set->{$name_key} || $set->{fasta_name};
				#print STDERR Dumper $set,$seq_type;stop;
				if (! ($set->{$seq_type} || $set->{nuc_seq})){
					warn "no seq: '$name'\n";
					if (!$drop_entries_w_no_seq_flag){
						print OUT ">", $name, "\n";
						print OUT $set->{$seq_type} || $set->{nuc_seq}, "\n$ret";
					}
				}else{
					print OUT ">", $name, "\n";
					print OUT $set->{$seq_type} || $set->{nuc_seq}, "\n$ret";
					#print OUT ">", $name, "\n", $set->{$seq_type} || $set->{nuc_seq}, "\n\n";
				}
			}
		}
		close (OUT) or confess "cannot close $outputfile_or_OUT: $!";
	}
	elsif (ref $outputfile_or_OUT eq 'GLOB' or ref $outputfile_or_OUT eq 'FileHandle'){
		if ($ret_as_is_flag){
			foreach my $set (@$fasta_set){
				print $outputfile_or_OUT ">", $set->{$name_key} || $set->{fasta_name}, "\n", 
					$set->{$seq_type} || $set->{nuc_seq};
			}
		}else{
			foreach my $set (@$fasta_set){
				print $outputfile_or_OUT ">", $set->{$name_key} || $set->{fasta_name}, "\n", 
					$set->{$seq_type} || $set->{nuc_seq}, "\n$ret";
					#$set->{$seq_type} || $set->{nuc_seq}, "\n\n";
			}
		}
	}
}

sub _____________________Seq___________________{
}
sub complement{
		#$sec_c = complement($seq);
		#$sec_c_ref = complement($seq_ref);
	my $seq = shift;
	
	my $seq_c;
	if (ref $seq eq 'SCALAR'){	#scalar_ref is given.
#		$$seq =~ tr/[A-Z]/[a-z]/;
		(my $seq_c0 = reverse $$seq) 
			=~ tr/[acgtryswmkbdhvnACGTRYSWMKBDHVN]/[tgcayrswkmvhdbnTGCAYRSWKMVHDBN]/;
		$seq_c = \$seq_c0;
	}
	elsif (ref $seq eq ''){	#scalar
#		$seq =~ tr/[A-Z]/[a-z]/;
		($seq_c = reverse $seq) 
			=~ tr/[acgtryswmkbdhvnACGTRYSWMKBDHVN]/[tgcayrswkmvhdbnTGCAYRSWKMVHDBN]/;
	}
	else {
		die "cannot transform ref type ", ref $seq, ".",
				" only 'SCALAR' and '' are allowed.:$!";
	}
	
	$seq_c;
}
sub ___________________Usage___________________{
}
sub printUsage {
	my $msg = shift;
	
	my $prog = basename($0);
	print STDERR << "MSG";
Requirement: 
	SINA should be in PATH variable.

Usage: 
	$prog fa=<fasta_file> arb=<arb_file> (lca=<0,1; default 1>) (lca_field=<lca field of arb; default tax_slv>) (out=<out; default fa.arb.fasta>)
	
Options:
	fa: multiple fasta file
	arb: arb file
	lca: do lca tax assinment of sina or not
		default: 1 (do lca search)
	lca_field: field of arb to be used for lca assignment
		default: tax_slv (for arb file of silva) can be self set field (e.g. tax_yh)
	out: output fasta file
		default: <fa>.<arb>.fasta
		<out>.log_sina.txt for lca assigned result (parsed from the sina record in the header of each entry in the output fasta file)
		<out>.tax.txt: tab-separated file for entry and lca-assigned tax
		<out>.tax.summary.txt: summary of lca-assigned tax result
		<out>.tax.summary.ident.80.txt: summary of lca-assigned tax result
		<out>.tax.summary.ident.90.txt: summary of lca-assigned tax result
		<out>.tax.summary.ident.95.txt: summary of lca-assigned tax result
		<out>.arb.fa: ARB-compatible as pre-aligned sequences (gap characters at 5' and 3' ends are '.' insted of '-'. 
Note:
	If you are just executing multiple fa2sina.pl simulutaneiously, check the arb file already has pt server. 
MSG
	# Note for developers:
	#  Be careful not to use this script itself as a command called by systemcall in other script that use CGI.pm. 
	#  This script use CGI.pm to parse options. 
	#  Nested usage of CGI.pm does not work in general. 
	#  To use CGI.pm in nested way, you need to call the systemcall by the way below. 
	#   systemcall "QUERY_STRING='';QUERY_STRING='param1=blabla&param2=hogehoge';cmd_to_call.cgi";
	
	warn "###\n", $msg, "\n" if $msg;
	
	exit(0);
}
sub ___________________Main___________________{
}
if ($0 eq __FILE__) {
	my ($cgi) = CGI->new();
	
	my $pwd = cwd();
	
	my $method_str = defined $cgi->param('methods') ? $cgi->param('methods') : 'fa2sina_ver_sh';
	my $fa = defined $cgi->param('fa') ? $cgi->param('fa') : printUsage("fa=<fasta_file> is required.");
	my $arb = defined $cgi->param('arb') ? $cgi->param('arb') : printUsage("arb=<arb_file> is required.");
	my $lca = defined $cgi->param('lca') ? $cgi->param('lca') : 1;
	my $lca_field = defined $cgi->param('lca_field') ? $cgi->param('lca_field') : 'tax_slv';	#tax_yh for hongoh-arb
	my $out = defined $cgi->param('out') ? $cgi->param('out') : $fa . '.' . basename($arb) . '.fasta';
	#my $is_sum_size = defined $cgi->param('is_sum_size') ? $cgi->param('is_sum_size') : 1;
	my $is_name_till_semicolon = defined $cgi->param('is_name_till_semicolon') ? $cgi->param('is_name_till_semicolon') : 0;
	#my $is_name_till_semicolon = defined $cgi->param('is_name_till_semicolon') ? $cgi->param('is_name_till_semicolon') : 1;
		#trim usearch size string
		#e.g.
		#Hak16S1_1;size=9279; -> Hak16S1_1
	
	my $methods = [split /,/, $method_str];
	
	foreach my $cmd_name ('sina', 'sina_fas2arb_fas.pl'){
		chomp (my $cmd = `which $cmd_name`);
		if (!$cmd){
			printUsage("cmd '$cmd' not found. it must be in your PATH variable.");
		}
	}
	
	my $obj = {
		cgi => $cgi, 
		
		fa => $fa,
		
		arb => $arb, 
		lca => $lca, #0, sina wo lca, then sina_fas2arb_fas.pl
		lca_field => $lca_field,
		out => $out, 
		
		#is_sum_size => $is_sum_size, 
		is_name_till_semicolon => $is_name_till_semicolon, 
	};
	bless $obj, 'main';
	
#	my $method = "fa2sina";
#	$obj->$method();
	foreach my $method (@$methods){
		if ($obj->can($method)){
			warn "$method..\n";
			$obj->$method();
		}else{
			warn "method '$method' unavailable. typo?\n";
		}
	}
}
1;

# fa2sina

This script does (1) sina search against arb file, (2) parse the output, and (3) make summary reports.  
SINA program and ARB ribosomal database files are distributed by SILVA. 

## Requirement: 
	SINA should be in PATH variable.

## Usage: 
	fa2sina.pl fa=<fasta_file> arb=<arb_file> (lca=<0,1; default 1>) (lca_field=<lca field of arb; default tax_slv>) (out=<out; default fa.arb.fasta>)
	
## Options:
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

## Note:
	If you are just executing multiple fa2sina.pl simulutaneiously, check the arb file already has pt server. 
      

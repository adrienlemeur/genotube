class myFunctions {

	public static boolean checkFile(String path, String name, String extension) {

		def path_ = path.toString()
		def query = "find $path_/ -name $name*$extension"
		def answer = query.execute().text.count("\n") == 0 ? false : true

		return(answer)
	}

	public static boolean createSymLink(String source, String destination){
		def baseName = source.replaceFirst(~/\.[^\.]+$/, '')
		"ln -s $source $destination".execute()
		return(1)
	}

	public static boolean checkFORCE(String selected_param, String OPTION) {
		return selected_param in OPTION.tokenize(',')
	}

	public static void checkParameters(variant_caller, contam_check, species, taxonomy, treebuild, FORCE, SKIP, REMOVE, krakenDB){
		def quit = false;

		if( ! krakenDB in ['8G', '16G', 'none'] ){
			println "\n'$krakenDB' is not a valid argument for --download_K2_DB. Supported values are 'none' (default), '8G' and '16G'.\n"
			quit = true
		}

		if( ! variant_caller in ['gatk', 'samtools', 'freebayes'] ){
			println "\n'$variant_caller' is not a valid argument for --variant_caller. Supported values are 'gatk' (default), 'samtools' and 'freebayes'.\n"
			quit = true
		}

		if( ! contam_check in ['taxo_class', 'compet_mapping', 'none' ] ){
			println "\n'$contam_check' is not a valid argument for --contam_check. Supported values are 'taxo_class' (default), 'compet_mapping' and 'none'.\n"
			quit = true
		}

		if( ! treebuild in ['none', 'classic', 'fast', 'fastest' ] ){
			println "\n'$treebuild' is not a valid argument for --build_tree. Supported values are 'none' (default), 'classic', 'fast', 'fastest'.\n"
			quit = true
		}

		if( ! species in ['MTBC', 'other'] ) {
			println "\n'$species' is not a valid argument for --species. Supported values are 'MTBC' (default) and 'other'.\n"
			quit = true
		}

		if( ! taxonomy in ['barcode', 'placement', 'none'] ){
			println "\n'$taxonomy' is not a valid argument for --taxonomy. Supported values are 'barcode' (default), 'placement' or 'none'.\n"
			quit = true
		}

		if( FORCE == true ){ println "\nNo argument following FORCE, set to 'none'" ; FORCE = ['none'] ; } else { FORCE = FORCE?.tokenize(',') ; }

		def forceParameters = ['DOWNLOAD', 'TRIM', 'MAP', 'CALL', 'ANN', 'none'] as String[]
		if(! FORCE.every{it -> it in forceParameters}) {
			println "\n${FORCE.join(', ')} is not a valid argument for --FORCE. Should be a a comma separated list of 'DOWNLOAD', 'TRIM', 'MAP', 'CALL', 'ANN' and/or 'none' (default)."
			quit = true
		}

 		if( SKIP == true ){ println "\nNo argument following --SKIP, set to 'none'" ; SKIP = ['none'] ; } else { SKIP = SKIP?.tokenize(',') ; }

		def skipParameters = ['COVERAGE_CHECK', 'PROFILING', 'none'] as String[]
		if(SKIP != 'empty' && ! SKIP.every{it -> it in skipParameters}) {
			println "\n${SKIP.join(', ')} is not a valid argument for --SKIP. Should be a a comma separated list of 'COVERAGE_CHECK', 'PROFILING', and/or 'none' (default) ."
			quit = true
		}

		if( REMOVE == true ){ println "\nNo argument following --REMOVE, set to 'none'" ; REMOVE = ['none'] ; } else { REMOVE = REMOVE?.tokenize(',') ; }

		def removeParameters = ['FASTQ', 'BAM', 'VCF', 'none'] as String[]
		if(REMOVE != 'empty' && ! REMOVE.every{it -> it in removeParameters}) {
			println "\n${REMOVE.join(', ')} is not a valid argument for --REMOVE. Should be a a comma separated list of 'FASTQ', 'BAM', 'VCF' and/or 'none' (default) ."
			quit = true
		}

		if(quit){
			print "\n"
			System.exit(0)
		}

	}
	public static checkFolders(){} // TODO : one overloaded function to rule them all

}

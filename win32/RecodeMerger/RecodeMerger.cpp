// RecodeMerger.cpp : Defines the entry point for the console application.
//

#ifndef ENABLE_MULTIPLE_COMPRESSIONS
#define ENABLE_MULTIPLE_COMPRESSIONS 0
#endif

#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <assert.h>

#include "zlib.h"

#if ENABLE_MULTIPLE_COMPRESSIONS
#include "bzlib.h"
#include "lzma.h"
#include "snappy-c.h"
#include "lz4.h"
#endif

#include <tchar.h>
#include <wchar.h>
#include <SimpleOpt.h>

#ifdef _WIN32
#include <win32\dirent.h>
#include <windows.h>
#else
#include <dirent.h>
#endif

#include "recodefs.h"
#include "initparser.h"
#include "recode_utils.h"
#include "mrchandler.h"
#include "sequence_reader.h"
#include "commons.h"
#include "compressor.h"
#include "logger.h"
#include "fileutils.h"
/*
#include "mrc_header.h"
*/
#include "argparser.h"
#include "recode_header.h"
/*
#include "ReCoDeConfig.h"
*/
#include "L1.h"
#include "L2.h"
#include "L3.h"
#include "L4.h"

void mergePartFiles(InitParams *init_params, InputParams *params, RCHeader *header) {

	const char *compressed_filename = makeCompressedFilename(init_params->image_filename, params->reduction_level);
	char** part_filenames = (char**)malloc((params->num_threads) * sizeof(char*));
	int part_num;
	for (part_num = 0; part_num < params->num_threads; part_num++) {
		part_filenames[part_num] = makePartFilename(part_num, init_params->run_name, params->reduction_level);
	}
	merge_RC1_Parts(init_params->output_directory, part_filenames, params->num_threads, header, compressed_filename);
}

#ifdef _WIN32
int _tmain(int argc, TCHAR * argv[]) {

	omp_set_nested(1);

	InitParams *init_params = (InitParams *)malloc(sizeof(InitParams));
	parse_recode_server_init_params(argc, argv, &init_params);

	// read and validate input parameters
	InputParams *params = (InputParams *)malloc(sizeof(InputParams));
	get_input_params(init_params->params_filename, &params);

	// if MRC file, replace params from MRC header
	if (params->source_file_type == 1) {
		MRCHeader *mrc_header = (MRCHeader *)malloc(sizeof(MRCHeader));
		parseMRCHeader(init_params->image_filename, &mrc_header);
		compile_missing_params(mrc_header, &params);
		if (params->source_bit_depth != 16) {
			printf("ReCoDe currently only supports 16-bit unsigned data. The provided data has bit-depth %d.", params->source_bit_depth);
			exit(1);
		}
		free(mrc_header);
	}

	/*
	==========================================================================================
	Create ReCoDE header
	==========================================================================================
	*/
	RCHeader *header = (RCHeader *)malloc(sizeof(RCHeader));
	const char *image_name = getFilenameFromPath(init_params->image_filename);
	const char *dark_name = getFilenameFromPath(init_params->dark_filename);
	create_recode_header(params, -1, image_name, dark_name, &header);
	print_recode_header(header);

	mergePartFiles(init_params, params, header);

	return 1;
}
#else

#endif

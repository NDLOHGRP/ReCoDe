
/**************************************************************************************************************
to compile:	gcc recode.c -o ../bin/recode -I../include -O3 -lz -lbz2 -llzma -lsnappy -llz4 -lm -fopenmp
to run:		export OMP_NUM_THREADS=11
			./recode <input image file name> <dark image filename> <input parameter filename> <output directory>

-lz: links zlib
-lm: links math.h NOTE: it is important to keep -lm as the last option
-o:	 write build output to file
-O3: optimizer flag
**************************************************************************************************************/

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

#include "zlib.h"

#if ENABLE_MULTIPLE_COMPRESSIONS
#include "bzlib.h"
#include "lzma.h"
#include "snappy-c.h"
#include "lz4.h"
#endif

#ifdef _WIN32
	#include <win32\dirent.h>
#else
	#include <dirent.h>
#endif

#include "recodefs.h"
#include "recode_utils.h"
#include "commons.h"
#include "compressor.h"
#include "logger.h"
#include "fileutils.h"
/*
#include "mrc_header.h"
*/
#include "mrchandler.h"
#include "argparser.h"
#include "recode_header.h"

/*
#include "ReCoDeConfig.h"
*/

#include "L1.h"
#include "L2.h"
#include "L3.h"
#include "L4.h"




int de (int argc, char *argv[]) {
	
	/*
	========================================================================================== 
	Process Inputs
	==========================================================================================
	*/
	if (argc != 4) {
		printf("Failed to initialize. When mode is DE, expects 3 arguments, but found %d\n", argc-1);
		fprintf(stderr, "Usage::\n");
		fprintf(stderr, "./recode -de <recode filename> <output directory>\n");
		exit(0);
	}
	
	
	/*
	========================================================================================== 
	Decompress-Expand
	==========================================================================================
	*/
    printf("Decompressing and expanding %s ", argv[2]);
	uint16_t *de_frameBuffer16;
	RCHeader *de_header = (RCHeader *)malloc(sizeof(RCHeader));
	decompressExpand_L1 (argv[2], &de_frameBuffer16, &de_header);

	
	/*
	========================================================================================== 
	Save decoded-expanded image
	==========================================================================================
	*/
	char *image_name = getFilenameFromPath (argv[2]);
	char *s = get_filename_sans_extension (image_name);
    char *outDir = format_directory_path(argv[3]);
    char *out_fname = concat(outDir, concat("Recoded_", s));
    
    printf("to %s\n", out_fname);
	FILE *fp = fopen (out_fname, "w+");
	uint32_t n_pixels = de_header->nx * de_header->ny * de_header->nz;
	fwrite (de_frameBuffer16, sizeof(uint16_t), n_pixels, fp);
	fclose(fp);
	printf("Done.\n");
	
	/*
	========================================================================================== 
	Clean-up
	==========================================================================================
	*/
	free(s);
	free(de_header);
	free(de_frameBuffer16);
	
}


int rc (int argc, char *argv[]) {
	
	/*
	========================================================================================== 
	Process Inputs
	==========================================================================================
	*/
	if (argc != 6) {
		printf("Failed to initialize. When mode is RC, expects 5 arguments, but found %d\n", argc-1);
		fprintf(stderr, "Usage::\n");
		fprintf(stderr, "./recode -rc <input image filename> <dark image filename> <input parameter filename> <output directory>\n");
		exit(0);
	}
	
    char *darkFile = argv[3];
	char *imageFile = argv[2];
    char *outDir = format_directory_path(argv[5]);
    
	// read and validate input parameters
	InputParams *params = (InputParams *)malloc(sizeof(InputParams));
	get_input_params (argv[4], &params);
	
	// if MRC file replace missing params from MRC header
	if (params->source_file_type == 1) {
		MRCHeader *mrc_header = (MRCHeader *)malloc(sizeof(MRCHeader));
        parseMRCHeader (imageFile, &mrc_header);
		compile_missing_params (mrc_header, &params);
        free(mrc_header);
	}
	
	
	/*
	========================================================================================== 
	Create ReCoDE header
	==========================================================================================
	*/
	RCHeader *header = (RCHeader *)malloc(sizeof(RCHeader));
	uint8_t *image_name = (uint8_t*)getFilenameFromPath (imageFile);
	uint8_t *dark_name = (uint8_t*)getFilenameFromPath (darkFile);
	create_recode_header (params, -1, image_name, dark_name, &header);
	print_recode_header (header);

	
	/*
	========================================================================================== 
	Load image data
	==========================================================================================
	*/
	DataSize b = { header->nx, header->ny, params->num_frames, params->source_bit_depth };
	uint64_t sz_frameBuffer = b.nx * b.ny * b.nz * sizeof(uint16_t);
	uint16_t *frameBuffer = (uint16_t *)malloc(sz_frameBuffer);
	
	loadData (imageFile, frameBuffer, header->frame_offset, header->nz, b);
	recode_print("Image data loaded.\n");
	
	/*
	========================================================================================== 
	Dark Noise Estimation
	==========================================================================================
	*/
	uint32_t n_pixels_in_frame = header->nx * header->ny;
	uint16_t *darkFrame = (uint16_t *)calloc(n_pixels_in_frame, sizeof(uint16_t));
	
	recode_print("header->num_dark_frames = %lu\n", header->num_dark_frames);
	
	if (header->num_dark_frames > 1) {
	
		DataSize d = { header->nx, header->ny, params->num_dark_frames, params->source_bit_depth };
		uint64_t sz_darkBuffer = d.nx * d.ny * d.nz * sizeof(uint16_t);
		uint16_t *darkBuffer = (uint16_t *)malloc(sz_darkBuffer);
		
		loadData (darkFile, darkBuffer, 0, header->num_dark_frames, d);
		printf("Dark data loaded.\n");
		
		getDarkMax (darkBuffer, d, darkFrame, params->num_threads);
		free(darkBuffer);
		
		printf("Dark noise estimation complete.\n");
		
		
	} else {
		
		DataSize d1 = { header->nx, header->ny, 1, params->source_bit_depth };
		loadData (darkFile, darkFrame, 0, 1, d1);
		printf("Pre-computed dark noise estimates loaded.\n");
		
	}
	
	/*
	========================================================================================== 
	Allocate space to hold part file names returned by do_L1_Reduce_Compress
	Allocate space to hold the number of frames, frame start indices, and timing results for each thread
	==========================================================================================
	*/
	int part_num;
	char** part_filenames = (char**)malloc((params->num_threads)*sizeof(char*));
	uint32_t* n_frames_in_thread = (uint32_t*)malloc((params->num_threads)*sizeof(uint32_t));
	uint32_t* frame_start_indices = (uint32_t*)malloc((params->num_threads)*sizeof(uint32_t));
	
	int rem = params->num_frames % params->num_threads;
	int min_frames_per_thread = floor(((params->num_frames)*1.0)/((params->num_threads)*1.0));
	unsigned long frame_start_index = 0;
	
	for (part_num = 0; part_num < params->num_threads; part_num++) {
		part_filenames[part_num] = (char*)malloc(MAX_PART_FILE_NAME_LENGTH*sizeof(char));
		if (part_num < rem) {
			n_frames_in_thread[part_num] = min_frames_per_thread + 1;
		} else {
			n_frames_in_thread[part_num] = min_frames_per_thread;
		}
		
		frame_start_indices[part_num] = frame_start_index;
		recode_print("RCT %d: n_frames_in_thread = %lu, frame_start_index = %lu\n", part_num, n_frames_in_thread[part_num], frame_start_index);
		frame_start_index += n_frames_in_thread[part_num];
	}
	
	float *compress_times_sizes_best_speed = (float*) malloc (params->num_threads * 2 * sizeof(float));
	
	
	/*
	========================================================================================== 
	Reduce and Compress into part-files
	==========================================================================================
	*/
	unsigned long nz;
	
	double omp_start = omp_get_wtime();
	//struct timeval start, end;
	//gettimeofday(&start, NULL);

	double start = _gettimeofday();
    
    printf("Reducing compressing data with %d OpenMP thread(s).\n", params->num_threads);

	#pragma omp parallel for num_threads(params->num_threads)
	for (part_num = 0; part_num < params->num_threads; part_num++) {
		
		recode_print("RCT %d: NZ = %u\n", part_num, n_frames_in_thread[part_num]);
		
		DataSize part_b = { header->nx, header->ny, n_frames_in_thread[part_num], params->source_bit_depth };
		
		unsigned long chunk_offset = frame_start_indices[part_num] * part_b.nx * part_b.ny;
		
		recode_print("RCT %d: frame_start_index = %lu, Chunk Offset = %lu\n", part_num, frame_start_indices[part_num], chunk_offset);
		part_filenames[part_num] = reduceCompress_L1 ( 	part_num, 
														(const char*)image_name,
														outDir,
														frameBuffer + chunk_offset, 
														darkFrame, 
														frame_start_indices[part_num], 
														part_b, 
														params,
														compress_times_sizes_best_speed
													);		// best speed
													
		recode_print("RCT %d: Reduction Time = %f, Compression Time = %f, Total Time = %f \n", 
				part_num, 
				compress_times_sizes_best_speed[part_num*2 + 0], 
				compress_times_sizes_best_speed[part_num*2 + 1], 
				compress_times_sizes_best_speed[part_num*2 + 0] + compress_times_sizes_best_speed[part_num*2 + 1]);
		
		//logResult (concat("GZIP-S L1_d12 ", "0.7"), compress_times_sizes_best_speed);
	}
    printf("Done.\n");
	
	double total_time = omp_get_wtime()-omp_start;
	//gettimeofday(&end, NULL);
	double end = _gettimeofday();
	double delta = ((end  - start) * 1000000u + end - start) / 1.e6;
	printf("Wall Time: %f, Total OMP Time: %f, Effective Avg. Time per Frame: %f\n", delta, total_time, total_time/(part_num*1.0));
	
	/*
	========================================================================================== 
	Merge part-files
	==========================================================================================
	*/
	printf("Merging part files.\n");
	char* fname = merge_RC1_Parts(outDir, part_filenames, params->num_threads, header, concat(outDir, (const char*)image_name));
	printf("Done.\n");
	
	/*
	========================================================================================== 
	Clean-up
	==========================================================================================
	*/
	free(header);
    free(darkFrame);
	free(frameBuffer);
    free(n_frames_in_thread);
    free(frame_start_indices);
	for (part_num = 0; part_num < params->num_threads; part_num++) {
		free(part_filenames[part_num]);
	}
	free(params);
    free(compress_times_sizes_best_speed);
	free(part_filenames);

	return 0;
}

int main(int argc, char *argv[]) {

	if (argc < 2) {

		printf("Failed to initialize. Expected 3 or 5 arguments, but found %d\n", argc - 1);
		fprintf(stderr, "Usage::\n");
		fprintf(stderr, "./recode -rc <input image filename> <dark image filename> <input parameter filename> <output directory>\n");
		fprintf(stderr, "./recode -de <recode filename> <output directory>\n");
		exit(0);

	}
	else {

		if (strcmp_ignore_case(argv[1], "-rc")) {
			rc(argc, argv);
		}
		else if (strcmp_ignore_case(argv[1], "-de")) {
			de(argc, argv);
		}
		else {
			printf("Failed to initialize. Unknown mode: %s. Acceptable modes are: -rc and -de", argv[1]);
			fprintf(stderr, "Usage::\n");
			fprintf(stderr, "./recode -rc <input image filename> <dark image filename> <input parameter filename> <output directory>\n");
			fprintf(stderr, "./recode -de <recode filename> <output directory>\n");
			exit(0);
		}

	}
}
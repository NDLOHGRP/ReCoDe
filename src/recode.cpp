
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

#if defined(_MSC_VER)
# include <windows.h>
# include <tchar.h>
#else
# define TCHAR		char
# define _T(x)		x
# define _tprintf	printf
# define _tmain		main
#endif
#include <SimpleOpt.h>

#ifdef _WIN32
	#include <win32\dirent.h>
#else
	#include <dirent.h>
	#include <ctype.h>
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
#include "argparser.h"
#include "recode_header.h"

/*
#include "ReCoDeConfig.h"
*/

#include "L1.h"
#include "L2.h"
#include "L3.h"
#include "L4.h"

#include "recode_tests.c"

char* reduceCompress(
	uint8_t	    process_id,
	const char  *original_filename,
	const char  *out_foldername,
	uint16_t    *frameBuffer, 				// points to the start of this thread's chunk
	uint16_t    *darkFrame,
	uint32_t    frame_start_index,
	DataSize    h,
	InputParams *params,
	RCHeader    *header,
	float 		*compression_time
) {

	if (params->reduction_level == 1) {
		return reduceCompress_L1(process_id, original_filename, out_foldername, frameBuffer, darkFrame, frame_start_index, h, params, header, compression_time);
	}
	else if (params->reduction_level == 2) {
		return reduceCompress_L2(process_id, original_filename, out_foldername, frameBuffer, darkFrame, frame_start_index, h, params, compression_time);
	}
	else if (params->reduction_level == 3) {
		return reduceCompress_L3(process_id, original_filename, out_foldername, frameBuffer, darkFrame, frame_start_index, h, params, compression_time);
	}
	else if (params->reduction_level == 4) {
		return reduceCompress_L4(process_id, original_filename, out_foldername, frameBuffer, darkFrame, frame_start_index, h, params, header, compression_time);
	}
	else {
		printf("Error in function reduceCompress() in recode.cpp. ReCoDe has reduction levels 1 to 4. %d is not a supported level.\n", params->reduction_level);
	}

}

//void decompressExpand(const char* compressed_filename, uint16_t **frameBuffer, RCHeader **header) {
void decompressExpand(FILE* rc_fp, const char* out_filename, RCHeader *header) {

	//FILE* fp = fopen(compressed_filename, "rb");
	//parse_recode_header(fp, out_fp);

	if (header->reduction_level == 1) {
		if (header->recode_operation_mode == 0) {
			return decompressExpand_L1_Reduced_Only(rc_fp, out_filename, header);
		}
		else if (header->recode_operation_mode == 1) {
			//return decompressExpand_L1_Reduced_Compressed(rc_fp, out_filename, header);
			return decompressExpand_L1_Reduced_Compressed_Sparse(rc_fp, out_filename, header);
			//return _h_decompressExpand_L1_Reduced_Compressed_Sparse(rc_fp, out_filename, header);
		}
	}
	/*
	else if ((*header)->reduction_level == 2) {
		return decompressExpand_L2(fp, frameBuffer, header);
	}
	else if ((*header)->reduction_level == 3) {
		return decompressExpand_L3(fp, frameBuffer, header);
	}
	else if ((*header)->reduction_level == 4) {
		if ((*header)->recode_operation_mode == 0) {
			return decompressExpand_L4_Reduced_Only(fp, frameBuffer, header);
		}
		else if ((*header)->recode_operation_mode == 1) {
			//return decompressExpand_L4_Reduced_Compressed(fp, frameBuffer, header);
			return _decompressExpand_L4(fp, frameBuffer, header);
		}
	}
	*/
	else {
		printf("Error in function decompressExpand() in recode.cpp. ReCoDe has reduction levels 1 to 4. %d is not a supported level. This could indicate that the source file is not a ReCoDe file.\n", header->reduction_level);
	}
}

int de (InitParams *init_params) {
	
	/*
	========================================================================================== 
	Decompress-Expand
	==========================================================================================
	*/
    printf("Decompressing and expanding %s ", init_params->image_filename);
	RCHeader *de_header = (RCHeader *)malloc(sizeof(RCHeader));
	FILE* fp = fopen(init_params->image_filename, "rb");
	parse_recode_header(fp, &de_header);
	//uint16_t *de_frameBuffer16;
	//decompressExpand (init_params->image_filename, &de_frameBuffer16, &de_header);
	
	/*
	========================================================================================== 
	Create decoded-expanded image file
	==========================================================================================
	*/
	//char *image_name = getFilenameFromPath (init_params->image_filename);
	char *s = get_filename_sans_extension ((char*)&(de_header->source_file_name));
    char *outDir = format_directory_path(init_params->output_directory);
    char *out_fname = concat(outDir, concat("Recoded_", concat(s,".bin")));
    printf("to %s\n", out_fname);

	/*
	==========================================================================================
	Decompress-Expand
	==========================================================================================
	*/
	decompressExpand(fp, out_fname, de_header);

	//FILE *fp = fopen (out_fname, "wb+");
	//uint32_t n_pixels = de_header->nx * de_header->ny * de_header->nz;
	//fwrite (de_frameBuffer16, sizeof(uint16_t), n_pixels, fp);
	fclose(fp);
	printf("Done.\n");
	
	/*
	========================================================================================== 
	Clean-up
	==========================================================================================
	*/
	free(s);
	free(de_header);
	//free(de_frameBuffer16);

	return 0;
}

int merge (InitParams *init_params) {

	/*
	==========================================================================================
	Read and Validate Input Parameters
	==========================================================================================
	*/
	InputParams *params = (InputParams *)malloc(sizeof(InputParams));
	get_input_params(init_params->params_filename, &params);

	/*
	==========================================================================================
	Process Inputs
	==========================================================================================
	*/
	char *imageFile = init_params->image_filename;
	const char *image_name = getFilenameFromPath(imageFile);
	char *outDir = format_directory_path(init_params->output_directory);
	
	/*
	==========================================================================================
	Make Intermediate File Names
	==========================================================================================
	*/
	char** part_filenames = (char**)malloc((params->num_threads) * sizeof(char*));
	int part_num;
	for (part_num = 0; part_num < params->num_threads; part_num++) {
		part_filenames[part_num] = (char*)malloc(MAX_PART_FILE_NAME_LENGTH * sizeof(char));
		part_filenames[part_num] = makePartFilename(part_num, init_params->run_name, 1);
	}

	/*
	==========================================================================================
	Merge part-files
	==========================================================================================
	*/
	printf("Merging part files.\n");
	if (params->reduction_level == 1) {
		char* fname = merge_RC1_Parts(outDir, part_filenames, params);
	}
	else if (params->reduction_level == 2) {
		//char* fname = merge_RC2_Parts(outDir, part_filenames, params->num_threads, header, concat(outDir, (const char*)image_name));
	}
	else if (params->reduction_level == 3) {
		//char* fname = merge_RC3_Parts(outDir, part_filenames, params->num_threads, header, concat(outDir, (const char*)image_name));
	}
	else if (params->reduction_level == 4) {
		char* fname = merge_RC4_Parts(outDir, part_filenames, params, concat(outDir, (const char*)image_name));
	}
	printf("Done.\n");


	for (part_num = 0; part_num < params->num_threads; part_num++) {
		free(part_filenames[part_num]);
	}
	free(params);
	free(part_filenames);

	return 0;
}

int tests (InitParams *init_params) {
	return _sequence_file_read_test();
}

int rc (InitParams *init_params) {
	
	/*
	========================================================================================== 
	Process Inputs
	==========================================================================================
	*/
	char *darkFile = init_params->dark_filename;
	char *imageFile = init_params->image_filename;
    char *outDir = format_directory_path(init_params->output_directory);
    
	// read and validate input parameters
	InputParams *params = (InputParams *)malloc(sizeof(InputParams));
	get_input_params (init_params->params_filename, &params);
	
	// if MRC file, replace params from MRC header
	if (params->source_file_type == 1) {
		MRCHeader *mrc_header = (MRCHeader *)malloc(sizeof(MRCHeader));
        parseMRCHeader (imageFile, &mrc_header);
		compile_missing_params (mrc_header, &params);
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
	const char *image_name = getFilenameFromPath (imageFile);
	const char *dark_name = getFilenameFromPath (darkFile);
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
	
	loadData(params->source_file_type, imageFile, frameBuffer, header->frame_offset, header->nz, b);
	//serializeFrames(frameBuffer, header->nx, header->ny, "image_frame.txt", 1);
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
		
		loadData (params->dark_file_type, darkFile, darkBuffer, 0, header->num_dark_frames, d);
		//serializeFrames(darkBuffer, header->nx, header->ny, "dark_frame.txt", 1);
		printf("Dark data loaded.\n");
		
		getDarkMax (darkBuffer, d, darkFrame, params->num_threads);
		free(darkBuffer);
		
		printf("Dark noise estimation complete.\n");
		
	} else {
		
		DataSize d1 = { header->nx, header->ny, 1, params->source_bit_depth };
		loadData (params->dark_file_type, darkFile, darkFrame, 0, 1, d1);
		//serializeFrames(darkFrame, header->nx, header->ny, "dark_frame.txt", 1);
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
		part_filenames[part_num] = reduceCompress ( 	part_num, 
														(const char*)image_name,
														outDir,
														frameBuffer + chunk_offset, 
														darkFrame, 
														frame_start_indices[part_num], 
														part_b, 
														params,
														header,
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
	if (params->reduction_level == 1) {
		char* fname = merge_RC1_Parts(outDir, part_filenames, params);
	} else if (params->reduction_level == 2) {
		//char* fname = merge_RC2_Parts(outDir, part_filenames, params->num_threads, header, concat(outDir, (const char*)image_name));
	} else if (params->reduction_level == 3) {
		//char* fname = merge_RC3_Parts(outDir, part_filenames, params->num_threads, header, concat(outDir, (const char*)image_name));
	} else if (params->reduction_level == 4) {
		char* fname = merge_RC4_Parts(outDir, part_filenames, params, concat(outDir, (const char*)image_name));
	}
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


int _tmain(int argc, TCHAR * argv[]) {

	InitParams *init_params = (InitParams *)malloc(sizeof(InitParams));
	parse_init_params(argc, argv, &init_params);

	VERBOSITY = init_params->verbosity;

	if (init_params->mode == 0) {
		rc(init_params);
	}
	else if (init_params->mode == 1) {
		de(init_params);
	}
	else if (init_params->mode == 2) {
		merge(init_params);
	}
	else if (init_params->mode == 3) {
		tests(init_params);
	}
	else {
		printf("Failed to initialize. Unknown mode: %d. Acceptable modes are: -rc and -de", init_params->mode);
		fprintf(stderr, "Usage::\n");
		fprintf(stderr, "./recode -rc <input image filename> <dark image filename> <input parameter filename> <output directory>\n");
		fprintf(stderr, "./recode -de <recode filename> <output directory>\n");
		exit(0);
	}
}

/*
extern "C"
{
	__declspec(dllexport) int py_rc(const char *image_filename, const char *dark_filename, const char *params_filename, const char *output_directory, uint8_t verbosity);
	__declspec(dllexport) uint16_t* py_de(const char *recode_filename, const char *output_directory, uint8_t make_sparse, uint8_t verbosity);
}

int py_rc(const char *image_filename, const char *dark_filename, const char *params_filename, const char *output_directory, uint8_t verbosity) {
	printf("py called me.\n");
	printf("Image Filename: %s\n", image_filename);
	printf("Dark Filename: %s\n", dark_filename);
	printf("Params Filename: %s\n", params_filename);
	printf("Output Directory: %s\n", output_directory);
	printf("Verbosity: %d\n", verbosity);

	InitParams *init_params = (InitParams *)malloc(sizeof(InitParams));
	init_params->mode = 0;
	init_params->image_filename = (char*)image_filename;
	init_params->dark_filename = (char*)dark_filename;
	init_params->params_filename = (char*)params_filename;
	init_params->output_directory = (char*)output_directory;
	init_params->verbosity = verbosity;

	return rc(init_params);
}

uint16_t* py_de(const char *recode_filename, const char *output_directory, uint8_t make_sparse, uint8_t verbosity) {
	printf("py called C.de.\n");
	printf("Image Filename: %s\n", recode_filename);
	printf("Output Directory: %s\n", output_directory);
	printf("Sparse: %d\n", make_sparse);
	printf("Verbosity: %d\n", verbosity);

	InitParams *init_params = (InitParams *)malloc(sizeof(InitParams));
	init_params->mode = 1;
	init_params->image_filename = (char*)recode_filename;
	init_params->output_directory = (char*)output_directory;
	init_params->verbosity = verbosity;

	uint16_t N = 1500;
	uint16_t *darkBuffer = (uint16_t *)malloc((2*N + 3) * sizeof(uint16_t));
	darkBuffer[0] = N;
	darkBuffer[1] = 4096;
	darkBuffer[2] = 4096;
	for (uint16_t i = 1; i < N + 1; i++) {
		darkBuffer[i*2+1] = 200 + i;
		darkBuffer[i*2+2] = 200 + i+1;
	}
	return darkBuffer;
}

*/

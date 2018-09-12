
/**************************************************************************************************************
to compile:	gcc recode_comp_alg_comparisons.c -lz -lbz2 -llzma -lsnappy -llz4 -O3 -o recode_comp_alg_comparisons -lm -fopenmp
to run:		export OMP_NUM_THREADS=11
			./recode_comp_alg_comparisons

-lz: 	links zlib
-lbz2: 	links bzip2
-lm: links math.h NOTE: it is important to keep -lm as the last option
-o:	 write build output to file
-O3: optimizer flag
**************************************************************************************************************/

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
#include <bzlib.h>
#include <lzma.h>
#include "snappy-c.h"
#include "lz4.h"

#include <dirent.h>

#include "recodefs.h"
#include "recode_utils.h"
#include "commons.h"
#include "compressor.h"
#include "logger.h"
#include "mrc_header.h"
#include "argparser.h"
#include "recode_header.h"

#include "fileutils.h"

#include "L1_next.h"
#include "L2_next.h"
#include "L3_next.h"
#include "L4.h"



char* getDarkFilename (int i) {

	switch (i) {

		case 0:
			return "/home/abhik/code/DataCompression/simulated_images_2/0.0025/Simulated_Dark_Frame.bin";

		case 1:
			return "/home/abhik/code/DataCompression/simulated_images_2/0.005/Simulated_Dark_Frame.bin";

		case 2:
			return "/home/abhik/code/DataCompression/simulated_images_2/0.01/Simulated_Dark_Frame.bin";

		case 3:
			return "/home/abhik/code/DataCompression/simulated_images_2/0.02/Simulated_Dark_Frame.bin";

		case 4:
			return "/home/abhik/code/DataCompression/simulated_images_2/0.03/Simulated_Dark_Frame.bin";

		case 5:
			return "/home/abhik/code/DataCompression/simulated_images_2/0.04/Simulated_Dark_Frame.bin";

		case 6:
			return "/home/abhik/code/DataCompression/simulated_images_2/0.05/Simulated_Dark_Frame.bin";

		case 7:
			return "/home/abhik/code/DataCompression/simulated_images_2/0.06/Simulated_Dark_Frame.bin";

		case 8:
			return "/home/abhik/code/DataCompression/simulated_images_2/0.07/Simulated_Dark_Frame.bin";

		case 9:
			return "/home/abhik/code/DataCompression/simulated_images_2/0.08/Simulated_Dark_Frame.bin";

		case 10:
			return "/home/abhik/code/DataCompression/simulated_images_2/0.09/Simulated_Dark_Frame.bin";

		case 11:
			return "/home/abhik/code/DataCompression/simulated_images_2/0.1/Simulated_Dark_Frame.bin";

		default:
			return "";

	}

}


char* getImageFilename (int i) {

	switch (i) {

		case 0:
			return "/home/abhik/code/DataCompression/simulated_images_2/0.0025/L0_Simulated.bin";

		case 1:
			return "/home/abhik/code/DataCompression/simulated_images_2/0.005/L0_Simulated.bin";

		case 2:
			return "/home/abhik/code/DataCompression/simulated_images_2/0.01/L0_Simulated.bin";

		case 3:
			return "/home/abhik/code/DataCompression/simulated_images_2/0.02/L0_Simulated.bin";

		case 4:
			return "/home/abhik/code/DataCompression/simulated_images_2/0.03/L0_Simulated.bin";

		case 5:
			return "/home/abhik/code/DataCompression/simulated_images_2/0.04/L0_Simulated.bin";

		case 6:
			return "/home/abhik/code/DataCompression/simulated_images_2/0.05/L0_Simulated.bin";

		case 7:
			return "/home/abhik/code/DataCompression/simulated_images_2/0.06/L0_Simulated.bin";

		case 8:
			return "/home/abhik/code/DataCompression/simulated_images_2/0.07/L0_Simulated.bin";

		case 9:
			return "/home/abhik/code/DataCompression/simulated_images_2/0.08/L0_Simulated.bin";

		case 10:
			return "/home/abhik/code/DataCompression/simulated_images_2/0.09/L0_Simulated.bin";

		case 11:
			return "/home/abhik/code/DataCompression/simulated_images_2/0.1/L0_Simulated.bin";

		default:
			return "";

	}

}

char* getAlgorithmName (int i) {
	
	switch (i) {
		case 0:
			return "DEFLATE";
			
		case 1:
			return "BZIP";
			
		case 2:
			return "LZMA";
			
		case 3:
			return "SNAPPY";
			
		case 4:
			return "LZ4";
			
		default:
			return "";
	}
	
}

void evaluate_compression_algorithm (float dose_rate, int reduction_level, int compression_scheme, int compression_level, char* darkFile, char* imageFile) {
	
	// read and validate input parameters
	InputParams *params = (InputParams *)malloc(sizeof(InputParams));
	get_input_params ("in/recode_params_2.txt", &params);

	params->compression_scheme = compression_scheme;
	params->compression_level = compression_level;
	
	// replace these with argvs later
	//char* darkFile = "/home/abhik/code/DataCompression/simulated_images_2/0.0025/Simulated_Dark_Frame.bin";
	//char* imageFile = "/home/abhik/code/DataCompression/simulated_images_2/0.0025/L0_Simulated.bin";
	
	/*
	========================================================================================== 
	Create ReCoDE header
	==========================================================================================
	*/
	RCHeader *header = (RCHeader *)malloc(sizeof(RCHeader));
	char *image_name = getFilenameFromPath (imageFile);
	char *dark_name = getFilenameFromPath (darkFile);
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
	printf("Image data loaded.\n");
	
	//serializeFrames (frameBuffer, header->nx, header->ny, "Frame_0.txt", 0);
	
	/*
	========================================================================================== 
	Load Dark Data
	==========================================================================================
	*/
	uint32_t n_pixels_in_frame = header->nx * header->ny;
	uint16_t *darkFrame = (uint16_t *)calloc(n_pixels_in_frame, sizeof(uint16_t));
	
	//DataSize d1 = { header->nx, header->ny, 1, params->source_bit_depth };
	//loadData (darkFile, darkFrame, 0, 1, d1);
	printf("Pre-computed dark noise estimates loaded.\n");
		
	
	
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
		printf("RCT %d: n_frames_in_thread = %lu, frame_start_index = %lu\n", part_num, n_frames_in_thread[part_num], frame_start_index);
		frame_start_index += n_frames_in_thread[part_num];
	}
	
	float *run_metrics = (float*) calloc (params->num_threads * 5, sizeof(float));
	
	
	/*
	========================================================================================== 
	Reduce and Compress into part-files
	==========================================================================================
	*/
	uint64_t original_frame_size = 4096*4096*2;
	unsigned long nz;
	
	double omp_start = omp_get_wtime();
	struct timeval start, end;
	gettimeofday(&start, NULL);

	#pragma omp parallel for num_threads(params->num_threads)
	for (part_num = 0; part_num < params->num_threads; part_num++) {
		
		printf("RCT %d: NZ = %u\n", part_num, n_frames_in_thread[part_num]);
		
		DataSize part_b = { header->nx, header->ny, n_frames_in_thread[part_num], params->source_bit_depth };
		
		unsigned long chunk_offset = frame_start_indices[part_num] * part_b.nx * part_b.ny;
		
		printf("RCT %d: frame_start_index = %lu, Chunk Offset = %lu\n", part_num, frame_start_indices[part_num], chunk_offset);

		if (reduction_level == 1) {

			part_filenames[part_num] = reduceCompress_L1 ( 	part_num, 
															image_name, 
															frameBuffer + chunk_offset, 
															darkFrame, 
															frame_start_indices[part_num], 
															part_b, 
															params,
															run_metrics + part_num*5
														);

		} else if (reduction_level == 2) {

			part_filenames[part_num] = reduceCompress_L2 ( 	part_num, 
															image_name, 
															frameBuffer + chunk_offset, 
															darkFrame, 
															frame_start_indices[part_num], 
															part_b, 
															params,
															run_metrics + part_num*5
														);

		} else if (reduction_level == 3) {

			part_filenames[part_num] = reduceCompress_L3 ( 	part_num, 
															image_name, 
															frameBuffer + chunk_offset, 
															darkFrame, 
															frame_start_indices[part_num], 
															part_b, 
															params,
															run_metrics + part_num*5
														);

		} else if (reduction_level == 4) {

			part_filenames[part_num] = reduceCompress_L4 ( 	part_num, 
															image_name, 
															frameBuffer + chunk_offset, 
															darkFrame, 
															frame_start_indices[part_num], 
															part_b, 
															params,
															run_metrics + part_num*5
														);

		}
		
	
		float RSz = run_metrics[part_num*5 + 2]/50;
		float CSz = run_metrics[part_num*5 + 3]/50;

		printf("ReT = %f, CoT = %f, RSz = %f, CSz = %f, DeT = %f, ReR = %f, CoR = %f, RCR = %f \n", 
				run_metrics[part_num*5 + 0]/50, 
				run_metrics[part_num*5 + 1]/50, 
				RSz,
				CSz,
				run_metrics[part_num*5 + 4]/50,
				original_frame_size/RSz,
				RSz/CSz,
				original_frame_size/CSz);
		
		/*
		========================================================================================== 
		Log Results
		========================================================================================== 
		*/
		FILE *logfile2 = fopen("comp_algo_comparisons.txt", "a");
		
		char* compression_scheme_name = getAlgorithmName(compression_scheme);
		fprintf(logfile2, "%4f \t L%d \t %s \t %d \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f\n", 
							dose_rate,
							reduction_level,
							compression_scheme_name,
							compression_level,
							run_metrics[part_num*5 + 0]/50, 
							run_metrics[part_num*5 + 1]/50, 
							RSz,
							CSz,
							run_metrics[part_num*5 + 4]/50,
							original_frame_size/RSz,
							RSz/CSz,
							original_frame_size/CSz);
		fclose(logfile2);


	}
	
	double total_time = omp_get_wtime()-omp_start;
	gettimeofday(&end, NULL);
	double delta = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
	printf("Wall Time: %f, Total OMP Time: %f, Effective Avg. Time per Frame: %f\n", delta, total_time, total_time/(part_num*1.0));

	
	/*
	========================================================================================== 
	Clean-up (ReCo)
	==========================================================================================
	*/
	free(darkFrame);
	free(frameBuffer);
	free(params);
	for (part_num = 0; part_num < params->num_threads; part_num++) {
		free(part_filenames[part_num]);
	}
	free(part_filenames);

	
	return;
}



int main (int argc, char *argv[]) {


	/*
	========================================================================================== 
	Create Log
	==========================================================================================
	*/
	time_t timer;
	char buffer[27];
	struct tm* tm_info;

	time(&timer);
	tm_info = localtime(&timer);
	strftime(buffer, 27, "%d-%b-%Y %H:%M:%S", tm_info);
		
	FILE *logfile = fopen("comp_algo_comparisons.txt", "a");
	
	fprintf(logfile, "\nTest Start Time: %s\n", buffer);
	fprintf(logfile, "=======================================\n");
	fprintf(logfile, "DR \t RL \t CS \t CL \t ReT \t CoT \t RSz \t CSz \t DeT \t ReR \t CoR \t RCR\n");
	fclose(logfile);

	float dose_rates[12] = {0.0025, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1};

    int compression_algorithm = 0;
	int dr, rl;
	for(dr = 0; dr <= 12; dr++) {
		for (rl = 1; rl < 2; rl++) {
			for (compression_algorithm = 0; compression_algorithm < 1; compression_algorithm++) {
				if (compression_algorithm == 0 || compression_algorithm == 1) {
					evaluate_compression_algorithm (dose_rates[dr], rl, compression_algorithm, 1, getDarkFilename(dr), getImageFilename(dr));       // fastest
					//evaluate_compression_algorithm (dose_rates[dr], rl, compression_algorithm, 9, getDarkFilename(dr), getImageFilename(dr));     // optimal compression
				} else {
					//evaluate_compression_algorithm (dose_rates[dr], rl, compression_algorithm, 1, getDarkFilename(dr), getImageFilename(dr));
				}
			}
		}
	}

}




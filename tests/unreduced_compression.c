#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include "zlib.h"
#include <dirent.h>

#include "recodefs.h"
#include "commons.h"
#include "compressor.h"
#include "fileutils.h"
#include "logger.h"
#include "mrc_header.h"
#include "argparser.h"
#include "recode_header.h"


uint64_t compress_serialize_frame (uint8_t process_id, uint16_t *imageBuffer, uint64_t nFrames) {

	const char* filename = makePartFilename(process_id, "unreduced_compressed_data", 1);

	//FILE *fp = fopen(concat("out/",filename), "wb");
	FILE *fp = fopen(concat("../../../../../../scratch/loh/abhik/code_data/DataCompression/on_the_fly_output/",filename), "wb");

	uint32_t n_pixels_in_frame 		= 4096 * 512;
	uint32_t n_bytes_in_frame		= n_pixels_in_frame * 2;
	uint8_t  *compressed_frame 		= (uint8_t*) malloc (n_bytes_in_frame * sizeof(uint8_t));
	uint32_t n_compressed_bytes;
	uint64_t total_compressed_bytes	= 0;
	uint64_t frame_start_index;

	uint64_t frame_no;
	for (frame_no = 0; frame_no < nFrames; frame_no++) {

		n_compressed_bytes			= 0;
		frame_start_index			= frame_no * n_pixels_in_frame;

		//printf("RCT %d: Starting frame %"PRIu64".\n", process_id, frame_no);

		gzip_compress_stream_2 ((uint8_t *)(imageBuffer + frame_start_index), n_pixels_in_frame*2, 1, 0, &n_compressed_bytes, compressed_frame);
		fwrite (compressed_frame, sizeof(uint8_t), n_compressed_bytes, fp);

		total_compressed_bytes += (uint64_t)n_compressed_bytes;

		//printf("RCT %d: Finished frame %"PRIu64".\n", process_id, frame_no);
	}

	fclose(fp);
	free(compressed_frame);
	return total_compressed_bytes;
}


int main(int argc, char *argv[]) {


	uint64_t nImageFrames;
	int replicates = 1;
	
	
	double dose_rates[4] = {0.8, 1.6, 3.2, 6.4};
	
	int dr;
	for (dr = 0; dr < 4; dr++) {

		double dose_rate = dose_rates[dr];

		/*
		========================================================================================== 
		Load Image Data
		==========================================================================================
		*/
		const char* imageFile;
		if (dose_rate == 0.8) {
			imageFile = "/home/abhik/code/DataCompression/beam_blanker_binary_data/12-06-23.598.bin";
			nImageFrames = 700*1;
		} else if (dose_rate == 1.6) {
			imageFile = "/home/abhik/code/DataCompression/beam_blanker_binary_data/12-13-59.436.bin";
			nImageFrames = 700*1;
		} else if (dose_rate == 3.2) {
			//imageFile = "/home/abhik/code/DataCompression/beam_blanker_binary_data/12-12-04.498.bin";
			imageFile = "../../../../../../scratch/loh/abhik/code_data/DataCompression/beam_blanker_binary_data/12-12-04.498.bin";
			nImageFrames = 700*1;
		} else if (dose_rate == 6.4) {
			imageFile = "/home/abhik/code/DataCompression/beam_blanker_binary_data/12-19-00.764.bin";
			nImageFrames = 700*1;
		}


		DataSize b = { 4096, 512, nImageFrames, 16 };
		uint64_t sz_frameBuffer = b.nx * b.ny * b.nz * sizeof(uint16_t);
		uint16_t *frameBuffer = (uint16_t *)malloc(sz_frameBuffer);
		loadData(imageFile, frameBuffer, 0, nImageFrames, b);
		printf("Image data loaded.\n");

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
			
		FILE *logfile = fopen("unreduced_comression_results", "a");
			
		fprintf(logfile, "\nRun Start Time: %s\n", buffer);
		fprintf(logfile, "Data Saved to Network Drive over IPoIB\n");
		fprintf(logfile, "======================================\n");
		fprintf(logfile, "DR \t NT \t n_Fr \t C_Sz \t T \n");
		fclose(logfile);

		
		uint8_t NT;
		for (NT = 5; NT < 51; NT = NT + 5) { 
		
			/*
			========================================================================================== 
			Compute and initialize thread variables
			==========================================================================================
			*/
			uint64_t* n_frames_in_thread 		= (uint64_t*)malloc((NT)*sizeof(uint64_t));
			uint64_t* frame_start_indices 		= (uint64_t*)malloc((NT)*sizeof(uint64_t));
			uint64_t* compressed_sz_per_thread 	= (uint64_t*)malloc((NT)*sizeof(uint64_t));
			
			uint64_t rem = nImageFrames % NT;
			uint64_t min_frames_per_thread = floor(((nImageFrames)*1.0)/((NT)*1.0));
			uint64_t frame_start_index = 0;
			
			uint8_t part_num;
			for (part_num = 0; part_num < NT; part_num++) {
				
				if (part_num < rem) {
					n_frames_in_thread[part_num] = min_frames_per_thread + 1;
				} else {
					n_frames_in_thread[part_num] = min_frames_per_thread;
				}
				
				frame_start_indices[part_num] = frame_start_index;
				printf("RCT %d: n_frames_in_thread = %lu, frame_start_index = %lu\n", part_num, n_frames_in_thread[part_num], frame_start_index);
				frame_start_index += n_frames_in_thread[part_num];
			}

			/*
			========================================================================================== 
			Compress
			==========================================================================================
			*/
			struct timeval start, end;
			gettimeofday(&start, NULL);
			double omp_start = omp_get_wtime();

			#pragma omp parallel for num_threads(NT)
			for (part_num = 0; part_num < NT; part_num++) {
				compressed_sz_per_thread[part_num] = 
						compress_serialize_frame (part_num, frameBuffer + frame_start_indices[part_num], n_frames_in_thread[part_num]);
			}

			double total_time = omp_get_wtime()-omp_start;
			gettimeofday(&end, NULL);
			double delta = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
			printf("Wall Time: %f, Total OMP Time: %f, Effective Avg. Time per Frame: %f\n", delta, total_time, total_time/(part_num*1.0));

			/*
			========================================================================================== 
			Log Results
			========================================================================================== 
			*/
			uint64_t total_sz = 0;
			for (part_num = 0; part_num < NT; part_num++) {
				total_sz += compressed_sz_per_thread[part_num];
			}

			FILE *logfile2 = fopen("unreduced_comression_results", "a");
			
			fprintf(logfile2, "%.1f \t %d \t %"PRIu64" \t %"PRIu64" \t %.4f\n", 
								dose_rate, NT, nImageFrames, total_sz, total_time);
			fclose(logfile2);

		}

		/*
		========================================================================================== 
		Clean-up
		==========================================================================================
		*/
		free(frameBuffer);
		

	}

}